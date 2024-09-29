import datetime
import json
from io import BufferedReader
from typing import Union, Tuple

from celery import shared_task

from configs import dify_config
from configs.websocket_config import websocket_handler
from extensions.ext_database import db
from extensions.ext_storage import storage
from models.account import Account
from models.model import EndUser, UploadFile
from models.molecular_docking import MolecularDockingTask, Status
import requests
import logging
import click
import requests

from services.molecular_docking.tool_center_position_service import ToolCenterPositionService

ALLOWED_EXTENSIONS = ["pdb", "sdf", "mol"]

if dify_config.MOLECULAR_DOCKING_API_URL == "":
    logging.error(
        click.style(
            "分子对接API URL不能为空, 请在配置文件.env中设置MOLECULAR_DOCKING_API_URL, 否则分子对接功能无法正常使用",
            fg='red', bold=True))


class MolecularDockingService:
    @classmethod
    def start_task(cls, task_name: str, pdb_file_id: str, center_x: float, center_y: float, center_z: float,
                   size_x: float, size_y: float, size_z: float, ligand_file_ids: list[str], out_pose_num: int,
                   user: Union[Account, EndUser], start_celery: bool = True) -> MolecularDockingTask:
        """
        Start molecular docking task.
        :param task_name: task name
        :param pdb_file_id: pdb file id
        :param center_x: center x
        :param center_y: center y
        :param center_z: center z
        :param size_x: size x
        :param size_y: size y
        :param size_z: size z
        :param ligand_file_ids: ligand file ids
        :param out_pose_num: output pose number
        :param user: user
        :param start_celery: 是否启动消息队列进行任务处理，true表示使用消息队列一个一个处理任务，false表示不启动直接调用api接口
        :return: molecular docking task object
        """
        status = Status.PROCESSING.status
        if start_celery:
            status = Status.PENDING.status
            logging.info(click.style(f"分子对接使用消息队列进行任务处理", fg='blue', bold=True))
        # 初始化任务。并且默认任务状态是正在处理中，保存到数据表中
        molecular_docking_task = MolecularDockingTask(
            task_name=task_name,
            pdb_file_id=pdb_file_id,
            center_x=center_x,
            center_y=center_y,
            center_z=center_z,
            size_x=size_x,
            size_y=size_y,
            size_z=size_z,
            ligand_file_ids=ligand_file_ids,
            out_pose_num=out_pose_num,
            status=status,
            created_by=user.id,
        )
        db.session.add(molecular_docking_task)
        db.session.commit()

        if start_celery:
            # 启动消息队列进行任务处理
            molecular_docking_celery_task.apply_async(
                args=[user.serialize, center_x, center_y, center_z,
                      size_x, size_y, size_z, out_pose_num, molecular_docking_task.serialize],
            )
            return molecular_docking_task
        else:
            pdb_file_buffer = cls.get_upload_file_buffer(pdb_file_id, user)
            ligand_file_buffer_list = cls.get_upload_file_buffer(ligand_file_ids, user, return_list=True)
            return cls.main_processor(user, task_name, pdb_file_buffer, ligand_file_buffer_list, center_x, center_y,
                                      center_z,
                                      size_x, size_y, size_z, out_pose_num, molecular_docking_task)

    @classmethod
    def get_center_position(cls, pdb_file_id: str, user: Union[Account, EndUser]) -> dict:
        """
        获取pdb文件中心点坐标
        :param pdb_file_id: pdb文件id
        :param user: 用户
        :return: 坐标
        """
        upload_file = UploadFile.query.filter_by(id=pdb_file_id, created_by=user.id).first()
        file_bytes = storage.load_once(upload_file.key)
        return ToolCenterPositionService.get_box_center(file_bytes)

    @classmethod
    def get_upload_file_buffer(cls, upload_file_ids: Union[str, list[str]], user: Union[Account, EndUser],
                               return_list: bool = False) -> Union[BufferedReader, list[BufferedReader]]:
        """
        Get upload file buffer.
        :param upload_file_ids: upload file ids
        :param user: user
        :param return_list: 是否返回数组形式
        :return: pdb file buffer
        """
        if return_list:
            upload_file_buffer_list = []
            for upload_file_id in upload_file_ids:
                upload_file = UploadFile.query.filter_by(id=upload_file_id, created_by=user.id).first()
                upload_file_buffer_list.append(storage.load_buffer(upload_file.key))
            return upload_file_buffer_list
        else:
            upload_file = UploadFile.query.filter_by(id=upload_file_ids, created_by=user.id).first()
            return storage.load_buffer(upload_file.key)

    @classmethod
    def main_processor(cls, user: Union[Account, EndUser], task_name: str, pdb_file_buffer: BufferedReader,
                       ligand_file_buffer_list: list[BufferedReader], center_x: float, center_y: float, center_z: float,
                       size_x: float, size_y: float, size_z: float, out_pose_num: int,
                       molecular_docking_task: MolecularDockingTask) -> MolecularDockingTask:
        # 调用DockingProcessorAPI进行docking
        try:
            logging.info(click.style(f"{user.name} 开始分子对接任务：{task_name}", fg='blue'))
            result, success_flag = cls.docking_processor(pdb_file_buffer, ligand_file_buffer_list, center_x, center_y,
                                                         center_z,
                                                         size_x, size_y,
                                                         size_z, out_pose_num)
            if success_flag:
                logging.info(click.style(f"分子对接成功", fg='blue'))
                status = Status.SUCCESS.status
                result = json.dumps(result)
            else:
                logging.info(click.style(f"分子对接失败", fg='red', bold=True))
                status = Status.FAILURE.status
        except Exception as e:
            logging.info(click.style(f"分子对接任务失败：{e}", fg='red', bold=True))
            status = Status.FAILURE.status
            result = e

        # 更新任务状态和结果
        molecular_docking_task = MolecularDockingTask.query.filter_by(id=molecular_docking_task.id,
                                                                      created_by=user.id).first()
        molecular_docking_task.status = status
        molecular_docking_task.result = result
        molecular_docking_task.updated_at = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        db.session.commit()

        # todo 使用websocket通知前端任务状态变化为成功或失败
        data = {
            "user_id": user.id,
            "message": json.dumps({
                "task_id": molecular_docking_task.id,
                "status": molecular_docking_task.status
            }),
        }
        requests.post(f"http://127.0.0.1:5001/console/api/molecular-docking/send-websocket-message", data=data)

        return molecular_docking_task

    @classmethod
    def docking_processor(cls, pdb_file_buffer, ligand_file_buffer_list, center_x, center_y, center_z, size_x, size_y,
                          size_z, out_pose_num) -> Tuple[Union[list, str], bool]:
        """
        Docking processor.
        :param pdb_file_buffer: pdb file buffer
        :param ligand_file_buffer_list: ligand file buffer list
        :param center_x: center x
        :param center_y: center y
        :param center_z: center z
        :param size_x: size x
        :param size_y: size y
        :param size_z: size z
        :param out_pose_num: output pose number
        :return: result data and success or not
        """
        docking_params = {
            "center_x": center_x,
            "center_y": center_y,
            "center_z": center_z,
            "size_x": size_x,
            "size_y": size_y,
            "size_z": size_z,
            "num_modes": out_pose_num,
            "out": "out.sdf"
        }
        files = {
            'receptor_file': pdb_file_buffer,
            # todo 这里后期可能需要支持多个ligand文件，现在默认只支持一个
            'ligand_file': ligand_file_buffer_list[0] if len(ligand_file_buffer_list) == 1 else None
        }

        response = requests.post(dify_config.MOLECULAR_DOCKING_API_URL, files=files, data=docking_params)
        result_data = response.json()

        if 'error' in result_data:
            return result_data['error'], False

        return result_data['output'], True


# acks_late 设置为 True 时，任务的消息确认（acknowledgement）会在任务执行完成后才发送，确保任务在失败或 worker 崩溃时能重新被执行。
# time_limit 设置为 120 秒，任务的执行时间不能超过 120 秒，超过这个时间，任务会被自动取消。
# bind 如果设置为 True，任务将绑定到当前任务实例（self），从而允许你在任务中访问 self（即任务对象本身）。这在需要访问任务元数据（例如任务ID、重试次数）时非常有用
@shared_task(queue='molecular_docking', bind=True, time_limit=120, acks_late=True)
def molecular_docking_celery_task(self, user_dict: dict, center_x: float, center_y: float,
                                  center_z: float,
                                  size_x: float, size_y: float, size_z: float, out_pose_num: int,
                                  molecular_docking_task_dict: dict):
    logging.info(click.style(f"molecular_docking_celery_task 开始执行，任务ID：{self.request.id}", fg='blue'))
    # 因为celery需要序列化之后才能传递到该参数，所以现在的user和molecular_docking_task都是json类型的，需要进行反序列化
    molecular_docking_task = MolecularDockingTask(**molecular_docking_task_dict)
    user = Account(**user_dict)
    # 更新任务状态为正在处理中
    MolecularDockingTask.query.filter_by(id=molecular_docking_task.id, created_by=user.id).update(
        {'status': Status.PROCESSING.status})
    db.session.commit()

    data = {
        "user_id": user.id,
        "message": json.dumps({
            "task_id": molecular_docking_task.id,
            "status": molecular_docking_task.status
        }),
    }
    requests.post(f"http://127.0.0.1:5001/console/api/molecular-docking/send-websocket-message", data=data)


    # 获取pdb和ligand文件buffer
    pdb_file_buffer = MolecularDockingService.get_upload_file_buffer(molecular_docking_task.pdb_file_id, user)
    ligand_file_buffer_list = MolecularDockingService.get_upload_file_buffer(
        molecular_docking_task.ligand_file_ids, user, return_list=True)
    MolecularDockingService.main_processor(
        user,
        molecular_docking_task.task_name,
        pdb_file_buffer,
        ligand_file_buffer_list,
        center_x,
        center_y,
        center_z,
        size_x,
        size_y,
        size_z,
        out_pose_num,
        molecular_docking_task
    )
    logging.info(click.style(f"molecular_docking_celery_task 执行完成，任务ID：{self.request.id}", fg='blue'))
