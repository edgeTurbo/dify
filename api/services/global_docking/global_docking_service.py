"""
@Time    : 2024/10/16 下午1:36
@Author  : bigboss
@Description : 全局分子对接的服务类
"""
import io
import json
import logging
import os.path
from io import BufferedReader
from typing import Union, Tuple

import click
import requests
from celery import shared_task
from rdkit import Chem

from configs import dify_config
from controllers.inner_api.websocket.websocket import calling_websocket_internal_send
from extensions.ext_database import db
from extensions.ext_storage import storage
from models.account import Account
from models.model import EndUser, UploadFile
from models.molecular_docking import Status
from models.sciminer import SciminerHistoryTask
from models.sciminer_models.global_docking import GlobalDockingTask

if dify_config.GLOBAL_DOCKING_API_URL == "" or dify_config.GLOBAL_DOCKING_API_URL is None:
    logging.error(
        click.style(
            "全局对接API URL不能为空, 请在配置文件.env中设置MOLECULAR_DOCKING_API_URL, 否则分子对接功能无法正常使用",
            fg='red', bold=True))

if dify_config.INNER_API is None or dify_config.INNER_API_KEY is None:
    logging.error(
        click.style(
            "内网API地址或API KEY未设置, 请在配置文件.env中设置INNER_API为true和INNER_API_KEY, 否则分子对接功能无法正常使用，如websocket",
            fg='red', bold=True)
    )


class GlobalDockingService:
    @classmethod
    def start_task(cls, task_name: str, fasta_file_id: str, ligand_file_ids: list[str], out_pose_num: int,
                   user: Union[Account, EndUser],
                   start_celery: bool = True) -> GlobalDockingTask:
        """
        Start molecular docking task.
        :param task_name: task name
        :param fasta_file_id: fasta file id
        :param ligand_file_ids: ligand file id list
        :param out_pose_num: output pose number
        :param user: user
        :param start_celery: 是否启动消息队列进行任务处理，true表示使用消息队列一个一个处理任务，false表示不启动直接调用api接口
        :return: global docking task object
        """
        status = Status.PROCESSING.status
        if start_celery:
            status = Status.PENDING.status
            logging.info(click.style(f"全局对接使用消息队列进行任务处理", fg='green', bold=True))
        # 初始化任务，保存到数据表中
        global_docking_task = GlobalDockingTask(
            fasta_file_id=fasta_file_id,
            ligand_file_ids=ligand_file_ids,
            out_pose_num=out_pose_num,
            created_by=user.id,
            status=status,
        )
        db.session.add(global_docking_task)
        db.session.commit()

        # 新增任务历史记录,记录任务id和任务名称
        sciminer_history_task = SciminerHistoryTask(
            task_id=global_docking_task.id,
            task_name=task_name,
            task_type="GLOBAL_DOCKING",
            label="Global docking",
            status=status,
            created_by=user.id,
        )
        db.session.add(sciminer_history_task)
        db.session.commit()

        if start_celery:
            # 启动消息队列进行任务处理
            global_docking_celery_task.apply_async(
                args=[user.serialize, out_pose_num, global_docking_task.serialize],
            )
            return global_docking_task
        else:
            fasta_file_buffer = cls.get_upload_file_buffer(fasta_file_id, user)
            ligand_file_buffer_list = cls.get_upload_file_buffer(ligand_file_ids, user, return_list=True)
            return cls.main_processor(user, task_name, fasta_file_buffer, ligand_file_buffer_list, out_pose_num,
                                      global_docking_task, start_celery)

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
                upload_file_buffer_list.append(storage.load_buffer(upload_file.key, upload_file.name))
            return upload_file_buffer_list
        else:
            upload_file = UploadFile.query.filter_by(id=upload_file_ids, created_by=user.id).first()
            return storage.load_buffer(upload_file.key, upload_file.name)

    @classmethod
    def main_processor(cls, user: Union[Account, EndUser], task_name: str, fasta_file_buffer: BufferedReader,
                       ligand_file_buffer_list: list[BufferedReader], out_pose_num: int,
                       global_docking_task: GlobalDockingTask, start_celery: bool = True) -> GlobalDockingTask:
        # 调用DockingProcessorAPI进行docking
        try:
            logging.info(click.style(f"{user.name} 开始全局分子对接任务：{task_name}", fg='green'))
            result, success_flag = cls.docking_processor(fasta_file_buffer, ligand_file_buffer_list, out_pose_num)
            if success_flag:
                logging.info(click.style(f"全局分子对接成功", fg='green'))
                status = Status.SUCCESS.status
                result = json.dumps(result)
            else:
                logging.info(click.style(f"全局分子对接失败", fg='red', bold=True))
                status = Status.FAILURE.status
        except Exception as e:
            import traceback
            traceback.print_exc()
            logging.info(click.style(f"全局分子对接任务失败：{e}", fg='red', bold=True))
            status = Status.FAILURE.status
            result = e
        finally:
            # 关闭所有文件流
            fasta_file_buffer.close()
            for ligand_file_buffer in ligand_file_buffer_list:
                ligand_file_buffer.close()

        # 更新global_docking_task数据表中的结果信息
        global_docking_task = GlobalDockingTask.query.filter_by(id=global_docking_task.id,
                                                                created_by=user.id).first()
        global_docking_task.result = result
        global_docking_task.status = status
        db.session.commit()

        # 在sciminer_history_task数据表中更新任务状态
        SciminerHistoryTask.query.filter_by(task_id=global_docking_task.id, created_by=user.id).update(
            {'status': status}
        )
        db.session.commit()

        if start_celery:
            # 当启动消息队列的时候才进行发送websocket消息
            calling_websocket_internal_send(channel='global_docking', user_id=user.id, message={
                "task_id": global_docking_task.id,
                "status": status
            })

        return global_docking_task

    @classmethod
    def docking_processor(cls, fasta_file_buffer, ligand_file_buffer_list, out_pose_num) -> Tuple[
        Union[list, str], bool]:
        """
        Docking processor.
        :param fasta_file_buffer: fasta file buffer
        :param ligand_file_buffer_list: ligand file buffer list
        :param out_pose_num: output pose number
        :return: result data and success or not
        """

        protein_content = fasta_file_buffer.read().decode('utf-8')
        # 如果开头是>，则去掉>所在的那一行
        protein_lines = protein_content.splitlines()
        # 使用列表推导式来过滤掉以">"开头的行
        filtered_protein_lines = [line for line in protein_lines if not line.startswith(">")]
        protein_content = "\n".join(filtered_protein_lines)

        smiles_list = []

        for ligand_file_buffer in ligand_file_buffer_list:
            mol = Chem.MolFromMolBlock(ligand_file_buffer.read().decode('utf-8'))
            smiles = Chem.MolToSmiles(mol)
            smiles_list.append(smiles)

        docking_params = {
            "protein_content": protein_content,
            "smiles_list": smiles_list,
            # "num_modes": out_pose_num,
        }

        response = requests.post(dify_config.GLOBAL_DOCKING_API_URL, json=docking_params, timeout=600)
        result_data = response.json()

        if 'error' in result_data:
            return result_data['error'], False

        return result_data['message'], True

    # @classmethod
    # def download_task_result(cls, task_id, _range, current_user):
    #     """
    #     下载分子对接任务结果
    #     :param task_id: 分子对接任务id
    #     :param _range: 下载文件范围
    #     :param current_user: 当前用户
    #     :return:
    #     """
    #     molecular_docking_task = MolecularDockingTask.query.filter_by(id=task_id, created_by=current_user.id).first()
    #     if molecular_docking_task is None:
    #         return None
    #     if molecular_docking_task.status == Status.SUCCESS.status:
    #         if _range == 'all':
    #             # 下载全部结果
    #             # 解析json数据，将mol提取出来，写入到sdf文件
    #             result_list = json.loads(molecular_docking_task.result)
    #             sdf_content = ""
    #             for result in result_list:
    #                 sdf_content += result['mol'] + "$$$$\n"
    #             return sdf_content
    #         else:
    #             try:
    #                 range_list = [int(x) for x in _range.split(',')]
    #                 # 下载指定结果
    #                 result_list = json.loads(molecular_docking_task.result)
    #                 filter_result_list = list(map(lambda i: result_list[i], range_list))
    #                 sdf_content = ""
    #                 for result in filter_result_list:
    #                     sdf_content += result['mol'] + "$$$$\n"
    #                 return sdf_content
    #             except Exception as e:
    #                 logging.error(f"下载指定结果失败：{e}")
    #                 return None
    #     else:
    #         return None
    #
    # # todo 这是临时性的方法，后面可能还需要再深究
    # @classmethod
    # def start_task_for_custom_tool(cls, task_name: str, pdb_file_url: str, center_x: float, center_y: float,
    #                                center_z: float,
    #                                size_x: float, size_y: float, size_z: float, ligand_file_urls: list[str],
    #                                out_pose_num: int,
    #                                ) -> dict:
    #     """
    #     Start molecular docking task.
    #     :param task_name: task name
    #     :param pdb_file_url: pdb file url
    #     :param center_x: center x
    #     :param center_y: center y
    #     :param center_z: center z
    #     :param size_x: size x
    #     :param size_y: size y
    #     :param size_z: size z
    #     :param ligand_file_urls: ligand file urls
    #     :param out_pose_num: output pose number
    #     :return: molecular docking task object
    #     """
    #     pdb_file_buffer = cls.download_file(pdb_file_url)
    #     ligand_file_buffer_list: list[BufferedReader] = []
    #     for ligand_file_url in ligand_file_urls:
    #         ligand_file_buffer_list.append(cls.download_file(ligand_file_url))
    #     return cls.main_processor_for_custom_tool(task_name, pdb_file_buffer, ligand_file_buffer_list, center_x,
    #                                               center_y,
    #                                               center_z,
    #                                               size_x, size_y, size_z, out_pose_num)
    #
    # # todo 这是临时性的方法，后面可能还需要再深究
    # @classmethod
    # def main_processor_for_custom_tool(cls, task_name: str, pdb_file_buffer: BufferedReader,
    #                                    ligand_file_buffer_list: list[BufferedReader], center_x: float, center_y: float,
    #                                    center_z: float,
    #                                    size_x: float, size_y: float, size_z: float, out_pose_num: int,
    #                                    ) -> dict:
    #     # 调用DockingProcessorAPI进行docking
    #     try:
    #         logging.info(click.style(f"自定义工具 开始分子对接任务：{task_name}", fg='blue'))
    #         result, success_flag = cls.docking_processor(pdb_file_buffer, ligand_file_buffer_list, center_x, center_y,
    #                                                      center_z,
    #                                                      size_x, size_y,
    #                                                      size_z, out_pose_num)
    #         if success_flag:
    #             logging.info(click.style(f"分子对接成功", fg='blue'))
    #         else:
    #             logging.info(click.style(f"分子对接失败", fg='red', bold=True))
    #     except Exception as e:
    #         logging.info(click.style(f"分子对接任务失败：{e}", fg='red', bold=True))
    #         result = {"error": str(e)}
    #
    #     return result

    @classmethod
    def download_file(cls, file_url):
        try:
            response = requests.get(file_url)
            byte_io = io.BytesIO(response.content)
            byte_io.name = os.path.basename(os.path.basename(file_url))
            buffered_reader = io.BufferedReader(byte_io)
            return buffered_reader
        except Exception as e:
            import traceback
            traceback.print_exc()
            raise ValueError("Download file from url failed: " + str(file_url))


# acks_late 设置为 True 时，任务的消息确认（acknowledgement）会在任务执行完成后才发送，确保任务在失败或 worker 崩溃时能重新被执行。
# time_limit 设置为 120 秒，任务的执行时间不能超过 120 秒，超过这个时间，任务会被自动取消。
# bind 如果设置为 True，任务将绑定到当前任务实例（self），从而允许你在任务中访问 self（即任务对象本身）。这在需要访问任务元数据（例如任务ID、重试次数）时非常有用
@shared_task(queue='global_docking', bind=True, time_limit=120, acks_late=True)
def global_docking_celery_task(self, user_dict: dict, out_pose_num: int,
                               global_docking_task_dict: dict):
    logging.info(click.style(f"global_docking_celery_task 开始执行，任务ID：{self.request.id}", fg='green'))
    # 因为celery需要序列化之后才能传递到该参数，所以现在的user和global_docking_task都是json类型的，需要进行反序列化
    global_docking_task = GlobalDockingTask(**global_docking_task_dict)
    user = Account(**user_dict)
    # 在sciminer_history_task数据表中更新任务状态为处理中
    sciminer_history_task = SciminerHistoryTask.query.filter_by(task_id=global_docking_task.id, user_id=user.id).first()
    if sciminer_history_task is not None and isinstance(sciminer_history_task, SciminerHistoryTask):
        sciminer_history_task.status = Status.PROCESSING.status
        db.session.commit()
    # 设置global_docking_task的状态为处理中
    global_docking_task.status = Status.PROCESSING.status
    calling_websocket_internal_send(channel='global_docking', user_id=user.id, message={
        "task_id": sciminer_history_task.task_id,
        "status": sciminer_history_task.status
    })

    # 获取fasta文件和ligand文件buffer
    fasta_file_buffer = GlobalDockingService.get_upload_file_buffer(global_docking_task.fasta_file_id, user)
    ligand_file_buffer_list = GlobalDockingService.get_upload_file_buffer(
        global_docking_task.ligand_file_ids, user, return_list=True)
    GlobalDockingService.main_processor(
        user,
        sciminer_history_task.task_name,
        fasta_file_buffer,
        ligand_file_buffer_list,
        out_pose_num,
        global_docking_task
    )
    logging.info(click.style(f"global_docking_celery_task 执行完成，任务ID：{self.request.id}", fg='green'))
