import csv
import datetime
import hashlib
import io
import json
import logging
import os.path
import uuid
import zipfile
from io import BufferedReader
from typing import Union, Tuple

import click
import requests
from celery import shared_task

from configs import dify_config
from controllers.inner_api.websocket.websocket import calling_websocket_internal_send
from core.tools.utils.upload_file_utils import UploadFileUtils
from extensions.ext_database import db
from extensions.ext_storage import storage
from models.account import Account
from models.model import EndUser, UploadFile
from models.sciminer_models.molecular_docking import MolecularDockingTask, Status
from models.sciminer_models.sciminer import SciminerHistoryTask
from services.sciminer_services.molecular_docking.tool_center_position_service import ToolCenterPositionService
from services.sciminer_services.molecular_docking.tool_delete_special_ligand_service import ToolDeleteSpecialLigandService
from services.sciminer_services.molecular_docking.tool_ligand_info_service import ToolLigandInfoService
from services.sciminer_services.sciminer_base_service import SciminerBaseService

ALLOWED_EXTENSIONS = ["pdb", "sdf", "mol"]

if dify_config.MOLECULAR_DOCKING_API_URL == "" or dify_config.MOLECULAR_DOCKING_API_URL is None:
    logging.error(
        click.style(
            "口袋对接API URL不能为空, 请在配置文件.env中设置MOLECULAR_DOCKING_API_URL, 否则分子对接功能无法正常使用",
            fg='red', bold=True))

if dify_config.INNER_API is None or dify_config.INNER_API_KEY is None:
    logging.error(
        click.style(
            "内网API地址或API KEY未设置, 请在配置文件.env中设置INNER_API为true和INNER_API_KEY, 否则分子对接功能无法正常使用，如websocket",
            fg='red', bold=True)
    )


class MolecularDockingService(SciminerBaseService):
    service_type = "POCKET_DOCKING"
    task_label = "Pocket docking"

    @classmethod
    def start_task(cls, task_name: str, pdb_file_id: str, center_x: float, center_y: float, center_z: float,
                   size_x: float, size_y: float, size_z: float, ligand_file_ids: list[str], out_pose_num: int,
                   chain: str, residue_number: int, user: Union[Account, EndUser],
                   start_celery: bool = True) -> MolecularDockingTask:
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
        :param chain: pdb chain
        :param residue_number: pdb residue number
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

        # 新增任务历史记录
        sciminer_history_task = SciminerHistoryTask(
            task_id=molecular_docking_task.id,
            task_name=task_name,
            task_type=cls.service_type,
            label=cls.task_label,
            status=status,
            created_by=user.id,
        )
        db.session.add(sciminer_history_task)
        db.session.commit()

        if start_celery:
            # 启动消息队列进行任务处理
            molecular_docking_celery_task.apply_async(
                args=[user.serialize, center_x, center_y, center_z,
                      size_x, size_y, size_z, out_pose_num, chain, residue_number, molecular_docking_task.serialize],
            )
            return molecular_docking_task
        else:
            pdb_file_buffer = UploadFileUtils.get_upload_file_buffer(pdb_file_id, user)
            ligand_file_buffer_list = UploadFileUtils.get_upload_file_buffer(ligand_file_ids, user, return_list=True)
            return cls.main_processor(user, task_name, pdb_file_buffer, ligand_file_buffer_list, center_x, center_y,
                                      center_z,
                                      size_x, size_y, size_z, out_pose_num, chain, residue_number,
                                      molecular_docking_task, start_celery)

    @classmethod
    def get_center_position_and_residue_number_and_chain(cls, pdb_file_id: str, user: Union[Account, EndUser]) -> dict:
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
    def main_processor(cls, user: Union[Account, EndUser], task_name: str, pdb_file_buffer: BufferedReader,
                       ligand_file_buffer_list: list[BufferedReader], center_x: float, center_y: float, center_z: float,
                       size_x: float, size_y: float, size_z: float, out_pose_num: int,
                       chain: str, residue_number: int,
                       molecular_docking_task: MolecularDockingTask, start_celery: bool = True) -> MolecularDockingTask:
        remove_ligand_file_id = None
        remove_ligand_file = None
        # 调用DockingProcessorAPI进行docking
        try:
            # 在进行分子对接任务之前，先通过前端提交的chain和residue_number判断是不是包含在全部配体信息里面，如果包含在配体信息里，则需要将其删除
            # 1. 获取全部配体信息
            if chain == "" or residue_number == "":
                logging.info(click.style(f"前端未提交chain和residue_number，不进行配体信息的删除", fg='blue'))
            else:
                residue_number = int(residue_number)
                logging.info(click.style(f"前端提交的chain和residue_number：{chain}和{residue_number}", fg='blue'))
                unique_ligands = ToolLigandInfoService.extract_unique_ligands(pdb_file_buffer)
                logging.info(click.style(f"全部配体信息：{unique_ligands}", fg='blue'))
                # 2. 判断是否包含在配体信息里
                if (chain, residue_number) in unique_ligands:
                    logging.info(click.style(f"配体信息中包含{chain}和{residue_number}，将其从pdb文件中删除", fg='blue'))
                    # 3. 在配体信息中，进行部分删除
                    pdb_file_buffer = ToolDeleteSpecialLigandService.remove_ligand(pdb_file_buffer, chain, residue_number)
                    # 4. 将新的pdb文件进行保存
                    new_pdb_file_file_key, current_tenant_id = cls.get_pocket_docking_file_path(
                        file_name=str(uuid.uuid4()) + "." + pdb_file_buffer.name.split('.')[-1],
                        user=user
                    )
                    file_content = pdb_file_buffer.read()
                    file_size = len(file_content)
                    pdb_file_buffer.seek(0)
                    storage.save(new_pdb_file_file_key, pdb_file_buffer.read())
                    pdb_file_buffer.seek(0)
                    # 5. 在文件数据表中新增一条修改过后的pdb文件记录
                    new_pdb_file = UploadFile(
                        tenant_id=current_tenant_id,
                        storage_type=dify_config.STORAGE_TYPE,
                        key=new_pdb_file_file_key,
                        name=f"{pdb_file_buffer.name.split('.')[0]}_modified.{pdb_file_buffer.name.split('.')[-1]}",
                        size=file_size,
                        extension=pdb_file_buffer.name.split('.')[-1],
                        mime_type="application/octet-stream",
                        created_by_role=("account" if isinstance(user, Account) else "end_user"),
                        created_by=user.id,
                        created_at=datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None),
                        used=False,
                        hash=hashlib.sha3_256(file_content).hexdigest(),
                    )
                    db.session.add(new_pdb_file)
                    db.session.commit()
                    # 6. 更新任务中的remove_ligand_file_id
                    remove_ligand_file_id = new_pdb_file.id
                    remove_ligand_file = new_pdb_file
                else:
                    logging.info(click.style(f"{chain}和{residue_number}不包含在配体信息中，不做任何处理", fg='blue'))
                    # 还原文件坐标
                    pdb_file_buffer.seek(0)

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
            import traceback
            traceback.print_exc()
            logging.info(click.style(f"分子对接任务失败：{e}", fg='red', bold=True))
            status = Status.FAILURE.status
            result = e
        finally:
            # 关闭所有文件流
            pdb_file_buffer.close()
            for ligand_file_buffer in ligand_file_buffer_list:
                ligand_file_buffer.close()

        # todo 更新任务状态和结果,后期这个MolecularDockingTask的任务状态需要删除，使用sciminer_history_task数据表来接管，现在先全部数据表更新
        molecular_docking_task = MolecularDockingTask.query.filter_by(id=molecular_docking_task.id,
                                                                      created_by=user.id).first()
        molecular_docking_task.status = status
        molecular_docking_task.result = result
        molecular_docking_task.updated_at = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if remove_ligand_file_id is not None:
            molecular_docking_task.remove_ligand_file_id = remove_ligand_file_id
            molecular_docking_task.remove_ligand_file = remove_ligand_file
        db.session.commit()

        # 在sciminer_history_task数据表中更新任务状态
        SciminerHistoryTask.query.filter_by(task_id=molecular_docking_task.id, created_by=user.id).update(
            {'status': status}
        )
        db.session.commit()

        if start_celery:
            # 当启动消息队列的时候才进行发送websocket消息
            calling_websocket_internal_send(channel='molecular_docking', user_id=user.id, message={
                "id": molecular_docking_task.id,
                "task_name": molecular_docking_task.task_name,
                "result": molecular_docking_task.result,
                "status": molecular_docking_task.status,
                "remove_ligand_file": molecular_docking_task.remove_ligand_file.serialize
            })

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

        try:
            response = requests.post(dify_config.MOLECULAR_DOCKING_API_URL, files=files, data=docking_params, timeout=120)
            result_data = response.json()

            if 'error' in result_data:
                return result_data['error'], False
            if 'message' in result_data:
                return result_data['message'], False
            return result_data['output'], True
        except Exception as e:
            error_message = str(e)
            logging.error(
                click.style(
                    error_message,
                    fg='red', bold=True)
            )
            return error_message, False

    @classmethod
    def download_task_result(cls, task_id, _range, current_user, zip_csv_file: bool = True):
        """
        下载分子对接任务结果
        :param task_id: 分子对接任务id
        :param _range: 下载文件范围
        :param current_user: 当前用户
        :param zip_csv_file: 是否要压缩csv文件，默认压缩
        :return:
        """
        molecular_docking_task = MolecularDockingTask.query.filter_by(id=task_id, created_by=current_user.id).first()
        if molecular_docking_task is None:
            return None
        if molecular_docking_task.status == Status.SUCCESS.status:
            zip_buffer = io.BytesIO()
            if _range == 'all':
                with zipfile.ZipFile(zip_buffer, 'w') as _zip:
                    # 下载全部结果
                    # 解析json数据，将mol提取出来，写入到sdf文件
                    result_list = json.loads(molecular_docking_task.result)
                    sdf_content = ""
                    for result in result_list:
                        sdf_content += result['mol'] + "$$$$\n"
                    _zip.writestr(f"{task_id}.sdf", sdf_content)

                    if zip_csv_file:
                        csv_content = cls.get_csv_data(data=result_list)
                        _zip.writestr('result.csv', csv_content)
                zip_buffer.seek(0)
                return zip_buffer
            else:
                try:
                    range_list = [x for x in _range.split(',')]
                    # 下载指定结果
                    result_list = json.loads(molecular_docking_task.result)

                    filter_result_list = []
                    for __range in range_list:
                        for result in result_list:
                            if result['mode'] == __range:
                                filter_result_list.append(result)

                    with zipfile.ZipFile(zip_buffer, 'w') as _zip:
                        sdf_content = ""
                        for result in filter_result_list:
                            sdf_content += result['mol'] + "$$$$\n"
                        _zip.writestr(f"{task_id}.sdf", sdf_content)

                        if zip_csv_file:
                            csv_content = cls.get_csv_data(data=filter_result_list)
                            _zip.writestr('result.csv', csv_content)

                    zip_buffer.seek(0)
                    return zip_buffer
                except Exception as e:
                    import traceback
                    traceback.print_exc()
                    logging.error(f"下载指定结果失败：{e}")
                    return None
        else:
            return None

    # todo 这是临时性的方法，后面可能还需要再深究
    @classmethod
    def start_task_for_custom_tool(cls, task_name: str, pdb_file_url: str, center_x: float, center_y: float,
                                   center_z: float,
                                   size_x: float, size_y: float, size_z: float, ligand_file_urls: list[str],
                                   out_pose_num: int,
                                   ) -> dict:
        """
        Start molecular docking task.
        :param task_name: task name
        :param pdb_file_url: pdb file url
        :param center_x: center x
        :param center_y: center y
        :param center_z: center z
        :param size_x: size x
        :param size_y: size y
        :param size_z: size z
        :param ligand_file_urls: ligand file urls
        :param out_pose_num: output pose number
        :return: molecular docking task object
        """
        pdb_file_buffer = cls.download_file(pdb_file_url)
        ligand_file_buffer_list: list[BufferedReader] = []
        for ligand_file_url in ligand_file_urls:
            ligand_file_buffer_list.append(cls.download_file(ligand_file_url))
        return cls.main_processor_for_custom_tool(task_name, pdb_file_buffer, ligand_file_buffer_list, center_x,
                                                  center_y,
                                                  center_z,
                                                  size_x, size_y, size_z, out_pose_num)

    # todo 这是临时性的方法，后面可能还需要再深究
    @classmethod
    def main_processor_for_custom_tool(cls, task_name: str, pdb_file_buffer: BufferedReader,
                                       ligand_file_buffer_list: list[BufferedReader], center_x: float, center_y: float,
                                       center_z: float,
                                       size_x: float, size_y: float, size_z: float, out_pose_num: int,
                                       ) -> dict:
        # 调用DockingProcessorAPI进行docking
        try:
            logging.info(click.style(f"自定义工具 开始分子对接任务：{task_name}", fg='blue'))
            result, success_flag = cls.docking_processor(pdb_file_buffer, ligand_file_buffer_list, center_x, center_y,
                                                         center_z,
                                                         size_x, size_y,
                                                         size_z, out_pose_num)
            if success_flag:
                logging.info(click.style(f"分子对接成功", fg='blue'))
            else:
                logging.info(click.style(f"分子对接失败", fg='red', bold=True))
        except Exception as e:
            logging.info(click.style(f"分子对接任务失败：{e}", fg='red', bold=True))
            result = {"error": str(e)}

        return result

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

    @classmethod
    def get_pocket_docking_file_path(cls, file_name: str, user: Union[Account, EndUser]):
        """
        获取口袋对接文件保存路径
        :param file_name: 文件名
        :param user:     用户
        :return: 文件路径
        """
        date_path = datetime.datetime.now().strftime('%Y/%m/%d')
        if isinstance(user, Account):
            current_tenant_id = user.current_tenant_id
        else:
            current_tenant_id = user.tenant_id
        return f"upload_files/{current_tenant_id}/pocket_docking/{date_path}/{file_name}", current_tenant_id

    @classmethod
    def get_csv_data(cls, data: list) -> str:
        if data:
            csv_io = io.StringIO()
            try:
                # 删除mol字段
                data[0].pop('mol', None)
                headers = data[0].keys()
                headers = ['mode'] + [x for x in headers if x != 'mode']

                writer = csv.DictWriter(csv_io, fieldnames=headers)
                # 写入表头
                writer.writeheader()
                for row in data:
                    row.pop('mol', None)
                    row['mode'] = f'="{row["mode"]}"'
                    writer.writerow(row)
                csv_content = csv_io.getvalue()
                return csv_content
            except Exception as e:
                logging.error(click.style(f"生成csv文件失败：{e}", fg='red', bold=True))
                raise ValueError("生成csv文件失败")
            finally:
                csv_io.close()

    @classmethod
    def get_service_result_data(cls, task_id: str, user: Union[Account, EndUser]):
        data = MolecularDockingTask.query.filter_by(id=task_id, created_by=user.id).first()
        return data.serialize


# acks_late 设置为 True 时，任务的消息确认（acknowledgement）会在任务执行完成后才发送，确保任务在失败或 worker 崩溃时能重新被执行。
# time_limit 设置为 120 秒，任务的执行时间不能超过 120 秒，超过这个时间，任务会被自动取消。
# bind 如果设置为 True，任务将绑定到当前任务实例（self），从而允许你在任务中访问 self（即任务对象本身）。这在需要访问任务元数据（例如任务ID、重试次数）时非常有用
@shared_task(queue='molecular_docking', bind=True, time_limit=120, acks_late=True)
def molecular_docking_celery_task(self, user_dict: dict, center_x: float, center_y: float,
                                  center_z: float,
                                  size_x: float, size_y: float, size_z: float, out_pose_num: int,
                                  chain: str, residue_number: int,
                                  molecular_docking_task_dict: dict):
    logging.info(click.style(f"molecular_docking_celery_task 开始执行，任务ID：{self.request.id}", fg='blue'))
    # 因为celery需要序列化之后才能传递到该参数，所以现在的user和molecular_docking_task都是json类型的，需要进行反序列化
    molecular_docking_task = MolecularDockingTask(**molecular_docking_task_dict)
    user = Account(**user_dict)

    MolecularDockingTask.query.filter_by(id=molecular_docking_task.id, created_by=user.id).update(
        {'status': Status.PROCESSING.status})
    db.session.commit()
    # 在sciminer_history_task数据表中更新任务状态
    SciminerHistoryTask.query.filter_by(task_id=molecular_docking_task.id, created_by=user.id).update(
        {'status': Status.PROCESSING.status}
    )
    db.session.commit()
    calling_websocket_internal_send(channel='molecular_docking', user_id=user.id, message={
        "id": molecular_docking_task.id,
        "task_name": molecular_docking_task.task_name,
        "result": None,
        "status": Status.PROCESSING.status,
        "remove_ligand_file": None
    })

    # 获取pdb和ligand文件buffer
    pdb_file_buffer = UploadFileUtils.get_upload_file_buffer(molecular_docking_task.pdb_file_id, user)
    ligand_file_buffer_list = UploadFileUtils.get_upload_file_buffer(
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
        chain, residue_number,
        molecular_docking_task
    )
    logging.info(click.style(f"molecular_docking_celery_task 执行完成，任务ID：{self.request.id}", fg='blue'))
