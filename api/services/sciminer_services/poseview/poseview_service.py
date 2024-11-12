"""
@Time    : 2024/11/4 上午11:02
@Author  : bigboss
@Description : poseview service
"""
import hashlib
import json
import logging
import os
from mimetypes import guess_type
from typing import Union, Tuple

import click
import requests
from celery import shared_task

from configs import dify_config
from controllers.console.sciminer_apps.poseview.request_fields import PoseviewTaskRequestFields
from controllers.inner_api.websocket.websocket import calling_websocket_internal_send
from core.tools.utils.result_file_utils import ResultFileUtils
from core.tools.utils.upload_file_utils import UploadFileUtils, UploadFileProperty
from extensions.ext_database import db
from libs import rdkit_lib
from models.account import Account
from models.model import EndUser
from models.sciminer_models.poseview import PoseViewTask
from models.sciminer_models.sciminer import Status, SciminerHistoryTask
from services.sciminer_services import SciminerBaseService
from services.sciminer_services.poseview.poseview_file_service import PoseViewFileService

if dify_config.POSEVIEW_TASK_API_URL == "" or dify_config.POSEVIEW_TASK_API_URL is None:
    logging.error(
        click.style(
            "poseview任务API URL不能为空, 请在配置文件.env中设置POSEVIEW_TASK_API_URL, 否则poseview功能无法正常使用",
            fg='red', bold=True))

if dify_config.POSEVIEW_RESULT_API_URL == "" or dify_config.POSEVIEW_RESULT_API_URL is None:
    logging.error(
        click.style(
            "poseview结果查询API URL不能为空, 请在配置文件.env中设置POSEVIEW_RESULT_API_URL, 否则poseview功能无法正常使用",
            fg='red', bold=True))

if dify_config.INNER_API is None or dify_config.INNER_API_KEY is None:
    logging.error(
        click.style(
            "内网API地址或API KEY未设置, 请在配置文件.env中设置INNER_API为true和INNER_API_KEY, 否则poseview功能无法正常使用，websocket消息队列无法正常使用",
            fg='red', bold=True)
    )


class PoseViewService(SciminerBaseService):
    service_type = "POSEVIEW"
    task_label = "PoseView"

    @classmethod
    def start_task(cls, poseview_request_fields: PoseviewTaskRequestFields, user: Union[Account, EndUser],
                   start_celery: bool = False) -> PoseViewTask:
        status = Status.PROCESSING.status
        if start_celery:
            status = Status.PENDING.status
            logging.info(click.style(f"poseview使用消息队列进行任务处理", fg='green', bold=True))

        # 初始化poseview任务。保存到数据表中
        poseview_task = PoseViewTask(
            receptor_file_id=poseview_request_fields.receptor_file_id,
            ligand_file_ids=poseview_request_fields.ligand_file_ids,
            status=status,
            created_by=user.id,
            task_name=poseview_request_fields.task_name,
        )
        db.session.add(poseview_task)
        db.session.commit()

        # 新增任务历史记录
        sciminer_history_task = SciminerHistoryTask(
            task_id=poseview_task.id,
            task_name=poseview_request_fields.task_name,
            task_type=cls.service_type,
            label=cls.task_label,
            status=status,
            created_by=user.id,
        )
        db.session.add(sciminer_history_task)
        db.session.commit()

        if start_celery:
            # 启动消息队列进行任务处理
            poseview_celery_task.apply_async(
                args=[user.serialize, poseview_task.serialize],
            )
            return poseview_task
        else:
            receptor_file_property = UploadFileUtils.get_upload_file_bytes(poseview_request_fields.receptor_file_id,
                                                                           user)
            ligand_file_property_list = UploadFileUtils.get_upload_file_bytes(poseview_request_fields.ligand_file_ids,
                                                                              user, return_list=True)
            return cls.main_processor(user, poseview_request_fields.task_name, receptor_file_property,
                                      ligand_file_property_list, poseview_task, start_celery)

    @classmethod
    def get_service_result_data(cls, task_id: str, user: Union[Account, EndUser]):
        pass

    @classmethod
    def main_processor(cls, user: Union[Account, EndUser], task_name: str, receptor_file_property: UploadFileProperty,
                       ligand_file_property_list: list[UploadFileProperty], poseview_task: PoseViewTask,
                       start_celery: bool) -> PoseViewTask:
        # 调用poseview的API进行任务处理
        try:
            logging.info(click.style(f"{user.name} 开始poseview任务：{task_name}", fg='green'))
            result, success_flag = cls.poseview_processor(user, receptor_file_property, ligand_file_property_list)
            if isinstance(result, list):
                result = json.dumps(result)
            if success_flag:
                logging.info(click.style(f"poseview任务成功", fg='green'))
                status = Status.SUCCESS.status
            else:
                logging.info(click.style(f"poseview任务失败", fg='red', bold=True))
                status = Status.FAILURE.status
        except Exception as e:
            import traceback
            traceback.print_exc()
            logging.info(click.style(f"poseview任务失败：{e}", fg='red', bold=True))
            status = Status.FAILURE.status
            result = e

        # 更新utility_poseview_task数据表中的结果信息
        poseview_task = PoseViewTask.query.filter_by(id=poseview_task.id, created_by=user.id).first()
        poseview_task.result = result
        poseview_task.status = status
        poseview_task.task_name = task_name
        db.session.commit()

        # 在sciminer_history_task数据表中更新任务状态
        SciminerHistoryTask.query.filter_by(task_id=poseview_task.id, created_by=user.id).update(
            {'status': status}
        )
        db.session.commit()

        if start_celery:
            # 当启动消息队列的时候才进行发送websocket消息
            calling_websocket_internal_send(channel='poseview', user_id=user.id, message={
                "id": poseview_task.id,
                "task_name": task_name,
                "result": result,
                "status": status,
            })

        return poseview_task

    @classmethod
    def poseview_processor(cls, user: Union[Account, EndUser], receptor_file_property: UploadFileProperty,
                           ligand_file_property_list: list[UploadFileProperty]) -> Tuple[
        Union[list, str], bool]:
        # poseview的API调用代码
        try:
            result_file_list = []
            for ligand_file_property in ligand_file_property_list:
                # 用来保存所有小分子bytes的列表
                sdf_bytes_list = []
                mol_list = rdkit_lib.read_sdf_bytes(ligand_file_property.file_bytes)
                ligand_mol_num = 0
                for mol_object in mol_list:
                    if mol_object is not None:
                        sdf_bytes = rdkit_lib.mol_object_to_sdf_bytes(mol_object)
                        if sdf_bytes is not None:
                            sdf_bytes_list.append(sdf_bytes)
                    ligand_mol_num += 1
                logging.info(click.style(
                    f"单个ligand文件中一共有{ligand_mol_num}个小分子，其中可被使用的有{len(sdf_bytes_list)}个",
                    fg='green'))

                # 调用poseview的API进行任务处理
                for sdf_bytes in sdf_bytes_list:
                    files = [
                        ('protein_file', (receptor_file_property.file_name, receptor_file_property.file_bytes,
                                          'application/octet-stream')),
                        ('ligand_file', (ligand_file_property.file_name, sdf_bytes, 'application/octet-stream'))
                    ]
                    task_response = requests.post(dify_config.POSEVIEW_TASK_API_URL, files=files)
                    task_json = task_response.json()
                    job_id = task_json.get('job_id')
                    logging.info(click.style(f"job_id={job_id}", fg='green'))
                    if job_id is None:
                        return "job_id为空，poseview任务失败", False
                    # 查询poseview任务结果
                    result_response = requests.get(dify_config.POSEVIEW_RESULT_API_URL.format(job_id))
                    result_json = result_response.json()
                    logging.info(click.style(f"poseview api返回的结果：{result_json}", fg='green'))
                    result_status = result_json.get('status')
                    if result_status == 'success':
                        result_image_url = result_json.get('image')
                        image_response = requests.get(result_image_url)
                        image_content = image_response.content
                        file_size = len(image_content)
                        file_name = os.path.basename(result_image_url)
                        # 保存poseview结果图片到服务器
                        file_key, current_tenant_id, lower_file_extension = PoseViewFileService.get_poseview_file_key(
                            file_name=file_name, user=user)
                        ResultFileUtils.result_file_to_storage(
                            file_key=file_key,
                            file_content=image_content,
                            extension=lower_file_extension,
                        )
                        result_file = ResultFileUtils.add_result_file_to_db(
                            tenant_id=current_tenant_id,
                            storage_type=dify_config.STORAGE_TYPE,
                            key=file_key,
                            name=file_name,
                            size=file_size,
                            extension=lower_file_extension,
                            mime_type=guess_type(file_name)[0] or "octet/stream",
                            created_by_role=("account" if isinstance(user, Account) else "end_user"),
                            created_by=user.id,
                            used=False,
                            hash=hashlib.sha3_256(image_content).hexdigest(),
                            source=cls.service_type,
                        )
                        result_file_list.append(
                            {
                                "file_id": result_file.id,
                            }
                        )
                    else:
                        return "不是成功状态，poseview任务失败", False
            return result_file_list, True
        except Exception as e:
            logging.info(click.style(f"poseview任务调用失败：{e}", fg='red', bold=True))
            return "poseview任务调用失败", False


# acks_late 设置为 True 时，任务的消息确认（acknowledgement）会在任务执行完成后才发送，确保任务在失败或 worker 崩溃时能重新被执行。
# time_limit 设置为 120 秒，任务的执行时间不能超过 120 秒，超过这个时间，任务会被自动取消。
# bind 如果设置为 True，任务将绑定到当前任务实例（self），从而允许你在任务中访问 self（即任务对象本身）。这在需要访问任务元数据（例如任务ID、重试次数）时非常有用
@shared_task(queue='poseview', bind=True, time_limit=600, acks_late=True)
def poseview_celery_task(self, user_dict: dict, poseview_task_dict: dict):
    logging.info(click.style(f"poseview_celery_task 开始执行，任务ID：{self.request.id}", fg='blue'))
    # 因为celery需要序列化之后才能传递到该参数，所以现在的user和molecular_docking_task都是json类型的，需要进行反序列化
    poseview_task = PoseViewTask(**poseview_task_dict)
    user = Account(**user_dict)

    # 在sciminer_history_task数据表中更新任务状态
    SciminerHistoryTask.query.filter_by(task_id=poseview_task.id, created_by=user.id).update(
        {'status': Status.PROCESSING.status}
    )
    db.session.commit()
    calling_websocket_internal_send(channel='poseview', user_id=user.id, message={
        "id": poseview_task.id,
        "task_name": poseview_task.task_name,
        "result": None,
        "status": Status.PROCESSING.status,
    })

    receptor_file_property = UploadFileUtils.get_upload_file_bytes(poseview_task.receptor_file_id,
                                                                   user)
    ligand_file_property_list = UploadFileUtils.get_upload_file_bytes(poseview_task.ligand_file_ids,
                                                                      user, return_list=True)

    PoseViewService.main_processor(user, poseview_task.task_name, receptor_file_property,
                                   ligand_file_property_list, poseview_task, start_celery=True)
    logging.info(click.style(f"poseview_celery_task 执行完成，任务ID：{self.request.id}", fg='blue'))
