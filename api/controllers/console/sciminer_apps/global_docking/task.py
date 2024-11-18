import json
import uuid

from flask import request, send_file
from flask_login import current_user
from flask_restful import Resource, marshal_with

from configs.websocket_config import websocket_handler
from controllers.console import api
from controllers.console.sciminer_apps.error import (
    IllegalParametersError,
)
from controllers.console.setup import setup_required
from controllers.console.wraps import account_initialization_required
from extensions.ext_database import db
from fields.global_docking_fields import global_docking_task_fields
from libs.login import login_required
from models.sciminer_models.global_docking import GlobalDockingTask
from models.sciminer_models.sciminer import Status, SciminerHistoryTask
from services.sciminer_services.global_docking.global_docking_service import GlobalDockingService


class GlobalDockingTaskApi(Resource):

    @setup_required
    @login_required
    @account_initialization_required
    @marshal_with(global_docking_task_fields)
    # @cloud_edition_billing_resource_check("documents")
    def post(self):
        data = request.get_json()
        task_name = data["task_name"] if data["task_name"] else "global_docking_task_" + str(
            uuid.uuid4())

        fasta_file_id = data["fasta_file_id"]

        ligand_file_ids = data["ligand_file_ids"].split(",") if data[
            "ligand_file_ids"] else None

        # 输出结果数目
        out_pose_num = int(data["out_pose_num"])

        if fasta_file_id is None or ligand_file_ids is None or out_pose_num is None:
            raise IllegalParametersError()

        global_docking_task = GlobalDockingService.start_task(task_name, fasta_file_id, ligand_file_ids,
                                                              out_pose_num, current_user,
                                                              start_celery=True)
        return global_docking_task, 201


class GlobalDockingTaskResultCallbackApi(Resource):
    """
    分子对接任务结果回调接口(供模型返回结果时调用)
    """

    def post(self):
        data = request.get_json()
        message = data.get("message", None)
        user_id = data.get("user_id", None)
        task_id = data.get("task_id", None)
        global_docking_task = GlobalDockingTask.query.filter_by(id=task_id, created_by=user_id).first()
        if global_docking_task.status != Status.FAILURE.status:
            # 只有当任务状态不是失败状态时才更新任务结果
            global_docking_task.result = message
            global_docking_task.status = Status.SUCCESS.status
            db.session.commit()

            # 在sciminer_history_task数据表中更新任务状态
            SciminerHistoryTask.query.filter_by(task_id=task_id, created_by=user_id).update(
                {'status': Status.SUCCESS.status}
            )
            db.session.commit()
            # 发送websocket消息通知前端任务状态变化
            websocket_handler.send_message_to_user(
                uid=user_id,
                message={
                    "global_docking": {
                        "id": task_id,
                        "task_name": global_docking_task.task_name,
                        "result": message,
                        "status": Status.SUCCESS.status,
                    }
                },
            )


class GlobalDockingTaskResultDownloadApi(Resource):
    """
    分子对接任务结果下载接口
    """

    @setup_required
    @login_required
    @account_initialization_required
    def get(self):
        task_id = request.args.get("task_id", default=None, type=str)
        _range = request.args.get("range", default="all", type=str)
        zip_buffer = GlobalDockingService.download_task_result(task_id, _range, current_user)
        if zip_buffer is None:
            return {"message": "Task not found"}, 500
        else:
            return send_file(zip_buffer, as_attachment=True, download_name=f"{task_id}.zip",
                             mimetype='application/octet-stream')


api.add_resource(GlobalDockingTaskApi, "/global-docking/task")
api.add_resource(GlobalDockingTaskResultDownloadApi, "/global-docking/download")

api.add_resource(GlobalDockingTaskResultCallbackApi, "/global-docking/task/callback")
