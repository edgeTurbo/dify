"""
@Time    : 2024/11/4 下午1:45
@Author  : bigboss
@Description : poseview任务相关接口
"""

from flask import request, send_file
from flask_login import login_required, current_user
from flask_restful import Resource, marshal_with

from controllers.console import api
from controllers.console.sciminer_apps.molcraft.request_fields import MolcraftTaskRequestFields
from controllers.console.setup import setup_required
from controllers.console.wraps import account_initialization_required
from fields.sciminer_utilities.poseview.poseview_fields import poseview_task_fields
from services.sciminer_services.molcraft.molcraft_service import MolcraftService


class MolcraftTaskApi(Resource):
    """
    Molcraft任务提交接口
    """

    @setup_required
    @login_required
    @account_initialization_required
    @marshal_with(poseview_task_fields)
    def post(self):
        data = request.get_json()
        molcraft_request_fields = MolcraftTaskRequestFields(**data)

        molcraft_task = MolcraftService.start_task(molcraft_request_fields, user=current_user, start_celery=False)

        return molcraft_task, 201


# class PoseviewTaskResultDownloadApi(Resource):
#     """
#     poseview任务结果下载接口
#     """
#     @setup_required
#     @login_required
#     @account_initialization_required
#     def get(self):
#         task_id = request.args.get("task_id", default=None, type=str)
#         if task_id is None:
#             return {"message": "Task id not found"}, 500
#         zip_buffer = PoseViewService.download_task_result(task_id, current_user)
#         if zip_buffer is None:
#             return {"message": "Task not found"}, 500
#         else:
#             return send_file(zip_buffer, as_attachment=True, download_name=f"{task_id}.zip",
#                              mimetype='application/octet-stream')


api.add_resource(MolcraftTaskApi, "/molcraft/task")
# api.add_resource(PoseviewTaskResultDownloadApi, "/poseview/download")
