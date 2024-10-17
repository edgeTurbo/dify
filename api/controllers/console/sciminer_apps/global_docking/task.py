import io
import uuid

from flask import request, send_file
from flask_login import current_user
from flask_restful import Resource, marshal_with

from controllers.console import api
from controllers.console.molecular_docking.error import (
    IllegalParametersError,
)
from controllers.console.setup import setup_required
from controllers.console.wraps import account_initialization_required, cloud_edition_billing_resource_check
from fields.global_docking_fields import global_docking_task_fields
from libs.login import login_required
from services.global_docking.global_docking_service import GlobalDockingService


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
                                                              start_celery=False)
        return global_docking_task, 201


# class GlobalDockingTaskResultDownloadApi(Resource):
#     """
#     分子对接任务结果下载接口
#     """
#
#     @setup_required
#     @login_required
#     @account_initialization_required
#     def get(self):
#         task_id = request.args.get("task_id", default=None, type=str)
#         _range = request.args.get("range", default="all", type=str)
#         sdf_content = MolecularDockingService.download_task_result(task_id, _range, current_user)
#         if sdf_content is None:
#             return {"message": "Task not found"}, 500
#         else:
#             buffer = io.BytesIO()
#             # 如果是文本内容，进行编码
#             buffer.write(sdf_content.encode('utf-8'))
#             # 重置指针到文件开头
#             buffer.seek(0)
#             return send_file(buffer, as_attachment=True, download_name=f"{task_id}.sdf",
#                              mimetype='application/octet-stream')


api.add_resource(GlobalDockingTaskApi, "/global-docking/task")
# api.add_resource(GlobalDockingTaskResultDownloadApi, "/global-docking/download")
