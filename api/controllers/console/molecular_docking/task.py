import io
import json
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
from fields.molecular_docking_fields import molecular_docking_task_fields
from libs.login import login_required
from services.molecular_docking.molecular_docking_service import MolecularDockingService


class MolecularDockingTaskApi(Resource):

    @setup_required
    @login_required
    @account_initialization_required
    @marshal_with(molecular_docking_task_fields)
    # @cloud_edition_billing_resource_check("documents")
    def post(self):
        data = request.get_json()
        task_name = data["task_name"] if data["task_name"] else "molecular_docking_task_" + str(
            uuid.uuid4())

        pdb_file_id = data["pdb_file_id"]

        center_x = float(data["center_x"])
        center_y = float(data["center_y"])
        center_z = float(data["center_z"])

        size_x = float(data["size_x"])
        size_y = float(data["size_y"])
        size_z = float(data["size_z"])

        ligand_file_ids = data["ligand_file_ids"].split(",") if data[
            "ligand_file_ids"] else None

        # 输出结果数目
        out_pose_num = int(data["out_pose_num"])

        if pdb_file_id is None or center_x is None or center_y is None or center_z is None or size_x is None or size_y is None or size_z is None or ligand_file_ids is None or out_pose_num is None:
            raise IllegalParametersError()

        molecular_docking_task = MolecularDockingService.start_task(task_name, pdb_file_id, center_x, center_y,
                                                                    center_z, size_x, size_y, size_z, ligand_file_ids,
                                                                    out_pose_num, current_user, start_celery=False)
        return molecular_docking_task, 201


class MolecularDockingCenterPositionApi(Resource):
    """
    分子对接获取中心点坐标接口
    """

    @setup_required
    @login_required
    @account_initialization_required
    def post(self):
        data = request.get_json()
        pdb_file_id = data["pdb_file_id"]
        if pdb_file_id is None:
            raise IllegalParametersError()
        center_position = MolecularDockingService.get_center_position(pdb_file_id, current_user)
        return center_position, 200


class MolecularDockingTaskResultDownloadApi(Resource):
    """
    分子对接任务结果下载接口
    """

    @setup_required
    @login_required
    @account_initialization_required
    def get(self):
        task_id = request.args.get("task_id", default=None, type=str)
        _range = request.args.get("range", default="all", type=str)
        sdf_content = MolecularDockingService.download_task_result(task_id, _range, current_user)
        if sdf_content is None:
            return {"message": "Task not found"}, 500
        else:
            buffer = io.BytesIO()
            # 如果是文本内容，进行编码
            buffer.write(sdf_content.encode('utf-8'))
            # 重置指针到文件开头
            buffer.seek(0)
            return send_file(buffer, as_attachment=True, download_name=f"{task_id}.sdf", mimetype='application/octet-stream')


class MolecularDockingTaskCustomToolsApi(Resource):
    """
    提供给工作流的自定义工具，分子对接任务接口
    """
    def post(self):
        data = request.get_json()
        task_name = data["task_name"] if data["task_name"] else "molecular_docking_task_" + str(
            uuid.uuid4())

        pdb_file_url = data["pdb_file_url"]

        center_x = float(data["center_x"])
        center_y = float(data["center_y"])
        center_z = float(data["center_z"])

        size_x = float(data["size_x"])
        size_y = float(data["size_y"])
        size_z = float(data["size_z"])

        ligand_file_urls = data["ligand_file_urls"].split(",") if data[
            "ligand_file_urls"] else None

        # 输出结果数目
        out_pose_num = int(data["out_pose_num"])

        if pdb_file_url is None or center_x is None or center_y is None or center_z is None or size_x is None or size_y is None or size_z is None or ligand_file_urls is None or out_pose_num is None:
            raise IllegalParametersError()

        result = MolecularDockingService.start_task_for_custom_tool(task_name, pdb_file_url, center_x, center_y,
                                                                    center_z, size_x, size_y, size_z, ligand_file_urls,
                                                                    out_pose_num)
        return result, 201


api.add_resource(MolecularDockingTaskApi, "/molecular-docking/task")
api.add_resource(MolecularDockingCenterPositionApi, "/molecular-docking/center-position")
api.add_resource(MolecularDockingTaskResultDownloadApi, "/molecular-docking/download")

api.add_resource(MolecularDockingTaskCustomToolsApi, "/custom-tools/molecular-docking/task")
