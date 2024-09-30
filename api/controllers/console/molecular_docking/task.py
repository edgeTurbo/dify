import json
import uuid

from flask import request
from flask_login import current_user
from flask_restful import Resource, marshal_with

import services
from configs import dify_config
from controllers.console import api
from controllers.console.molecular_docking.error import (
    IllegalParametersError,
)
from controllers.console.setup import setup_required
from controllers.console.wraps import account_initialization_required, cloud_edition_billing_resource_check
from fields.molecular_docking_fields import molecular_docking_task_fields
from libs.login import login_required
from services.molecular_docking.molecular_docking_service import MolecularDockingService
from configs.websocket_config import websocket_handler


class MolecularDockingTaskApi(Resource):

    @setup_required
    @login_required
    @account_initialization_required
    @marshal_with(molecular_docking_task_fields)
    # @cloud_edition_billing_resource_check("documents")
    def post(self):
        task_name = request.form.get("task_name") if request.form.get("task_name") else "molecular_docking_task_" + str(
            uuid.uuid4())

        pdb_file_id = request.form.get("pdb_file_id")

        center_x = float(request.form.get("center_x"))
        center_y = float(request.form.get("center_y"))
        center_z = float(request.form.get("center_z"))

        size_x = float(request.form.get("size_x"))
        size_y = float(request.form.get("size_y"))
        size_z = float(request.form.get("size_z"))

        ligand_file_ids = request.form.get("ligand_file_ids").split(",") if request.form.get(
            "ligand_file_ids") else None

        # 输出结果数目
        out_pose_num = int(request.form.get("out_pose_num"))

        if pdb_file_id is None or center_x is None or center_y is None or center_z is None or size_x is None or size_y is None or size_z is None or ligand_file_ids is None or out_pose_num is None:
            raise IllegalParametersError()

        molecular_docking_task = MolecularDockingService.start_task(task_name, pdb_file_id, center_x, center_y,
                                                                    center_z, size_x, size_y, size_z, ligand_file_ids,
                                                                    out_pose_num, current_user)

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


api.add_resource(MolecularDockingTaskApi, "/molecular-docking/task")
api.add_resource(MolecularDockingCenterPositionApi, "/molecular-docking/center-position")
