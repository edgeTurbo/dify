"""
@Time    : 2024/11/4 下午1:45
@Author  : bigboss
@Description : poseview任务相关接口
"""

from flask import request
from flask_login import login_required, current_user
from flask_restful import Resource, marshal_with

from controllers.console import api
from controllers.console.sciminer_apps.poseview.request_fields import PoseviewTaskRequestFields
from controllers.console.setup import setup_required
from controllers.console.wraps import account_initialization_required
from fields.sciminer_utilities.poseview.poseview_fields import poseview_task_fields
from services.sciminer_services.poseview.poseview_service import PoseViewService


class PoseViewTaskApi(Resource):

    @setup_required
    @login_required
    @account_initialization_required
    @marshal_with(poseview_task_fields)
    def post(self):
        data = request.get_json()
        poseview_request_fields = PoseviewTaskRequestFields(**data)

        poseview_task = PoseViewService.start_task(poseview_request_fields, user=current_user, start_celery=True)

        return poseview_task, 201


api.add_resource(PoseViewTaskApi, "/poseview/task")