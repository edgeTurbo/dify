"""
@Time    : 2024/10/24 下午1:25
@Author  : bigboss
@Description : 任务历史记录
"""
from flask import request
from flask_login import login_required, current_user
from flask_restful import Resource, marshal_with
from werkzeug.exceptions import Unauthorized

from controllers.console import api
from controllers.console.setup import setup_required
from controllers.console.wraps import account_initialization_required
from fields.file_fields import file_fields
from models.account import Account
from models.model import EndUser, UploadFile
from models.sciminer_models.sciminer import SciminerHistoryTask
from fields.sciminer.sciminer_fields import history_task_pagination_fields
from services.sciminer_services import service_classes


class SciminerHistoryTaskApi(Resource):
    """
    任务历史记录API
    """

    @setup_required
    @login_required
    @account_initialization_required
    @marshal_with(history_task_pagination_fields)
    def get(self):
        # 获取分页参数，默认为第 1 页，每页显示 10 条数据
        page = request.args.get('page', 1, type=int)
        page_size = request.args.get('page_size', 10, type=int)
        # 获取当前用户ID
        if isinstance(current_user, Account) or isinstance(current_user, EndUser):
            history_tasks_query = SciminerHistoryTask.query.filter_by(created_by=current_user.id)
            history_tasks_query = history_tasks_query.order_by(SciminerHistoryTask.created_at.desc())
            history_tasks = history_tasks_query.paginate(page=page, per_page=page_size, max_per_page=1000,
                                                         error_out=False)
            return history_tasks, 200
        else:
            raise Unauthorized()


class SciminerHistoryTaskDetailApi(Resource):
    """
    任务历史记录详情API
    """

    @setup_required
    @login_required
    @account_initialization_required
    def post(self):
        task_id = request.json.get("task_id")
        task_type = request.json.get("task_type")
        service_class = service_classes.get(task_type)
        if service_class is None:
            raise ValueError("Task type not found")
        data = service_class.get_service_result_data(task_id, current_user)
        return data, 200


class SciminerHistoryTskFileApi(Resource):
    """
    获取任务历史记录文件API
    """

    @setup_required
    @login_required
    @account_initialization_required
    @marshal_with(file_fields)
    def get(self):
        file_id = request.args.get("file_id")
        upload_file = UploadFile.query.filter_by(id=file_id, created_by=current_user.id).first()
        return upload_file, 200


api.add_resource(SciminerHistoryTaskApi, "/history_task")
api.add_resource(SciminerHistoryTaskDetailApi, "/history_task/detail")
api.add_resource(SciminerHistoryTskFileApi, "/history_task/file")
