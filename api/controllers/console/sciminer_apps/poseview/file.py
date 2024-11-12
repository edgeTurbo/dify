"""
@Time    : 2024/11/4 上午10:14
@Author  : bigboss
@Description : poseview 文件上传接口
"""
from flask import request
from flask_login import login_required, current_user
from flask_restful import Resource, marshal_with

from services import errors as service_errors
from controllers.console import api
from controllers.console.setup import setup_required
from controllers.console.wraps import account_initialization_required
from controllers.service_api.app.error import NoFileUploadedError, TooManyFilesError, UnsupportedFileTypeError, InvalidFileTypeError
from fields.file_fields import file_fields
from services.sciminer_services.poseview.poseview_file_service import PoseViewFileService, PoseViewSourceType


class PoseViewFileApi(Resource):
    """
    2D相互作用图APP 文件上传接口
    """

    @setup_required
    @login_required
    @account_initialization_required
    @marshal_with(file_fields)
    def post(self):
        # PoseView 使用文件上传
        file = request.files.get("file", None)
        # 来源 1: receptor 2: ligand 用来区分上传的蛋白质或分子来限制上传文件类型
        source = request.args.get("source", "receptor")
        try:
            source_type = PoseViewSourceType(source)
        except ValueError:
            raise InvalidFileTypeError()

        # check file
        if "file" not in request.files:
            raise NoFileUploadedError()

        if len(request.files) > 1:
            raise TooManyFilesError()
        try:
            upload_file = PoseViewFileService.upload_file(file, current_user, source_type)
        except service_errors.file.UnsupportedFileTypeError:
            raise UnsupportedFileTypeError()

        return upload_file, 201


class PoseViewTaskResultFileApi(Resource):
    """
    PoseView任务结果文件查询接口（生成出来的svg图片）
    """

    def get(self, file_name: str):
        # file_name参数不需要用到，只是为了前端能够直接使用img标签进行访问，因为src属性需要特定的文件名才能够进行缓存
        file_id = request.args.get("file_id", default=None, type=str)
        if file_id is None:
            raise ValueError("maybe file_id is none")
        return PoseViewFileService.read_result_file(file_id)


api.add_resource(PoseViewFileApi, "/poseview/files/upload")
api.add_resource(PoseViewTaskResultFileApi, "/poseview/files/read/<file_name>")
