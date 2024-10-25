"""
@Time    : 2024/10/15 上午9:25
@Author  : bigboss
@Description : 全局对接路由
"""
from flask import request
from flask_login import current_user
from flask_restful import Resource, marshal_with

import services
from configs import dify_config
from controllers.console import api
from controllers.console.datasets.error import (
    FileTooLargeError,
    NoFileUploadedError,
    UnsupportedFileTypeError, TooManyFilesError,
)
from controllers.console.setup import setup_required
from controllers.console.wraps import account_initialization_required, cloud_edition_billing_resource_check
from fields.file_fields import file_fields, upload_config_fields
from libs.login import login_required
from services.sciminer_services.global_docking.global_docking_file_service import GlobalDockingFileService, GlobalDockingSourceType


class GlobalDockingFileApi(Resource):
    @setup_required
    @login_required
    @account_initialization_required
    @marshal_with(upload_config_fields)
    def get(self):
        file_size_limit = dify_config.UPLOAD_FILE_SIZE_LIMIT
        batch_count_limit = dify_config.UPLOAD_FILE_BATCH_LIMIT
        image_file_size_limit = dify_config.UPLOAD_IMAGE_FILE_SIZE_LIMIT
        return {
            "file_size_limit": file_size_limit,
            "batch_count_limit": batch_count_limit,
            "image_file_size_limit": image_file_size_limit,
        }, 200

    @setup_required
    @login_required
    @account_initialization_required
    @marshal_with(file_fields)
    @cloud_edition_billing_resource_check("documents")
    def post(self):
        # 全局对接 可以接收文件或者字符串形式的文件内容，二者选一，优先使用文件上传
        file = request.files.get("file", None)
        file_content = request.form.get("file_content", None)
        source = request.args.get("source", None)
        if source is None or source not in GlobalDockingSourceType._value2member_map_:
            raise ValueError("maybe source is none or not in allowed source type")
        else:
            source = GlobalDockingSourceType(source)
        try:
            if file:
                if len(request.files) > 1:
                    raise TooManyFilesError()
                upload_file_list = GlobalDockingFileService.upload_file(file, source, current_user)
            elif file_content:
                upload_file_list = GlobalDockingFileService.upload_file(file_content, source, current_user)
            else:
                raise NoFileUploadedError()
        except services.errors.file.FileTooLargeError as file_too_large_error:
            raise FileTooLargeError(file_too_large_error.description)
        except services.errors.file.UnsupportedFileTypeError:
            raise UnsupportedFileTypeError()

        return upload_file_list, 201


class GlobalDockingFileReadApi(Resource):
    """
    全局对接文件读取接口，根据file_id获取文件内容
    """
    @setup_required
    @login_required
    @account_initialization_required
    def get(self):
        file_id = request.args.get("file_id", None)
        if file_id is None:
            raise ValueError("maybe file_id is none")
        file_content = GlobalDockingFileService.get_file_content(file_id)
        return {"file_content": file_content}, 200


# class GetFileApi(Resource):
#     """
#     这个路由为什么没有设置需要login_required？
#     是因为这个路由是用来获取pdb文件流提供给前端mol*的渲染，需要修改源代码携带token，比较麻烦，暂时没有登录验证的必要
#     """
#     @setup_required
#     def get(self, file_id):
#         mime_type = request.args.get("mime_type", default="text/plain", type=str)
#         file_id = str(file_id)
#         file_stream = MolecularDockingFileService.get_file(file_id, mime_type)
#         return file_stream
#
#
# class RenderingMoleculeFileApi(Resource):
#     @setup_required
#     @login_required
#     @account_initialization_required
#     def get(self):
#         file_id = request.args.get("file_id", default=None, type=str)
#         return MolecularDockingFileService.rendering_molecule_file(file_id, current_user), 200


# class FileSupportTypeApi(Resource):
#     @setup_required
#     @login_required
#     @account_initialization_required
#     def get(self):
#         etl_type = dify_config.ETL_TYPE
#         allowed_extensions = UNSTRUCTURED_ALLOWED_EXTENSIONS if etl_type == "Unstructured" else ALLOWED_EXTENSIONS
#         return {"allowed_extensions": allowed_extensions}


api.add_resource(GlobalDockingFileApi, "/global-docking/files/upload")
api.add_resource(GlobalDockingFileReadApi, "/global-docking/files/read")
# api.add_resource(GetFileApi, "/global-docking/files/<uuid:file_id>")
# api.add_resource(FileSupportTypeApi, "/files/support-type")

# api.add_resource(RenderingMoleculeFileApi, "/global-docking/files/rendering")
