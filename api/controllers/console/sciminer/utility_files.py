"""
@Time    : 2024/11/13 上午8:57
@Author  : bigboss
@Description : 用于sciminer工具文件预览的API
"""
import traceback
from enum import Enum
from io import BytesIO

from flask import request, send_file
from flask_restful import Resource

from controllers.console import api
from controllers.console.setup import setup_required
from extensions.ext_storage import storage
from libs.exception import BaseHTTPException
from models.sciminer_models.sciminer import ResultFile


class PreviewType(Enum):
    STREAM = "stream"
    TEXT = "text"


class PreviewFileSource(Enum):
    """
    文件来源, 目前支持上传文件和结果文件
    """
    UPLOAD_FILE = "upload"
    RESULT_FILE = "result"


class UtilityFilePreviewApi(Resource):
    """
    用于sciminer工具文件预览的API
    工具 - 文件预览 API
    """

    @setup_required
    def get(self, file_name: str):
        """
        preview_type 用来判断预览类型，返回给前端的预览内容，目前支持 stream 和 text
        """
        preview_type_str = request.args.get("preview_type", "stream")
        file_id = request.args.get("file_id", None)
        source_str = request.args.get("source", "upload")
        if file_id is None:
            raise UnsupportedFileIdError()
        try:
            preview_type = PreviewType(preview_type_str)
        except ValueError:
            raise UnsupportedPreviewTypeError()
        try:
            source = PreviewFileSource(source_str)
        except ValueError:
            raise UnsupportedSourceError()

        if source == PreviewFileSource.UPLOAD_FILE:
            print("preview upload file")
        elif source == PreviewFileSource.RESULT_FILE:
            if preview_type == PreviewType.STREAM:
                # 查询result_file数据表，返回给前端文件流
                result_file = ResultFile.query.filter_by(id=file_id).first()
                if result_file is None:
                    raise UnsupportedFileError()
                file_content = storage.load_once(result_file.key)
                return send_file(BytesIO(file_content), mimetype=result_file.mime_type)


api.add_resource(UtilityFilePreviewApi, "/files/utility/<string:file_name>")


class UnsupportedPreviewTypeError(BaseHTTPException):
    error_code = "unsupported_preview_type"
    description = "Preview type not allowed."
    code = 500


class UnsupportedFileIdError(BaseHTTPException):
    error_code = "unsupported_file_id"
    description = "File id not allowed."
    code = 500


class UnsupportedSourceError(BaseHTTPException):
    error_code = "unsupported_source"
    description = "Source not allowed."
    code = 500


class UnsupportedFileError(BaseHTTPException):
    error_code = "unsupported_file"
    description = "File not allowed."
    code = 500
