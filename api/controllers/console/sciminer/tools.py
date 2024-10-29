"""
@Time    : 2024/10/28 下午4:58
@Author  : bigboss
@Description : 
"""
import io
import mimetypes
from pathlib import Path

from flask import send_file
from flask_login import login_required
from flask_restful import Resource, marshal_with
from os import path


from configs import dify_config
from controllers.console import api
from controllers.console.setup import setup_required
from controllers.console.wraps import account_initialization_required
from core.tools.errors import ToolProviderNotFoundError
from services.sciminer_services import service_yaml_list, icon_dict


class SciminerToolsApi(Resource):
    """
    sciminer tools api
    """

    @setup_required
    @login_required
    @account_initialization_required
    def get(self):
        return service_yaml_list


class SciminerToolsIconApi(Resource):
    """
    sciminer tools icon api
    """

    @setup_required
    def get(self, util_name):
        if util_name in icon_dict:
            icon_path = icon_dict[util_name]
            if not path.exists(icon_path):
                raise ToolProviderNotFoundError(f"sciminer utility {util_name} icon not found")
            mime_type, _ = mimetypes.guess_type(icon_path)
            mime_type = mime_type or "application/octet-stream"
            icon_bytes = Path(icon_path).read_bytes()
            icon_cache_max_age = dify_config.TOOL_ICON_CACHE_MAX_AGE
            return send_file(io.BytesIO(icon_bytes), mimetype=mime_type, max_age=icon_cache_max_age)
        else:
            raise ToolProviderNotFoundError(f"sciminer utility {util_name} icon not found")


api.add_resource(SciminerToolsApi, '/sciminer/tools')
api.add_resource(SciminerToolsIconApi, "/sciminer/utils/<util_name>/icon")
