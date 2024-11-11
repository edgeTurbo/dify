"""
@Time    : 2024/10/28 下午4:58
@Author  : bigboss
@Description : 
"""
import io
from pathlib import Path

from flask import send_file
from flask_login import login_required
from flask_restful import Resource, marshal_with

from configs import dify_config
from controllers.console import api
from controllers.console.setup import setup_required
from controllers.console.wraps import account_initialization_required
from core.model_runtime.utils.encoders import jsonable_encoder
from core.tools.errors import ToolProviderNotFoundError
from services.sciminer_services import sciminer_util_manager_dict, sciminer_util_provider_list, SciminerUtilManager
from core.sciminer_utility.entities.values import default_utility_labels


class SciminerToolsApi(Resource):
    """
    sciminer tools api
    """

    @setup_required
    @login_required
    @account_initialization_required
    def get(self):
        return jsonable_encoder(
            [
                provider.to_dict()
                for provider in sciminer_util_provider_list
            ]
        )


class SciminerToolsIconApi(Resource):
    """
    sciminer tools icon api
    """

    @setup_required
    def get(self, util_name):
        if util_name in sciminer_util_manager_dict:
            sciminer_util_manager = sciminer_util_manager_dict[util_name]
            if isinstance(sciminer_util_manager, SciminerUtilManager):
                icon_path, mime_type = sciminer_util_manager.get_util_icon_path(util_name=util_name)
                icon_bytes = Path(icon_path).read_bytes()
                icon_cache_max_age = dify_config.TOOL_ICON_CACHE_MAX_AGE
                return send_file(io.BytesIO(icon_bytes), mimetype=mime_type, max_age=icon_cache_max_age)
            else:
                raise ToolProviderNotFoundError(f"sciminer utility {util_name} not found")
        else:
            raise ToolProviderNotFoundError(f"sciminer utility {util_name} icon not found")


class SciminerToolLabelsApi(Resource):
    @setup_required
    @login_required
    @account_initialization_required
    def get(self):
        return jsonable_encoder(default_utility_labels)


api.add_resource(SciminerToolsApi, '/sciminer/tools')
api.add_resource(SciminerToolsIconApi, "/sciminer/utils/<util_name>/icon")

api.add_resource(SciminerToolLabelsApi, "/sciminer/util-labels")
