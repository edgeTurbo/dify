"""
@Time    : 2024/10/25 下午3:45
@Author  : bigboss
@Description : 
"""
import importlib
import mimetypes
import os
import pkgutil
from os import path
from typing import Any, Type

from configs import dify_config
from core.helper.position_helper import get_position_map, sort_by_position_map
from core.tools.entities.api_entities import UserToolProvider
from core.tools.entities.common_entities import I18nObject
from core.tools.entities.tool_entities import ToolProviderType, ToolLabelEnum
from core.tools.entities.values import default_tool_label_dict
from core.tools.errors import ToolProviderNotFoundError
from core.tools.provider.tool_provider import ToolProviderController
from core.tools.tool.tool import Tool
from core.tools.utils.yaml_utils import load_yaml_file
from services.sciminer_services.sciminer_base_service import SciminerBaseService


class SciminerUtilManager(ToolProviderController):

    def __init__(self, **data: Any) -> None:
        if self.provider_type in {ToolProviderType.API, ToolProviderType.APP}:
            super().__init__(**data)
            return

        # load sciminer utils yaml
        service_name = data.get("service_name")
        service_yaml_file = data.get("service_yaml_file")
        try:
            service_yaml_info = load_yaml_file(service_yaml_file, ignore_error=False)
        except Exception as e:
            raise ToolProviderNotFoundError(f"can not load provider yaml for {service_name}: {e}")

        super().__init__(
            **{
                "identity": service_yaml_info["identity"],
                "credentials_schema": service_yaml_info.get("credentials_for_provider", None),
            }
        )

    def get_tools(self) -> list[Tool]:
        pass

    def get_tool(self, tool_name: str) -> Tool:
        pass

    @property
    def tool_labels(self) -> list[str]:
        """
        returns the labels of the provider

        :return: labels of the provider
        """
        label_enums = self._get_tool_labels()
        return [default_tool_label_dict[label].name for label in label_enums]

    def _get_tool_labels(self) -> list[ToolLabelEnum]:
        """
        returns the labels of the provider
        """
        return self.identity.tags or []

    def get_user_provider(self) -> UserToolProvider:
        icon_url = self.get_icon_url()
        user_provider = UserToolProvider(
            id=self.identity.name,
            author=self.identity.author,
            name=self.identity.name,
            frontend_url=self.identity.frontend_url,
            description=I18nObject(
                en_US=self.identity.description.en_US,
                zh_Hans=self.identity.description.zh_Hans,
                pt_BR=self.identity.description.pt_BR,
                ja_JP=self.identity.description.ja_JP,
            ),
            icon=icon_url,
            label=I18nObject(
                en_US=self.identity.label.en_US,
                zh_Hans=self.identity.label.zh_Hans,
                pt_BR=self.identity.label.pt_BR,
                ja_JP=self.identity.label.ja_JP,
            ),
            type=ToolProviderType.BUILT_IN,
            masked_credentials={},
            is_team_authorization=False,
            tools=[],
            labels=self.tool_labels,
        )

        return user_provider

    def get_icon_url(self):
        if self.identity.icon:
            return f"{dify_config.CONSOLE_API_URL}/console/api/sciminer/utils/{self.identity.name}/icon"

    def get_util_icon_path(self, util_name: str):
        """
        获取工具的icon路径
        :param util_name:
        :return: icon路径
        """
        absolute_path = path.join(
            path.dirname(path.realpath(__file__)),
            util_name,
            "_assets",
            self.identity.icon,
        )

        # check if the icon exists
        if not path.exists(absolute_path):
            raise ToolProviderNotFoundError(f"sciminer utility {util_name} icon not found")

        # get the mime type
        mime_type, _ = mimetypes.guess_type(absolute_path)
        mime_type = mime_type or "application/octet-stream"

        return absolute_path, mime_type

    @classmethod
    def sort_utils(cls, utils_list: list[UserToolProvider]):
        """
        排序工具列表
        """
        _position = get_position_map(os.path.dirname(__file__))

        def name_func(util: UserToolProvider):
            return util.name

        sorted_providers = sort_by_position_map(_position, utils_list, name_func)
        return sorted_providers


def generate_service_dict() -> (dict[Any, Type[SciminerBaseService]], dict):
    # 这个是给动态调用服务类用的
    _service_classes = {}
    # 下面两个是给前端查询服务列表用的
    _sciminer_util_manager_dict: dict[str, SciminerUtilManager] = {}
    _sciminer_util_provider_list: list[UserToolProvider] = []
    # 在services目录下，遍历所有子目录，当子目录的名称为sciminer_services时，导入该目录下的所有模块，并遍历模块中的所有类，
    package = 'services'
    for _, module_name, _ in pkgutil.iter_modules([package]):
        if module_name == 'sciminer_services':
            package_module = importlib.import_module(f"{package}.{module_name}")
            for _, service_name, is_pkg in pkgutil.walk_packages(list(package_module.__path__)):
                if is_pkg:
                    module = importlib.import_module(f"{package}.{module_name}.{service_name}.{service_name}_service")
                    for _name in dir(module):
                        _cls = getattr(module, _name)
                        if isinstance(_cls, type) and issubclass(_cls,
                                                                 SciminerBaseService) and _cls is not SciminerBaseService:
                            _service_classes[_cls.service_type] = _cls

                    service_yaml_file = path.join(path.dirname(module.__file__), f"{service_name}.yaml")
                    _sciminer_util_manager = SciminerUtilManager(
                        service_name=service_name,
                        service_yaml_file=service_yaml_file,
                    )
                    _sciminer_util_manager_dict[_sciminer_util_manager.identity.name] = _sciminer_util_manager

                    _sciminer_util_provider_list.append(_sciminer_util_manager.get_user_provider())

    # 排序工具列表
    return _service_classes, _sciminer_util_manager_dict, SciminerUtilManager.sort_utils(_sciminer_util_provider_list)


# 动态生成服务类字典
service_classes, sciminer_util_manager_dict, sciminer_util_provider_list = generate_service_dict()
