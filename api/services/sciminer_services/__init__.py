"""
@Time    : 2024/10/25 下午3:45
@Author  : bigboss
@Description : 
"""
import importlib
import pkgutil
from typing import Any, Type
from os import path

from configs import dify_config
from core.tools.utils.yaml_utils import load_yaml_file
from services.sciminer_services.sciminer_base_service import SciminerBaseService


def generate_service_dict() -> (dict[Any, Type[SciminerBaseService]], dict):
    _service_classes = {}
    _service_yaml_list = []
    _icon_dict = {}
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
                    service_yaml_info = load_yaml_file(service_yaml_file, ignore_error=False)
                    if "icon" in service_yaml_info["identity"]:
                        icon_name = service_yaml_info["identity"]["icon"]
                        service_yaml_info["identity"]["icon"] = (f"{dify_config.CONSOLE_API_URL}/console/api/sciminer/"
                                                                 f"utils/{service_yaml_info['identity']['name']}/icon")
                        _icon_dict[service_yaml_info["identity"]["name"]] = path.join(path.dirname(module.__file__),
                                                                                      "_assets", f"{icon_name}")
                    _service_yaml_list.append(service_yaml_info)
    return _service_classes, _service_yaml_list, _icon_dict


# 动态生成服务类字典
service_classes, service_yaml_list, icon_dict = generate_service_dict()
