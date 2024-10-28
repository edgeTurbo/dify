"""
@Time    : 2024/10/25 下午3:45
@Author  : bigboss
@Description : 
"""
import importlib
import pkgutil
from typing import Any, Type

from services.sciminer_services.sciminer_base_service import SciminerBaseService


def generate_service_dict() -> dict[Any, Type[SciminerBaseService]]:
    _service_classes = {}
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
    return _service_classes


# 动态生成服务类字典
service_classes = generate_service_dict()


