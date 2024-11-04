"""
@Time    : 2024/10/28 下午1:35
@Author  : bigboss
@Description : 
"""
from abc import ABC, abstractmethod
from typing import Union

from models.account import Account
from models.model import EndUser


class SciminerBaseService(ABC):
    service_type = None
    task_label = None

    @abstractmethod
    def get_service_result_data(self, task_id: str, user: Union[Account, EndUser]):
        """
        获取app服务的结果详细数据
        """
        pass
