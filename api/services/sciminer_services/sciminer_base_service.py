"""
@Time    : 2024/10/28 下午1:35
@Author  : bigboss
@Description : 
"""
from typing import Union

from models.account import Account
from models.model import EndUser


class SciminerBaseService:
    service_type = None

    def get_service_result_data(self, task_id: str, user: Union[Account, EndUser]):
        """
        获取app服务的结果详细数据
        """
        pass
