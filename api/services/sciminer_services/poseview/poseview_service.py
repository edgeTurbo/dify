"""
@Time    : 2024/11/4 上午11:02
@Author  : bigboss
@Description : poseview service
"""
from typing import Union

from controllers.console.sciminer_apps.poseview.request_fields import PoseviewTaskRequestFields
from models.account import Account
from models.model import EndUser
from services.sciminer_services import SciminerBaseService


class PoseViewService(SciminerBaseService):
    service_type = "POSEVIEW"
    task_label = "PoseView"

    @classmethod
    def start_task(cls, poseview_request_fields: PoseviewTaskRequestFields, user: Union[Account, EndUser],
                   start_celery: bool = False):
        pass

    @classmethod
    def get_service_result_data(cls, task_id: str, user: Union[Account, EndUser]):
        pass
