"""
@Time    : 2024/11/14 上午10:48
@Author  : bigboss
@Description : 
"""
import logging
from typing import Union

import click

from configs import dify_config
from controllers.console.sciminer_apps.molcraft.request_fields import MolcraftTaskRequestFields
from models.account import Account
from models.model import EndUser
from models.sciminer_models.sciminer import Status
from services.sciminer_services import SciminerBaseService

if dify_config.MOLCRAFT_API_URL == "" or dify_config.MOLCRAFT_API_URL is None:
    logging.error(
        click.style(
            "molcraft任务API URL不能为空, 请在配置文件.env中设置MOLCRAFT_API_URL, 否则molcraft任务功能无法正常使用",
            fg='red', bold=True))

if dify_config.INNER_API is None or dify_config.INNER_API_KEY is None:
    logging.error(
        click.style(
            "内网API地址或API KEY未设置, 请在配置文件.env中设置INNER_API为true和INNER_API_KEY, 否则相关sciminer应用功能无法正常使用，websocket消息队列无法正常使用",
            fg='red', bold=True)
    )


class MolcraftService(SciminerBaseService):
    service_type = "MOLCRAFT"
    task_label = "Molcraft"

    @classmethod
    def get_service_result_data(cls, task_id: str, user: Union[Account, EndUser]):
        pass

    @classmethod
    def start_task(cls, molcraft_request_fields: MolcraftTaskRequestFields, user: Union[Account, EndUser], start_celery: bool = False):
        status = Status.PROCESSING.status
