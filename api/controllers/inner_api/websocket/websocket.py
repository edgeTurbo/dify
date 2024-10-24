import json

import requests
from flask import request
from flask_restful import Resource

from configs import dify_config
from configs.websocket_config import websocket_handler
from controllers.console.setup import setup_required
from controllers.inner_api import api
from controllers.inner_api.wraps import inner_api_only


class WebsocketInternalSend(Resource):
    """
    用于内部调用，向指定用户发送websocket消息
    例如分子对接celery任务传递任务状态，因不在同一个服务中，所以需要通过内部调用的方式来实现
    """

    @setup_required
    @inner_api_only
    def post(self):
        user_id = request.form.get("user_id")
        message = request.form.get("message")
        channel = request.form.get("channel")
        # 发送websocket消息通知前端任务状态变化
        websocket_handler.send_message_to_user(
            uid=user_id,
            message={
                channel: json.loads(message)
            },
        )
        return True, 200


def calling_websocket_internal_send(channel: str, user_id: str, message: dict):
    """
    用于外部调用，向指定用户发送websocket消息
    channel: 消息频道，前端监听哪个频道的消息
    user_id: 用户id
    message: 消息内容
    """
    headers = {
        "X-Inner-Api-Key": f"{dify_config.INNER_API_KEY}"
    }
    data = {
        "channel": channel,
        "user_id": user_id,
        "message": json.dumps(message),
    }
    requests.post(f"{dify_config.CONSOLE_API_URL}/inner/api/websocket/send-message", headers=headers, data=data)


api.add_resource(WebsocketInternalSend, "/websocket/send-message")
