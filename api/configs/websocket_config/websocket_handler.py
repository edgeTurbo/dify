import json
import os

from flask import request
from flask_socketio import SocketIO, emit
from werkzeug.exceptions import Unauthorized

from libs.passport import PassportService
from services.account_service import AccountService, TenantService

client_query = {}
NAME_SPACE = '/websocket'

os.environ["OPENAI_API_KEY"] = "YOUR_OPENAI_API_KEY"


class WebSocketHandler:
    _instance = None

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(WebSocketHandler, cls).__new__(cls)
        return cls._instance

    def __init__(self, socketio: SocketIO):
        if not hasattr(self, 'initialized'):  # 防止重复初始化
            self.socketio = socketio
            self.register_events()
            self.initialized = True

    def register_events(self):

        @self.socketio.on('connect', namespace=NAME_SPACE)
        def on_connect():
            auth_header = request.headers.get("Authorization")
            current_user = get_tenant_by_header(auth_header)
            if current_user:
                # uid = request.args.get('uid')
                client_id = request.sid
                current_user_id = current_user.id
                token_list = client_query.get(current_user_id, [])
                token_list.append(client_id)
                client_query[current_user_id] = token_list
                print(client_query)
                print(f'New connection: id={client_id}, current connections={len(client_query)}')
            else:
                return False

        @self.socketio.on('disconnect', namespace=NAME_SPACE)
        def on_disconnect():
            client_sid = request.sid
            for key, value in client_query.items():
                for index, item_list in enumerate(value):
                    if item_list == client_sid:
                        client_query[key].pop(index)
                        print(client_query)
                        print(f'Connection closed: id={request.sid}, current connections={len(client_query)}')

        @self.socketio.on('message', namespace=NAME_SPACE)
        def on_message(message):
            print(f'Received message from id={request.sid}: {message}')
            emit('my_response_message', "Message received", broadcast=False, namespace=NAME_SPACE, room=request.sid)


# param uid 用户id
def get_client_id(uid):
    """获取特定 UID 的 client_id"""
    return client_query.get(str(uid))


# 封装的发送消息方法，带错误处理
def send_message_to_user(title: str, uid, message: dict):
    """向指定 UID 的客户端发送消息"""
    send_users = get_client_id(uid)
    if send_users is not None and len(send_users) != 0:
        try:
            for send_user in send_users:
                emit(title, json.dumps(message), broadcast=False, namespace=NAME_SPACE, room=send_user)
                print(f'Message sent to client with UID {uid}: {message}')
            return True
        except Exception as e:
            print(f'Error sending message to client with UID {uid}: {e}')
            return False
    else:
        print(f'Client with UID {uid} not found')
        return False


def send_message_to_all(title: str, message: str):
    """向所有的客户端发送消息"""
    try:
        emit(title, json.dumps(message), broadcast=True, namespace=NAME_SPACE)
        return True
    except Exception as e:
        return False


def get_tenant_by_header(auth_header):
    if not auth_header:
        auth_token = request.args.get("_token")
        if not auth_token:
            raise Unauthorized("Invalid Authorization token.")
    else:
        if " " not in auth_header:
            raise Unauthorized("Invalid Authorization header format. Expected 'Bearer <api-key>' format.")
        auth_scheme, auth_token = auth_header.split(None, 1)
        auth_scheme = auth_scheme.lower()
        if auth_scheme != "bearer":
            raise Unauthorized("Invalid Authorization header format. Expected 'Bearer <api-key>' format.")

    decoded = PassportService().verify(auth_token)
    user_id = decoded.get("user_id")
    try:
        account = AccountService.load_logged_in_account(account_id=user_id, token=auth_token)
        # account = AccountService.authenticate('1505938734@qq.com', '12yisan14')
        # tenants = TenantService.get_join_tenants(account)
        # if (len(tenants) == 0):
        #     return None
        return account
    except:
        return None


if __name__ == '__main__':

    import threading


    def job():
        print("任务正在执行...")
        print(client_query)


    def schedule_job():
        job()
        threading.Timer(5, schedule_job).start()


    schedule_job()
