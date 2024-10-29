from flask import request
from werkzeug.exceptions import Unauthorized

from libs.passport import PassportService

client_query = {}


def get_client_id(uid):
    """获取特定 UID 的 client_id"""
    return client_query.get(str(uid))


def send_message_to_user(uid, message: dict):
    """向指定 UID 的客户端发送消息"""
    send_users = get_client_id(uid)
    if send_users is not None and len(send_users) != 0:
        try:
            for client in send_users:
                client.send(message)
            return True
        except Exception as e:
            return False
    else:
        print(f'Client with UID {uid} not found')
        return False


def websocket_connection(ws):
    auth_header = request.args.get("token")

    current_user_id = get_tenant_by_header(auth_header)
    if current_user_id is not None:
        token_list = client_query.get(current_user_id, [])
        token_list.append(ws)
        client_query[current_user_id] = token_list
    try:
        while True:
            message = ws.receive()
            if message:
                if message == '{"ping":""}':
                    ws.send("pong")
            else:
                break
    finally:
        for key, value in client_query.items():
            for index, item_list in enumerate(value):
                if item_list == ws:
                    client_query[key].pop(index)


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
    if user_id is not None:
        return user_id
    else:
        return None
