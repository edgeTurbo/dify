"""
@Time    : 2024/10/24 下午1:38
@Author  : bigboss
@Description : sciminer的返回给前端的字段定义
"""
from datetime import datetime

from flask_restful import fields


# 自定义时间格式化字段
class DateTimeFormatted(fields.Raw):
    def format(self, value):
        if isinstance(value, datetime):
            return value.strftime('%Y-%m-%d %H:%M:%S')
        raise ValueError('Expected a datetime object')


history_task_fields = {
    "id": fields.String,
    "task_id": fields.String,
    "task_name": fields.String,
    "label": fields.String,
    "status": fields.String,
    "created_at": DateTimeFormatted,
    # "updated_at": DateTimeFormatted,
    "task_type": fields.String,
}

history_task_pagination_fields = {
    "page": fields.Integer,
    "limit": fields.Integer(attribute="per_page"),
    "total": fields.Integer,
    "has_more": fields.Boolean(attribute="has_next"),
    "data": fields.List(fields.Nested(history_task_fields), attribute="items"),
}
