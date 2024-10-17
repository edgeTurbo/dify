from datetime import datetime

from flask_restful import fields


# 自定义时间格式化字段
class DateTimeFormatted(fields.Raw):
    def format(self, value):
        if isinstance(value, datetime):
            return value.strftime('%Y-%m-%d %H:%M:%S')
        raise ValueError('Expected a datetime object')


global_docking_task_fields = {
    'id': fields.String,
    'task_name': fields.String,
    # 定义输出字段并自定义字段名称
    # 'result_dict' 将作为输出字段名用于前端显示获取，原始数据来自 'result'
    'result': fields.String(attribute="result"),
    'status': fields.String,
    'created_at': DateTimeFormatted,
    'updated_at': DateTimeFormatted,
}
