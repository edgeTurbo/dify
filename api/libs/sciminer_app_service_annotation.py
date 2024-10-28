"""
@Time    : 2024/10/28 上午9:29
@Author  : bigboss
@Description : 
"""


def sciminer_app_service(service_type: str):
    """
    注解装饰器，用于标记Sciminer应用服务类，可以通过反射来达到动态调用方法的目的
    service_type需要和sciminer_history_task数据表中的task_type字段对应，才能达到该效果
    """

    def decorator(cls):
        cls.service_type = service_type
        return cls

    return decorator
