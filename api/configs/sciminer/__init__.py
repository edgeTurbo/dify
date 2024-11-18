"""
@Time    : 2024/11/14 上午10:53
@Author  : bigboss
@Description : 
"""
from pydantic import Field
from pydantic_settings import BaseSettings


class MolecularDockingAPIConfig(BaseSettings):
    MOLECULAR_DOCKING_API_URL: str = Field(
        description="Pocket docking API URL",
        default=None,
    )


class GlobalDockingAPIConfig(BaseSettings):
    GLOBAL_DOCKING_API_URL: str = Field(
        description="Global docking API URL",
        default=None,
    )
    GLOBAL_DOCKING_CALLBACK_API_URL: str = Field(
        description="Global docking callback API URL,这是用来提供给全局对接模型的回调接口，用来更新任务状态和结果",
        default=None,
    )



class PoseViewAPIConfig(BaseSettings):
    POSEVIEW_TASK_API_URL: str = Field(
        description="PoseView task API URL",
        default=None,
    )
    POSEVIEW_RESULT_API_URL: str = Field(
        description="PoseView result API URL",
        default=None,
    )


class MolcraftAPIConfig(BaseSettings):
    MOLCRAFT_API_URL: str = Field(
        description="Molcraft API URL",
        default=None,
    )


class SciminerApiConfig(
    MolecularDockingAPIConfig,
    GlobalDockingAPIConfig,
    PoseViewAPIConfig,
    MolcraftAPIConfig,
):
    pass
