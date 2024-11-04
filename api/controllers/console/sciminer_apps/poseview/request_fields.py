"""
@Time    : 2024/11/4 下午2:02
@Author  : bigboss
@Description : 请求参数对象定义
"""
import uuid

from controllers.console.sciminer_apps.error import IllegalParametersError


class PoseviewTaskRequestFields:
    def __init__(self, task_name: str = None, pdb_file_id: str = None, ligand_file_id: str = None):
        if task_name is None:
            self.task_name = f"poseview_task_{str(uuid.uuid4())}"
        else:
            self.task_name = task_name
        self.pdb_file_id = pdb_file_id
        self.ligand_file_id = ligand_file_id

        if self.pdb_file_id is None or self.ligand_file_id is None:
            raise IllegalParametersError()
