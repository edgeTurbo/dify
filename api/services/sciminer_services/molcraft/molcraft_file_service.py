"""
@Time    : 2024/11/4 上午10:23
@Author  : bigboss
@Description : molcraft 文件服务
"""
import hashlib
import uuid
from datetime import datetime
from enum import Enum
from typing import Union

from werkzeug.datastructures import FileStorage

from configs import dify_config
from core.tools.utils.upload_file_utils import UploadFileUtils
from models.account import Account
from models.model import EndUser, UploadFile


class MolcraftSourceType(Enum):
    RECEPTOR = "receptor"
    LIGAND = "ligand"


class MolcraftFileService:
    ALLOWED_RECEPTOR_EXTENSIONS = ["pdb"]

    ALLOWED_LIGAND_EXTENSIONS = ["mol", "sdf", "mol2"]

    @classmethod
    def upload_file(cls, file: FileStorage, user: Union[Account, EndUser], source_type: MolcraftSourceType) -> UploadFile:
        """
        上传文件
        :param file: 文件对象
        :param user: 用户
        :param source_type: 原子来源类型
        :return: 上传文件对象
        """
        file_name = file.filename
        file_content = file.read()
        file_size = len(file_content)
        file_key, current_tenant_id, lower_file_extension = cls.get_molcraft_file_key(
            file_name=file_name,
            user=user
        )
        if source_type == MolcraftSourceType.RECEPTOR:
            UploadFileUtils.upload_file_to_storage(
                file_key=file_key,
                file_content=file_content,
                extension=lower_file_extension,
                allowed_extensions=cls.ALLOWED_RECEPTOR_EXTENSIONS
            )
        elif source_type == MolcraftSourceType.LIGAND:
            UploadFileUtils.upload_file_to_storage(
                file_key=file_key,
                file_content=file_content,
                extension=lower_file_extension,
                allowed_extensions=cls.ALLOWED_LIGAND_EXTENSIONS
            )
        upload_file = UploadFileUtils.add_upload_file_to_db(
            tenant_id=current_tenant_id,
            storage_type=dify_config.STORAGE_TYPE,
            key=file_key,
            name=file_name,
            size=file_size,
            extension=lower_file_extension,
            mime_type=file.mimetype,
            created_by_role=("account" if isinstance(user, Account) else "end_user"),
            created_by=user.id,
            used=False,
            hash=hashlib.sha3_256(file_content).hexdigest(),
        )
        return upload_file

    @classmethod
    def get_molcraft_file_key(cls, file_name: str, user: Union[Account, EndUser]) -> tuple[str, str, str]:
        """
        获取poseview文件保存路径
        """
        date_path = datetime.now().strftime('%Y/%m/%d')
        if isinstance(user, Account):
            current_tenant_id = user.current_tenant_id
        else:
            current_tenant_id = user.tenant_id
        file_extension = file_name.split(".")[-1]
        storage_file_name = f"{str(uuid.uuid4())}.{file_extension}"
        return f"upload_files/{current_tenant_id}/molcraft/{date_path}/{storage_file_name}", current_tenant_id, file_extension.lower()
