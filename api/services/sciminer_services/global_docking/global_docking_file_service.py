"""
@Time    : 2024/10/15 上午9:58
@Author  : bigboss
@Description : 全局对接的文件服务
"""
import io
import logging
import uuid
from datetime import datetime
from enum import Enum
from typing import Union, Tuple

import click
from rdkit import Chem
from rdkit.Chem import AllChem
from werkzeug.datastructures import FileStorage

from configs import dify_config
from core.tools.utils.upload_file_utils import UploadFileUtils
from extensions.ext_storage import storage
from models.account import Account
from models.model import EndUser, UploadFile
from services.errors.file import UnsupportedFileTypeError

ALLOWED_EXTENSIONS = ['fasta']
ALLOWED_LIGAND_EXTENSIONS = ['txt', 'smi']


class GlobalDockingSourceType(Enum):
    FASTA = 'fasta'
    LIGAND = 'ligand'


class GlobalDockingFileInfo:
    def __init__(self, file_name, file_content, file_size, mimetype, extension,
                 file_key, current_tenant_id):
        self.file_name = file_name
        self.file_content = file_content
        self.file_size = file_size
        self.mimetype = mimetype
        self.extension = extension
        self.file_key = file_key
        self.current_tenant_id = current_tenant_id


class GlobalDockingFileService:

    @classmethod
    def upload_file(cls, file: Union[str, FileStorage], source: GlobalDockingSourceType,
                    user: Union[Account, EndUser]) -> list[UploadFile]:
        if source == GlobalDockingSourceType.FASTA:
            logging.info(click.style(f"接收到的是fasta字符串，进行文件保存...", fg="green"))
            return cls.save_fasta_file(file, user=user)
        elif source == GlobalDockingSourceType.LIGAND:
            logging.info(click.style(f"接收到的是ligand字符串，进行文件保存...", fg="green"))
            return cls.save_ligand_file(file, user=user)

    @classmethod
    def save_fasta_file(cls, file: Union[str, FileStorage],
                        user: Union[Account, EndUser]) -> list[UploadFile]:
        if isinstance(file, str):
            file_name = f"{str(uuid.uuid4())}.fasta"
            file_content = file.encode('utf-8')
            file_size = len(file_content)
            mimetype = "text/plain"
        else:
            file_name = file.filename
            file_content = file.read()
            file_size = len(file_content)
            mimetype = file.mimetype

        extension = file_name.split('.')[-1]
        if extension.lower() not in ALLOWED_EXTENSIONS:
            raise UnsupportedFileTypeError()

        file_key, current_tenant_id = cls.get_global_docking_file_path(extension, user)

        storage.save(file_key, file_content)
        # save file to database
        upload_file = UploadFileUtils.add_upload_file(
            tenant_id=current_tenant_id,
            storage_type=dify_config.STORAGE_TYPE,
            key=file_key,
            name=file_name,
            size=file_size,
            extension=extension,
            mime_type=mimetype,
            created_by=user.id,
        )

        return [upload_file, ]

    @classmethod
    def save_ligand_file(cls, file: Union[str, FileStorage],
                         user: Union[Account, EndUser]) -> list[UploadFile]:
        """
        保存 ligand 文件
        先判断传递过来的file参数是字符串类型的还是文件流类型的
        字符串类型默认只支持一个smiles字符串，使用rdkit转换mol字符串再保存文件，再保存到数据库
        文件流的话，需要先判断文件格式是否正确，然后将文件内容分割成多个smiles字符串，使用rdkit转换mol字符串再保存文件，再保存到数据库
        :param file: 文件对象
        :param file_name: 文件名
        :param user: 用户
        :return: UploadFile 对象
        """
        # 先校验ligand文件格式是否正确
        global_docking_file_info_list = []
        mimetype = "text/plain"
        if isinstance(file, str):
            file_content_list = file.split('\n')
            file_name = f"{str(uuid.uuid4())}.sdf"
            # 使用 io.StringIO 作为内存中的文件
            sdf_block = io.StringIO()
            # 创建一个 SDWriter 用于写入内存中的文件
            writer = Chem.SDWriter(sdf_block)
            for file_content_temp in file_content_list:
                mol = Chem.MolFromSmiles(file_content_temp.strip())
                if mol is None:
                    logging.error(click.style(f"无法解析ligand字符串，请检查格式是否正确", fg="red"))
                    writer.close()
                    sdf_block.close()
                    raise UnsupportedFileTypeError()

                # 计算 2D 坐标，适用于可视化或 SDF 格式
                AllChem.Compute2DCoords(mol)
                writer.write(mol)
            writer.close()
            file_content = sdf_block.getvalue().encode('utf-8')
            file_size = len(file_content)
            sdf_block.close()

            extension = file_name.split('.')[-1]

            file_key, current_tenant_id = cls.get_global_docking_file_path(extension, user)

            file_info = GlobalDockingFileInfo(
                file_name=file_name,
                file_content=file_content,
                file_size=file_size,
                mimetype=mimetype,
                extension=extension,
                file_key=file_key,
                current_tenant_id=current_tenant_id,
            )
            global_docking_file_info_list.append(file_info)
        else:
            if file.filename.split('.')[-1] not in ALLOWED_LIGAND_EXTENSIONS:
                raise UnsupportedFileTypeError()
            file_content_list = file.read().decode('utf-8').split('\n')

            file_name = f"{str(uuid.uuid4())}.sdf"
            # 使用 io.StringIO 作为内存中的文件
            sdf_block = io.StringIO()
            # 创建一个 SDWriter 用于写入内存中的文件
            writer = Chem.SDWriter(sdf_block)

            for file_content in file_content_list:
                mol = Chem.MolFromSmiles(file_content.strip())
                if mol is None:
                    logging.error(click.style(f"无法解析ligand字符串，请检查格式是否正确", fg="red"))
                    writer.close()
                    sdf_block.close()
                    raise UnsupportedFileTypeError()

                # 计算 2D 坐标，适用于可视化或 SDF 格式
                AllChem.Compute2DCoords(mol)
                writer.write(mol)
            writer.close()
            file_content = sdf_block.getvalue().encode('utf-8')
            file_size = len(file_content)
            sdf_block.close()

            extension = file_name.split('.')[-1]

            file_key, current_tenant_id = cls.get_global_docking_file_path(extension, user)

            file_info = GlobalDockingFileInfo(
                file_name=file_name,
                file_content=file_content,
                file_size=file_size,
                mimetype=mimetype,
                extension=extension,
                file_key=file_key,
                current_tenant_id=current_tenant_id,
            )
            global_docking_file_info_list.append(file_info)

        upload_file_list = []
        # 保存文件
        for file_info in global_docking_file_info_list:
            storage.save(file_info.file_key, file_info.file_content)
            # save file to database
            upload_file = UploadFileUtils.add_upload_file(
                tenant_id=file_info.current_tenant_id,
                storage_type=dify_config.STORAGE_TYPE,
                key=file_info.file_key,
                name=file_info.file_name,
                size=file_info.file_size,
                extension=file_info.extension,
                mime_type=mimetype,
                created_by=user.id,
            )

            upload_file_list.append(upload_file)

        return upload_file_list

    @classmethod
    def get_global_docking_file_path(cls, extension: str, user: Union[Account, EndUser]) -> Tuple[str, str]:
        """
        获取全局对接文件保存路径
        :param extension: 文件后缀名
        :param user:     用户
        :return: 文件路径
        """
        date_path = datetime.now().strftime('%Y/%m/%d')
        if isinstance(user, Account):
            current_tenant_id = user.current_tenant_id
        else:
            current_tenant_id = user.tenant_id
        return f"upload_files/{current_tenant_id}/global_docking/{date_path}/{str(uuid.uuid4())}.{extension}", current_tenant_id