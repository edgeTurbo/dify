"""
@Time    : 2024/11/1 下午3:41
@Author  : bigboss
@Description : 针对数据表upload_files的各个操作工具类
"""
from io import BufferedReader
from typing import Union

from extensions.ext_database import db
from extensions.ext_storage import storage
from models.account import Account
from models.model import UploadFile, EndUser
from services.errors.file import UnsupportedFileTypeError


class UploadFileProperty:
    def __init__(self, file_name: str, file_bytes: bytes, extension: str):
        self.file_name = file_name
        self.file_bytes = file_bytes
        self.extension = extension


class UploadFileUtils(UploadFile):
    @classmethod
    def get_upload_file_content_by_id(cls, upload_file_id: str, user_id: str = None) -> str:
        """
        根据upload_file_id获取upload_file的content文本内容，用于提取txt,csv等文件内容
        :param upload_file_id:上传文件id
        :param user_id: 用户id
        :return: str
        """
        if user_id is None:
            upload_file = db.session.query(UploadFile).filter(
                UploadFile.id == upload_file_id,
            ).first()
        else:
            upload_file = db.session.query(UploadFile).filter(
                UploadFile.id == upload_file_id,
                UploadFile.created_by == user_id
            ).first()
        if upload_file is None:
            raise ValueError("文件查找不到")

        try:
            if isinstance(upload_file, UploadFile):
                file_bytes = storage.load_once(upload_file.key)
                return file_bytes.decode('utf-8')
        except Exception as e:
            raise ValueError("文件读取失败")

    @classmethod
    def add_upload_file_to_db(cls, **data) -> UploadFile:
        """
        新增上传文件到upload_files表中
        :param data: 上传文件数据
        :return: UploadFile
        """
        upload_file = UploadFile(**data)
        db.session.add(upload_file)
        db.session.commit()
        return upload_file

    @classmethod
    def upload_file_to_storage(cls, file_key: str, file_content,
                               extension: str = None,
                               allowed_extensions: list[str] = None) -> None:
        """
        将上传文件上传到云存储中
        :param extension: 文件后缀名
        :param file_key: 上传文件key路径
        :param file_content: 上传文件数据
        :param allowed_extensions: 允许的上传文件类型
        :return: None
        """
        if allowed_extensions is not None and extension not in allowed_extensions:
            raise UnsupportedFileTypeError()
        storage.save(file_key, file_content)

    @classmethod
    def get_upload_file_buffer(cls, upload_file_ids: Union[str, list[str]], user: Union[Account, EndUser],
                               return_list: bool = False) -> Union[BufferedReader, list[BufferedReader]]:
        """
        根据upload_file_id获取upload_file的BufferedReader内容，用于提取图片等二进制文件内容
        :param upload_file_ids: 上传文件id
        :param user: 用户
        :param return_list: 是否返回list
        :return: BufferedReader or list of BufferedReader
        """
        if return_list:
            upload_file_buffer_list = []
            for upload_file_id in upload_file_ids:
                upload_file = UploadFile.query.filter_by(id=upload_file_id, created_by=user.id).first()
                upload_file_buffer_list.append(storage.load_buffer(upload_file.key, upload_file.name))
            return upload_file_buffer_list
        else:
            upload_file = UploadFile.query.filter_by(id=upload_file_ids, created_by=user.id).first()
            return storage.load_buffer(upload_file.key, upload_file.name)

    @classmethod
    def get_upload_file_bytes(cls, upload_file_ids: Union[str, list[str]], user: Union[Account, EndUser],
                              return_list: bool = False) -> Union[UploadFileProperty, list[UploadFileProperty]]:
        """
        根据upload_file_id获取upload_file的bytes内容，用于提取图片等二进制文件内容
        :param upload_file_ids: 上传文件id
        :param user: 用户
        :param return_list: 是否返回list
        :return: bytes or list of bytes
        """
        if return_list:
            upload_file_property_list = []
            for upload_file_id in upload_file_ids:
                upload_file = UploadFile.query.filter_by(id=upload_file_id, created_by=user.id).first()
                upload_file_property = UploadFileProperty(upload_file.name, storage.load_once(upload_file.key), upload_file.extension)
                upload_file_property_list.append(upload_file_property)
            return upload_file_property_list
        else:
            upload_file = UploadFile.query.filter_by(id=upload_file_ids, created_by=user.id).first()
            upload_file_property = UploadFileProperty(upload_file.name, storage.load_once(upload_file.key),
                                                      upload_file.extension)
            return upload_file_property
