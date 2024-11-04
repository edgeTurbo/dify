"""
@Time    : 2024/11/1 下午3:41
@Author  : bigboss
@Description : 针对数据表upload_files的各个操作工具类
"""
from typing import Union

from extensions.ext_database import db
from extensions.ext_storage import storage
from models.model import UploadFile
from services.errors.file import UnsupportedFileTypeError


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
