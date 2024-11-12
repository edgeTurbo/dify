"""
@Time    : 2024/11/1 下午3:41
@Author  : bigboss
@Description : 针对数据表result_files的各个操作工具类
"""
from io import BufferedReader
from typing import Union

from core.tools.utils.upload_file_utils import UploadFileProperty
from extensions.ext_database import db
from extensions.ext_storage import storage
from models.account import Account
from models.model import EndUser
from models.sciminer_models.sciminer import ResultFile
from services.errors.file import UnsupportedFileTypeError


class ResultFileProperty(UploadFileProperty):
    def __init__(self, file_name: str, file_bytes: bytes, extension: str):
        super().__init__(file_name, file_bytes, extension)


class ResultFileUtils(ResultFile):
    @classmethod
    def get_result_file_content_by_id(cls, result_file_id: str, user_id: str = None) -> str:
        """
        根据result_file_id获取result_file的content文本内容，用于提取txt,csv等文件内容
        :param result_file_id:结果文件id
        :param user_id: 用户id
        :return: str
        """
        if user_id is None:
            result_file = db.session.query(ResultFile).filter(
                ResultFile.id == result_file_id,
            ).first()
        else:
            result_file = db.session.query(ResultFile).filter(
                ResultFile.id == result_file_id,
                ResultFile.created_by == user_id
            ).first()
        if result_file is None:
            raise ValueError("文件查找不到")

        try:
            if isinstance(result_file, ResultFile):
                file_bytes = storage.load_once(result_file.key)
                return file_bytes.decode('utf-8')
        except Exception as e:
            raise ValueError("文件读取失败")

    @classmethod
    def add_result_file_to_db(cls, **data) -> ResultFile:
        """
        新增上传文件到result_files表中
        :param data: 结果文件数据
        :return: ResultFiles
        """
        result_file = ResultFile(**data)
        db.session.add(result_file)
        db.session.commit()
        return result_file

    @classmethod
    def result_file_to_storage(cls, file_key: str, file_content,
                               extension: str = None,
                               allowed_extensions: list[str] = None) -> None:
        """
        将结果文件上传到云存储中
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
    def get_result_file_buffer(cls, result_file_ids: Union[str, list[str]], user: Union[Account, EndUser],
                               return_list: bool = False) -> Union[BufferedReader, list[BufferedReader]]:
        """
        根据result_file_id获取result_file的BufferedReader内容，用于提取图片等二进制文件内容
        :param result_file_ids: 结果文件id
        :param user: 用户
        :param return_list: 是否返回list
        :return: BufferedReader or list of BufferedReader
        """
        if return_list:
            result_file_buffer_list = []
            for result_file_id in result_file_ids:
                result_file = ResultFile.query.filter_by(id=result_file_id, created_by=user.id).first()
                result_file_buffer_list.append(storage.load_buffer(result_file.key, result_file.name))
            return result_file_buffer_list
        else:
            result_file = ResultFile.query.filter_by(id=result_file_ids, created_by=user.id).first()
            return storage.load_buffer(result_file.key, result_file.name)

    @classmethod
    def get_result_file_bytes(cls, result_file_ids: Union[str, list[str]], user: Union[Account, EndUser],
                              return_list: bool = False) -> Union[ResultFileProperty, list[ResultFileProperty]]:
        """
        根据result_file_id获取result_file的bytes内容，用于提取图片等二进制文件内容
        :param result_file_ids: 结果文件id
        :param user: 用户
        :param return_list: 是否返回list
        :return: bytes or list of bytes
        """
        if return_list:
            result_file_property_list = []
            for result_file_id in result_file_ids:
                result_file = ResultFile.query.filter_by(id=result_file_id, created_by=user.id).first()
                result_file_property = ResultFileProperty(result_file.name, storage.load_once(result_file.key), result_file.extension)
                result_file_property_list.append(result_file_property)
            return result_file_property_list
        else:
            result_file = ResultFile.query.filter_by(id=result_file_ids, created_by=user.id).first()
            result_file_property = ResultFileProperty(result_file.name, storage.load_once(result_file.key),
                                                      result_file.extension)
            return result_file_property
