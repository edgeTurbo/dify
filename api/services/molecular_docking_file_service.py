import datetime
import hashlib
import io
import uuid
from collections.abc import Generator
from typing import Union

from flask_login import current_user
from flask import send_file, Response
from werkzeug.datastructures import FileStorage
from werkzeug.exceptions import NotFound

from configs import dify_config
from core.file.upload_file_parser import UploadFileParser
from core.rag.extractor.extract_processor import ExtractProcessor
from extensions.ext_database import db
from extensions.ext_storage import storage
from models.account import Account
from models.model import EndUser, UploadFile
from services.errors.file import FileTooLargeError, UnsupportedFileTypeError
from services.file_service import UNSTRUCTURED_ALLOWED_EXTENSIONS, IMAGE_EXTENSIONS
from services.molecular_docking.tool_rendering_2d_structure_service import ToolRendering2DStructureService

ALLOWED_EXTENSIONS = ["pdb", "cif", "bcif", "mmcif", "mol", "sdf", "mol2", "xyz"]

PREVIEW_WORDS_LIMIT = 3000


class MolecularDockingFileService:

    @staticmethod
    def get_pocket_docking_file_path(file_name: str, user: Union[Account, EndUser]):
        """
        获取口袋对接文件保存路径
        :param file_name: 文件名
        :param user:     用户
        :return: 文件路径
        """
        date_path = datetime.datetime.now().strftime('%Y/%m/%d')
        if isinstance(user, Account):
            current_tenant_id = user.current_tenant_id
        else:
            current_tenant_id = user.tenant_id
        return f"upload_files/{current_tenant_id}/pocket_docking/{date_path}/{file_name}", current_tenant_id

    @staticmethod
    def upload_file(file: FileStorage, user: Union[Account, EndUser]) -> UploadFile:
        filename = file.filename
        extension = file.filename.split(".")[-1]
        if len(filename) > 200:
            filename = filename.split(".")[0][:200] + "." + extension
        if extension.lower() not in ALLOWED_EXTENSIONS:
            raise UnsupportedFileTypeError()
        else:
            file_size_limit = dify_config.UPLOAD_FILE_SIZE_LIMIT * 1024 * 1024

        # read file content
        file_content = file.read()

        # get file size
        file_size = len(file_content)

        if file_size > file_size_limit:
            message = f"File size exceeded. {file_size} > {file_size_limit}"
            raise FileTooLargeError(message)

        # user uuid as file name
        file_uuid = str(uuid.uuid4())

        file_key, current_tenant_id = MolecularDockingFileService.get_pocket_docking_file_path(
            file_name=file_uuid + "." + extension,
            user=user
        )

        # save file to storage
        storage.save(file_key, file_content)

        # save file to db
        upload_file = UploadFile(
            tenant_id=current_tenant_id,
            storage_type=dify_config.STORAGE_TYPE,
            key=file_key,
            name=filename,
            size=file_size,
            extension=extension,
            mime_type=file.mimetype,
            created_by_role=("account" if isinstance(user, Account) else "end_user"),
            created_by=user.id,
            created_at=datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None),
            used=False,
            hash=hashlib.sha3_256(file_content).hexdigest(),
        )

        db.session.add(upload_file)
        db.session.commit()

        return upload_file

    @staticmethod
    def get_file(file_id: str, mime_type: str) -> Response:
        upload_file = db.session.query(UploadFile).filter_by(id=file_id).first()
        if upload_file is None:
            raise NotFound("File not found")
        return send_file(io.BytesIO(storage.load_once(upload_file.key)), mimetype=mime_type, as_attachment=True,
                         download_name=upload_file.name)

    @staticmethod
    def upload_text(text: str, text_name: str) -> UploadFile:
        if len(text_name) > 200:
            text_name = text_name[:200]
        # user uuid as file name
        file_uuid = str(uuid.uuid4())
        file_key = "upload_files/" + current_user.current_tenant_id + "/" + file_uuid + ".txt"

        # save file to storage
        storage.save(file_key, text.encode("utf-8"))

        # save file to db
        upload_file = UploadFile(
            tenant_id=current_user.current_tenant_id,
            storage_type=dify_config.STORAGE_TYPE,
            key=file_key,
            name=text_name,
            size=len(text),
            extension="txt",
            mime_type="text/plain",
            created_by=current_user.id,
            created_at=datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None),
            used=True,
            used_by=current_user.id,
            used_at=datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None),
        )

        db.session.add(upload_file)
        db.session.commit()

        return upload_file

    @staticmethod
    def get_file_preview(file_id: str) -> str:
        upload_file = db.session.query(UploadFile).filter(UploadFile.id == file_id).first()

        if not upload_file:
            raise NotFound("File not found")

        # extract text from file
        extension = upload_file.extension
        etl_type = dify_config.ETL_TYPE
        allowed_extensions = UNSTRUCTURED_ALLOWED_EXTENSIONS if etl_type == "Unstructured" else ALLOWED_EXTENSIONS
        if extension.lower() not in allowed_extensions:
            raise UnsupportedFileTypeError()

        text = ExtractProcessor.load_from_upload_file(upload_file, return_text=True)
        text = text[0:PREVIEW_WORDS_LIMIT] if text else ""

        return text

    @staticmethod
    def get_image_preview(file_id: str, timestamp: str, nonce: str, sign: str) -> tuple[Generator, str]:
        result = UploadFileParser.verify_image_file_signature(file_id, timestamp, nonce, sign)
        if not result:
            raise NotFound("File not found or signature is invalid")

        upload_file = db.session.query(UploadFile).filter(UploadFile.id == file_id).first()

        if not upload_file:
            raise NotFound("File not found or signature is invalid")

        # extract text from file
        extension = upload_file.extension
        if extension.lower() not in IMAGE_EXTENSIONS:
            raise UnsupportedFileTypeError()

        generator = storage.load(upload_file.key, stream=True)

        return generator, upload_file.mime_type

    @staticmethod
    def get_public_image_preview(file_id: str) -> tuple[Generator, str]:
        upload_file = db.session.query(UploadFile).filter(UploadFile.id == file_id).first()

        if not upload_file:
            raise NotFound("File not found or signature is invalid")

        # extract text from file
        extension = upload_file.extension
        if extension.lower() not in IMAGE_EXTENSIONS:
            raise UnsupportedFileTypeError()

        generator = storage.load(upload_file.key)

        return generator, upload_file.mime_type

    @staticmethod
    def rendering_molecule_file(file_id: str, user: Union[Account, EndUser]) -> list[str | None]:
        """
        渲染分子结构文件
        :param file_id: 文件id
        :param user: 用户信息
        :return: base64编码的图片列表
        """
        upload_file = db.session.query(UploadFile).filter_by(id=file_id, created_by=user.id).first()
        if upload_file is None:
            raise NotFound("File not found")
        if upload_file.extension.lower() not in ["sdf", "mol"]:
            raise UnsupportedFileTypeError()
        file_bytes = storage.load(upload_file.key)
        image_base64_str_list = ToolRendering2DStructureService.generate_molecule_images_by_bytes(file_bytes,
                                                                                                  upload_file.name)
        return image_base64_str_list
