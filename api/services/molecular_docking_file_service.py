import datetime
import hashlib
import io
import uuid
from typing import Union

from flask import send_file, Response
from werkzeug.datastructures import FileStorage
from werkzeug.exceptions import NotFound

from configs import dify_config
from core.tools.utils.upload_file_utils import UploadFileUtils
from extensions.ext_database import db
from extensions.ext_storage import storage
from models.account import Account
from models.model import EndUser, UploadFile
from services.errors.file import FileTooLargeError, UnsupportedFileTypeError
from services.sciminer_services.molecular_docking.tool_rendering_2d_structure_service import \
    ToolRendering2DStructureService

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

        # read file content
        file_content = file.read()

        # get file size
        file_size = len(file_content)

        # user uuid as file name
        file_uuid = str(uuid.uuid4())

        file_key, current_tenant_id = MolecularDockingFileService.get_pocket_docking_file_path(
            file_name=file_uuid + "." + extension,
            user=user
        )

        UploadFileUtils.upload_file_to_storage(
            file_key=file_key,
            file_content=file_content,
            extension=extension,
            allowed_extensions=ALLOWED_EXTENSIONS
        )

        # save file to db
        upload_file = UploadFileUtils.add_upload_file_to_db(
            tenant_id=current_tenant_id,
            storage_type=dify_config.STORAGE_TYPE,
            key=file_key,
            name=filename,
            size=file_size,
            extension=extension,
            mime_type=file.mimetype,
            created_by_role=("account" if isinstance(user, Account) else "end_user"),
            created_by=user.id,
            used=False,
            hash=hashlib.sha3_256(file_content).hexdigest(),
        )

        return upload_file

    @staticmethod
    def get_file(file_id: str, mime_type: str) -> Response:
        upload_file = db.session.query(UploadFile).filter_by(id=file_id).first()
        if upload_file is None:
            raise NotFound("File not found")
        return send_file(io.BytesIO(storage.load_once(upload_file.key)), mimetype=mime_type, as_attachment=True,
                         download_name=upload_file.name)

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
        if upload_file.extension.lower() not in ["sdf", "mol", "pdb", "mol2"]:
            raise UnsupportedFileTypeError()
        file_bytes = storage.load(upload_file.key)
        image_base64_str_list = ToolRendering2DStructureService.generate_molecule_images_by_bytes(file_bytes,
                                                                                                  upload_file.name)
        return image_base64_str_list
