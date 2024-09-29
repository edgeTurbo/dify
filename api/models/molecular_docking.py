from enum import Enum

from extensions.ext_database import db
from models import StringUUID


class Status(Enum):
    """
    任务状态枚举类
    """
    PENDING = ("PENDING", "等待中")
    PROCESSING = ("PROCESSING", "处理中")
    SUCCESS = ("SUCCESS", "成功")
    FAILURE = ("FAILURE", "失败")

    @property
    def desc(self):
        return self.value[1]

    @property
    def status(self):
        return self.value[0]


class MolecularDockingTask(db.Model):
    """Molecular Docking Task result/status."""

    __tablename__ = "molecular_docking_task"

    __table_args__ = (
        # 主键
        db.PrimaryKeyConstraint("id", name="molecular_docking_task_pkey"),
    )

    id = db.Column(StringUUID, server_default=db.text("uuid_generate_v4()"), comment="分子对接的任务id")

    task_name = db.Column(db.String(155), nullable=True, comment="分子对接任务名称")
    pdb_file_id = db.Column(db.String(155), nullable=True, comment="分子对接任务所需的pdb文件id")

    center_x = db.Column(db.Float, nullable=True, comment="分子对接任务的中心点x坐标")
    center_y = db.Column(db.Float, nullable=True, comment="分子对接任务的中心点y坐标")
    center_z = db.Column(db.Float, nullable=True, comment="分子对接任务的中心点z坐标")

    size_x = db.Column(db.Float, nullable=True, comment="分子对接任务的大小x坐标")
    size_y = db.Column(db.Float, nullable=True, comment="分子对接任务的大小y坐标")
    size_z = db.Column(db.Float, nullable=True, comment="分子对接任务的大小z坐标")

    ligand_file_ids = db.Column(db.ARRAY(db.String(155)), nullable=True, comment="分子对接任务所需的ligand文件id数组，多个ligand文件以数组形式存储")

    out_pose_num = db.Column(db.Integer, nullable=True, comment="分子对接任务的输出结果pose数量")

    result = db.Column(db.Text, nullable=True, comment="分子对接的结果")
    status = db.Column(db.String(50), default=Status.PENDING.status, comment="分子对接的任务状态,默认等待中 PENDING，可选值：处理中 PROCESSING、成功 SUCCESS、失败 FAILURE")

    created_by = db.Column(StringUUID, nullable=False)
    created_at = db.Column(db.DateTime, nullable=False, server_default=db.text("CURRENT_TIMESTAMP(0)"),
                           comment="创建时间")
    updated_at = db.Column(db.DateTime, nullable=False, server_default=db.text("CURRENT_TIMESTAMP(0)"),
                           comment="更新时间")

    @property
    def serialize(self):
        return {
            "id": self.id,
            "task_name": self.task_name,
            "pdb_file_id": self.pdb_file_id,
            "center_x": self.center_x,
            "center_y": self.center_y,
            "center_z": self.center_z,
            "size_x": self.size_x,
            "size_y": self.size_y,
            "size_z": self.size_z,
            "ligand_file_ids": self.ligand_file_ids,
            "out_pose_num": self.out_pose_num,
            "result": self.result,
            "status": self.status,
            "created_by": self.created_by,
            "created_at": self.created_at.strftime("%Y-%m-%d %H:%M:%S"),
            "updated_at": self.updated_at.strftime("%Y-%m-%d %H:%M:%S")
        }
