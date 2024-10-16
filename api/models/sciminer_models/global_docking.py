"""
@Time    : 2024/10/16 下午1:36
@Author  : bigboss
@Description : 全局分子对接的模型类
"""
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


class GlobalDockingTask(db.Model):
    """Global Docking Task result/status."""

    __tablename__ = "global_docking_task"

    __table_args__ = (
        # 主键
        db.PrimaryKeyConstraint("id", name="global_docking_task_pkey"),
        {'comment': '全局分子对接任务结果数据表'}
    )

    id = db.Column(StringUUID, server_default=db.text("uuid_generate_v4()"), comment="全局分子对接的任务id")

    fasta_file_id = db.Column(db.String(155), nullable=True, comment="全局分子对接任务所需的fasta文件id")

    ligand_file_ids = db.Column(db.ARRAY(db.String(155)), nullable=True, comment="全局分子对接任务所需的ligand文件id数组，多个ligand文件以数组形式存储")

    out_pose_num = db.Column(db.Integer, nullable=True, comment="全局分子对接任务的输出结果pose数量")

    result = db.Column(db.Text, nullable=True, comment="全局分子对接的结果")

    created_by = db.Column(StringUUID, nullable=False)
    created_at = db.Column(db.DateTime, nullable=False, server_default=db.text("CURRENT_TIMESTAMP(0)"),
                           comment="创建时间")
    updated_at = db.Column(db.DateTime, nullable=False, server_default=db.text("CURRENT_TIMESTAMP(0)"),
                           onupdate=db.text("CURRENT_TIMESTAMP(0)"), comment="更新时间")

    status = None

    @property
    def serialize(self):
        return {
            "id": self.id,
            "fasta_file_id": self.fasta_file_id,
            "ligand_file_id": self.ligand_file_id,
            "out_pose_num": self.out_pose_num,
            "result": self.result,
            "created_by": self.created_by,
            "created_at": self.created_at.strftime("%Y-%m-%d %H:%M:%S"),
            "updated_at": self.updated_at.strftime("%Y-%m-%d %H:%M:%S")
        }
