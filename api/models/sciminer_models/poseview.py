"""
@Time    : 2024/11/4 下午2:38
@Author  : bigboss
@Description : 
"""
from extensions.ext_database import db
from models import StringUUID


class PoseViewTask(db.Model):
    """pose view Task result/status."""

    __tablename__ = "utility_poseview_task"

    __table_args__ = (
        # 主键
        db.PrimaryKeyConstraint("id", name="poseview_task_pkey"),
        {'comment': 'pose view 任务结果数据表'}
    )

    id = db.Column(StringUUID, server_default=db.text("uuid_generate_v4()"), comment="poseview的任务id")
    receptor_file_id = db.Column(db.String(155), nullable=True, comment="poseview任务所需的receptor文件id")
    ligand_file_ids = db.Column(db.ARRAY(db.String(155)), nullable=True, comment="poseview任务所需的ligand文件id数组，多个ligand文件以数组形式存储")
    result = db.Column(db.Text, nullable=True, comment="poseview任务的结果")

    created_by = db.Column(StringUUID, nullable=False)
    created_at = db.Column(db.DateTime, nullable=False, server_default=db.text("CURRENT_TIMESTAMP(0)"),
                           comment="创建时间")
    updated_at = db.Column(db.DateTime, nullable=False, server_default=db.text("CURRENT_TIMESTAMP(0)"),
                           onupdate=db.text("CURRENT_TIMESTAMP(0)"), comment="更新时间")

    status = None
    task_name = None

    @property
    def serialize(self):
        return {
            "id": self.id,
            "receptor_file_id": self.receptor_file_id,
            "ligand_file_ids": self.ligand_file_ids,
            "result": self.result,
            "status": self.status,
            "task_name": self.task_name,
            "created_by": self.created_by,
            "created_at": self.created_at,
            "updated_at": self.updated_at
        }

    @property
    def result_dict(self):
        return {
            "id": self.id,
            "receptor_file_id": self.receptor_file_id,
            "ligand_file_ids": self.ligand_file_ids,
            "result": self.result,
            "created_by": self.created_by,
            "created_at": self.created_at.strftime("%Y-%m-%d %H:%M:%S"),
            "updated_at": self.updated_at.strftime("%Y-%m-%d %H:%M:%S")
        }