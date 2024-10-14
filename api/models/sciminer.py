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


class SciminerHistoryTask(db.Model):
    """
    sciminer的工具调用任务历史记录
    """
    __tablename__ = "sciminer_history_task"

    __table_args__ = (
        # 主键
        db.PrimaryKeyConstraint("id", name="sciminer_history_task_pkey"),
    )

    id = db.Column(StringUUID, server_default=db.text("uuid_generate_v4()"), comment="历史记录数据表id")

    task_id = db.Column(StringUUID, nullable=False, comment="任务id", index=True)
    task_name = db.Column(db.String(155), nullable=True, comment="任务名称", index=True)

    task_type = db.Column(db.String(155), nullable=False, comment="任务类型", index=True)
    label = db.Column(db.String(155), nullable=False, comment="任务类别标签")

    status = db.Column(db.String(50), default=Status.PENDING.status,
                       comment="任务状态,默认等待中 PENDING，可选值：处理中 PROCESSING、成功 SUCCESS、失败 FAILURE")

    created_by = db.Column(StringUUID, nullable=False, comment="创建人id", index=True)
    created_at = db.Column(db.DateTime, nullable=False, server_default=db.text("CURRENT_TIMESTAMP(0)"),
                           comment="创建时间")
    updated_at = db.Column(db.DateTime, nullable=False, server_default=db.text("CURRENT_TIMESTAMP(0)"),
                           onupdate=db.text("CURRENT_TIMESTAMP(0)"), comment="更新时间")