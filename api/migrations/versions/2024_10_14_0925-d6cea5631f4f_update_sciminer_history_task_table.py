"""update sciminer_history_task table

Revision ID: d6cea5631f4f
Revises: 1a1a725a318a
Create Date: 2024-10-14 09:25:46.911569

"""
from alembic import op
import models as models
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'd6cea5631f4f'
down_revision = '1a1a725a318a'
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table('sciminer_history_task', schema=None) as batch_op:
        batch_op.add_column(sa.Column('task_type', sa.String(length=155), nullable=False, comment='任务类型'))
        batch_op.drop_index('sciminer_history_task_label_idx')
        batch_op.create_index(batch_op.f('sciminer_history_task_task_type_idx'), ['task_type'], unique=False)

    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table('sciminer_history_task', schema=None) as batch_op:
        batch_op.drop_index(batch_op.f('sciminer_history_task_task_type_idx'))
        batch_op.create_index('sciminer_history_task_label_idx', ['label'], unique=False)
        batch_op.drop_column('task_type')

    # ### end Alembic commands ###