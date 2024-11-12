"""
@Time    : 2024/11/6 上午8:48
@Author  : bigboss
@Description : 
"""
from io import BytesIO, StringIO
from typing import Optional

from rdkit import Chem
import traceback


def read_sdf_bytes(sdf_bytes: bytes):
    """
    Read sdf bytes and return a list of molecules.
    :param sdf_bytes: sdf bytes
    :return: a list of molecules
    """
    try:
        # 使用 SDMolSupplier 读取分子
        mol_list = Chem.ForwardSDMolSupplier(BytesIO(sdf_bytes))
        return mol_list
    except Exception as e:
        return None


def mol_object_to_smiles(mol: Chem.Mol) -> str | None:
    """
    mol 对象转化为 smiles 字符串
    :param mol: mol 对象
    :return: smiles 字符串
    """
    try:
        return Chem.MolToSmiles(mol)
    except Exception as e:
        return None


def mol_object_to_sdf_bytes(mol: Chem.Mol) -> bytes | None:
    """
    mol 对象转化为 sdf 字节流
    :param mol: mol 对象
    :return: sdf 字节流
    """
    try:
        sdf_buffer = StringIO()
        writer = Chem.SDWriter(sdf_buffer)
        writer.write(mol)
        writer.close()
        return sdf_buffer.getvalue().encode('utf-8')
    except Exception as e:
        return None
