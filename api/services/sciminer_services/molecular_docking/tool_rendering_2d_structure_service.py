import io
import os
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import logging
import base64


class ToolRendering2DStructureService:

    @classmethod
    def image_to_base64(cls, img):
        img_byte = io.BytesIO()
        img.save(img_byte, format='PNG')
        img_byte = img_byte.getvalue()
        base64_str = base64.b64encode(img_byte).decode()
        return base64_str

    @classmethod
    # 处理分子并生成2D图像
    def process_molecule(cls, mol, output_folder, idx=1):
        if mol is None:
            return

        AllChem.Compute2DCoords(mol)

        drawer = Draw.MolDraw2DSVG(300, 300)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg

    # 处理SDF文件
    @classmethod
    def process_sdf(cls, sdf_file, output_folder):
        if isinstance(sdf_file, bytes):
            suppl = Chem.ForwardSDMolSupplier(io.BytesIO(sdf_file))
        else:
            suppl = Chem.SDMolSupplier(sdf_file)
        list_of_images = []
        for idx, mol in enumerate(suppl):
            list_of_images.append(cls.process_molecule(mol, output_folder, idx + 1))
        return list_of_images

    # 处理PDB文件
    @classmethod
    def process_pdb(cls, pdb_file, output_folder):
        if isinstance(pdb_file, bytes):
            mol = Chem.MolFromPDBBlock(pdb_file.decode('utf-8'))
        else:
            mol = Chem.MolFromPDBFile(pdb_file)
        if mol is None:
            logging.error(f"无法读取该PDB文件")
            return None
        else:
            return [cls.process_molecule(mol, output_folder), ]

    # 处理MOL文件
    @classmethod
    def process_mol(cls, mol_file, output_folder):
        if isinstance(mol_file, bytes):
            mol = Chem.MolFromMolBlock(mol_file.decode('utf-8'))
        else:
            mol = Chem.MolFromMolFile(mol_file)
        if mol is None:
            logging.error(f"无法读取该MOL文件")
            return None
        else:
            return [cls.process_molecule(mol, output_folder), ]

    # 处理MOL2文件
    @classmethod
    def process_mol2(cls, mol2_file, output_folder):
        if isinstance(mol2_file, bytes):
            mol = Chem.MolFromMol2Block(mol2_file.decode('utf-8'))
        else:
            mol = Chem.MolFromMol2File(mol2_file)
        if mol is None:
            logging.error(f"无法读取该MOL2文件")
            return None
        else:
            return [cls.process_molecule(mol, output_folder), ]

    # 根据文件格式选择处理函数
    @classmethod
    def generate_molecule_images_by_path(cls, file: str, output_folder):
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        file_extension = os.path.splitext(file)[-1].lower()

        if file_extension == '.sdf':
            return cls.process_sdf(file, output_folder)
        elif file_extension == '.pdb':
            return cls.process_pdb(file, output_folder)
        elif file_extension == '.mol':
            return cls.process_mol(file, output_folder)
        elif file_extension == '.mol2':
            return cls.process_mol2(file, output_folder)
        else:
            print(f"不支持的文件格式: {file_extension}")

    @classmethod
    def generate_molecule_images_by_bytes(cls, file: bytes, file_name: str, output_folder=None) -> list[str | None]:
        """
        根据文件名和文件内容生成分子图片
        :param file: 文件内容
        :param file_name: 文件名
        :param output_folder: 保存图片的路径
        :return: 图片的base64字符串
        """
        lower_file_name = file_name.lower()

        if lower_file_name.endswith('.sdf'):
            return cls.process_sdf(file, output_folder)
        elif lower_file_name.endswith('.pdb'):
            return cls.process_pdb(file, output_folder)
        elif lower_file_name.endswith('.mol'):
            return cls.process_mol(file, output_folder)
        elif lower_file_name.endswith('.mol2'):
            return cls.process_mol2(file, output_folder)
        else:
            logging.error(f"该文件格式不支持渲染图片: {file_name}")
            raise ValueError(f"该文件格式不支持渲染图片: {file_name}")


if __name__ == '__main__':
    # input_file = r'/Users/bigboss/Desktop/17455-13-9.mol'
    save_folder = 'molecule_images'

    import time

    # start_time = time.time()
    # print(ToolRendering2DStructureService.generate_molecule_images_by_path(input_file, save_folder))
    # end_time = time.time()
    # print(f"程序运行时间：{end_time - start_time} s")

    with open("/Users/bigboss/Desktop/ligs.sdf", "rb") as f:
        file_bytes = f.read()

    print(ToolRendering2DStructureService.generate_molecule_images_by_bytes(file_bytes, "ligs.sdf", save_folder))
