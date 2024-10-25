from io import BufferedReader


class ToolLigandInfoService:

    @classmethod
    def read_pdb(cls, file_path):
        with open(file_path, 'r') as file:
            return file.readlines()

    @classmethod
    def is_ligand(cls, line):
        return line.startswith('HETATM') and line[17:20].strip() not in ['HOH', 'WAT']

    @classmethod
    def get_unique_ligands(cls, pdb_lines):
        ligands = set()
        for line in pdb_lines:
            if cls.is_ligand(line):
                chain = line[21]
                residue_number = int(line[22:26].strip())
                ligands.add((chain, residue_number))
        return ligands

    @classmethod
    def extract_unique_ligands(cls, pdb_file):
        if isinstance(pdb_file, BufferedReader):
            pdb_lines = pdb_file.read().decode('utf-8').split('\n')
        elif isinstance(pdb_file, str):
            pdb_lines = cls.read_pdb(pdb_file)
        elif isinstance(pdb_file, bytes):
            pdb_lines = pdb_file.decode().split('\n')
        else:
            raise ValueError("Invalid input type. Choose either a file path or a string of pdb lines.")
        unique_ligands = cls.get_unique_ligands(pdb_lines)
        return unique_ligands


if __name__ == '__main__':
    # Example usage:
    filename = '/Users/bigboss/pycharm_project/alphama-dify/api/storage/upload_files/1e74140f-ac64-4c6c-9cd6-489b2d326ee3/4968b675-3dd1-4be5-9451-b486af9d500c.pdb'
    unique_ligands_info = ToolLigandInfoService.extract_unique_ligands(filename)
    print(unique_ligands_info)

    dd = ('A', 1801)

    if dd in unique_ligands_info:
        print('Yes')
    else:
        print('No')
