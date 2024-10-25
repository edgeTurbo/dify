import io
from io import BufferedReader


class ToolDeleteSpecialLigandService:

    @classmethod
    def read_pdb(cls, file_path):
        with open(file_path, 'r') as file:
            return file.readlines()

    @classmethod
    def write_pdb(cls, file_path, lines):
        with open(file_path, 'w') as file:
            file.writelines(lines)

    @classmethod
    def is_ligand(cls, line, chain, residue_number):
        return (line.startswith('HETATM') and
                line[21] == chain and
                int(line[22:26]) == residue_number)

    @classmethod
    def remove_ligand(cls, pdb_file, chain, residue_number, output_file=None):
        file_name = None
        if isinstance(pdb_file, BufferedReader):
            file_name = pdb_file.name
            # 重置指针到开头
            pdb_file.seek(0)
            pdb_lines = pdb_file.read().decode('utf-8').split('\n')
        elif isinstance(pdb_file, str):
            pdb_lines = cls.read_pdb(pdb_file)
        elif isinstance(pdb_file, bytes):
            pdb_lines = pdb_file.decode().split('\n')
        else:
            raise ValueError("Invalid input type. Choose either a file path or a string of pdb lines.")

        new_lines = [line for line in pdb_lines if not cls.is_ligand(line, chain, residue_number)]
        if output_file is not None:
            cls.write_pdb(file_name, new_lines)
            return output_file
        else:
            # Return a buffered reader object
            binary_data = '\n'.join(new_lines).encode('utf-8')
            bytes_io = io.BytesIO(binary_data)
            bytes_io.name = file_name
            buffered_reader = io.BufferedReader(bytes_io)
            return buffered_reader


if __name__ == '__main__':
    # Example usage:
    input_filename = '/Users/bigboss/Desktop/rec.pdb'
    output_filename = '5v3x_modified.pdb'
    chain_to_remove = 'A'  # Specify the chain
    residue_number_to_remove = 1801  # Specify the residue number

    ToolDeleteSpecialLigandService.remove_ligand(input_filename, chain_to_remove, residue_number_to_remove, output_filename)
    print(f'Removed ligand from {input_filename} and saved to {output_filename}.')