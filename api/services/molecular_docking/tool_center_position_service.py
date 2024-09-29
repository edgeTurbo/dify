import numpy as np
from collections import defaultdict


class ToolCenterPositionService:

    @classmethod
    def read_pdb(cls, file_path):
        with open(file_path, 'r') as file:
            return file.readlines()

    @classmethod
    def is_valid_residue(cls, line):
        return line.startswith(('ATOM', 'HETATM')) and line[17:20].strip() not in ['HOH', 'WAT']

    @classmethod
    def is_ligand(cls, line):
        return line.startswith('HETATM') and line[17:20].strip() not in ['HOH', 'WAT']

    @classmethod
    def get_coordinates(cls, line):
        return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])

    @classmethod
    def calculate_center(cls, coordinates):
        return np.mean(coordinates, axis=0)

    @classmethod
    def get_residue_centers(cls, pdb_lines):
        residues = defaultdict(list)
        for line in pdb_lines:
            if cls.is_valid_residue(line):
                residue_name = line[17:20].strip()
                chain = line[21]
                residue_number = int(line[22:26])
                coordinates = cls.get_coordinates(line)
                residues[(residue_name, chain, residue_number)].append(coordinates)

        centers = {}
        for (residue_name, chain, residue_number), coords in residues.items():
            centers[(residue_name, chain, residue_number)] = cls.calculate_center(coords)

        return centers

    @classmethod
    def get_first_ligand_center(cls, pdb_lines, centers):
        for line in pdb_lines:
            if cls.is_ligand(line):
                residue_name = line[17:20].strip()
                chain = line[21]
                residue_number = int(line[22:26])
                key = (residue_name, chain, residue_number)
                if key in centers:
                    return centers[key]
        return None

    @classmethod
    def get_box_center(cls, pdb_file, method='default', **kwargs):
        if isinstance(pdb_file, str):
            pdb_lines = cls.read_pdb(pdb_file)
        elif isinstance(pdb_file, bytes):
            pdb_lines = pdb_file.decode().split('\n')
        else:
            raise ValueError("Invalid input type. Choose either a file path or a string of pdb lines.")
        centers = cls.get_residue_centers(pdb_lines)
        if method == 'default':
            info = cls.get_first_ligand_center(pdb_lines, centers)
            if info is not None:
                info.tolist()
                return {
                    'center_x': info[0],
                    'center_y': info[1],
                    'center_z': info[2]
                }
            return {
                    'center_x': 0,
                    'center_y': 0,
                    'center_z': 0
                }

        elif method == 'specific':
            specific_name = kwargs.get('specific_name')
            chain = kwargs.get('chain')
            residue_number = kwargs.get('residue_number')

            key = (specific_name, chain, residue_number)
            info = centers.get(key, None)
            if info is not None:
                info.tolist()
                return {
                    'center_x': info[0],
                    'center_y': info[1],
                    'center_z': info[2]
                }
            return {
                'center_x': 0,
                'center_y': 0,
                'center_z': 0
            }
        else:
            raise ValueError("Invalid method. Choose 'default' or 'specific'.")


if __name__ == '__main__':
    # Example usage:
    # filename = '5v3x_noligand.pdb'
    filename = '/Users/bigboss/Desktop/5v3x.pdb'
    default_center = ToolCenterPositionService.get_box_center(filename)
    print(type(default_center))
    print('default_center', default_center)
    specific_center = ToolCenterPositionService.get_box_center(filename, method='specific', specific_name='I28',
                                                               chain='A', residue_number=1801)
    print('specific_center_ligand', specific_center)
    # specific_center = get_box_center(filename, method='specific', specific_name='ASP', chain='B', residue_number=1661)
    # print('specific_center_residue',specific_center)
