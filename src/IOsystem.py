import numpy as np
import shutil
import json
import os


def _parse_xyz_str(fl: str) -> tuple:
    """
    Parse an xyz geom descriptor

    fl | str : string of the file splitted in a single geometry

    return | tuple : list of atoms, XYZ position of the atoms
    """
    fl = fl[2:]
    atoms, geom = [], []
    for line in fl:
        a, *g = line.split()
        atoms.append(a)
        geom.append(g)
    return np.array(atoms), np.array(geom, dtype=float)


def mkdir(directory: str) -> None:
    """
    Create a directory, raising an error if the directory already exists

    directory | str : directory name

    return None
    """
    if os.path.exists(directory):
        raise IOError(f"Directory {directory} already exists. Going to exit!")
        shutil.rmtree(directory)
    os.mkdir(directory)
    return None


class SerialiseEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj.__dict__
