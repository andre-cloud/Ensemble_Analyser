try:
    from src.IOsystem import mkdir
except ImportError as e:  # pragma: no cover
    print(e)
    from IOsystem import mkdir

import numpy as np
import random
from ase.atoms import Atoms


class Conformer:
    """
    Storing all the information on each conformer for all the parts of the protocol
    """

    def __init__(
        self,
        number: int,
        geom: np.ndarray,
        atoms: np.ndarray,
        charge: int = 0,
        mult: int = 1,
        raw=False,
    ) -> None:
        self.number = number
        self._initial_geometry = geom.copy()
        self.charge = charge
        self.mult = mult

        self.last_geometry = geom.copy()
        self.atoms = atoms.copy()
        self.energies = {}
        self.active = True
        self.color = "#%06x" % random.randint(0, 0xFFFFFF)

        self.cluster = None
        # IO
        self.folder = f"conf_{self.number}"
        if not raw:
            mkdir(self.folder)

    def get_ase_atoms(self, calc=None):
        """Return the atoms needed for the calculation

        :param calc: the type of calculator to use
        :type calc: ase.calculator
        :return: the ase instance ready to start the calculation
        :rtype: ase.Atoms
        """
        return Atoms(
            symbols="".join(list(self.atoms)),
            positions=self.last_geometry,
            calculator=calc,
        )

    @property
    def weight_mass(self):
        return np.sum(
            Atoms(
                symbols="".join(list(self.atoms)),
                positions=self.last_geometry,
            ).get_masses()
        )

    @property
    def rotatory(self):
        return self.energies[list(self.energies.keys())[-1]]["B"]

    @property
    def moment(self):
        return self.energies[list(self.energies.keys())[-1]]["m"]

    @property
    def get_energy(self):
        en = self._last_energy
        if en["G"]:
            return en["G"]
        return en["E"]

    @property
    def _last_energy(self):
        if self.energies:
            return self.energies[list(self.energies.keys())[-1]]
        return {"E": 0, "G": None}

    def write_xyz(self):
        """Write the XYZ string to be stored in a file

        :return: the string in the XYZ formatting
        :rtype: str
        """
        if not self.active:
            return ""
        txt = f'{len(self.atoms)}\nCONFORMER {self.number} {"G : {:.6f} kcal/mol".format(self._last_energy["G"]) if self._last_energy["G"] else "E : {:.6f} kcal/mol".format(self._last_energy["E"])}\n'
        for a, pos in zip(self.atoms, self.last_geometry):
            x, y, z = pos
            txt += f" {a}\t{x:14f}\t{y:14f}\t{z:14f}\n"
        return txt.strip()

    def create_log(self):
        """Generate all the information needed for the tabulation

        :return: a long tuple with all the information. (Number, E, G, B, Erel, Pop, Elapsed Time)
        :rtype: tuple
        """
        en = self._last_energy
        number, e, g, b, erel, time, pop = (
            self.number,
            en.get("E", float(0)),
            en.get("G", float(0)),
            en.get("B", float(0)),
            en.get("Erel", float(0)),
            en.get("time"),
            en.get("Pop", float(0)),
        )
        if g:
            g /= 627.51
        return number, e / 627.51, g, b, erel, pop, time

    @staticmethod
    def load_raw(json):
        """Load raw configuration. Restart purposes.

        :param json: All the information of the conformer
        :type json: dict
        :return: the conformer instance
        :rtype: Conformer
        """
        a = Conformer(
            number=json["number"],
            geom=json["last_geometry"],
            atoms=json["atoms"],
            charge=json["charge"],
            mult=json["mult"],
            raw=True,
        )
        a.energies = json["energies"]
        a.active = json["active"]
        # a.cluster = json["cluster"]
        return a

    def __str__(self) -> str:
        return str(self.number)

    def __repr__(self) -> str:
        return str(self.number)

    # Functions needed for sorting the conformers' ensemble

    def __lt__(self, other):
        if not self.active:
            return 0 < other.get_energy
        return self.get_energy < other.get_energy

    def __gt__(self, other):
        if not self.active:
            return 0 > other.get_energy
        return self.get_energy > other.get_energy

    def __eq__(self, other):
        if not self.active:
            return 0 == other.get_energy
        return self.get_energy == other.get_energy
