import json
import os
import sys
from ase.calculators.orca import ORCA


DEBUG = os.getenv("DEBUG")


def load_protocol(file: str):  # pragma: no cover
    default = "ensemble_analyser/parameters_file/default_protocol.json"
    return json.load(open(default if not file else file))


LEVEL_DEFINITION = {
    0: "SP".lower(),  # mere energy calculation
    1: "OPT".lower(),  # optimisation step
    2: "FREQ".lower(),  # single point and frequency analysis
    3: "OPT+FREQ".lower(),  # optimisation and frequency analysis
}


class Solvent:
    """
    Solvent class
    """

    def __init__(self, solv: dict):
        self.solvent = solv["solvent"]
        self.smd = solv["smd"]

    def __str__(self):  # pragma: no cover
        if self.smd:
            return f"SMD({self.solvent})"
        elif self.solvent:
            return f"CPCM({self.solvent})"
        else: 
            return "CPCM"

    def __repr__(self):  # pragma: no cover
        if self.smd:
            return f"SMD({self.solvent})"
        elif self.solvent:
            return f"CPCM({self.solvent})"
        else: 
            return "CPCM"
        
    def orca_input_smd(self):
        if self.smd:
            return f'%cpcm smd true smdsolvent "{self.solvent}" end'
        return ""


class Protocol:
    def __init__(
        self,
        number: int,
        functional: str,
        basis: str = "def2-svp",
        solvent: dict = {},
        opt: bool = False,
        freq: bool = False,
        add_input: str = "",
        freq_fact: float = 1,
        graph: bool = False,
        calculator: str = "orca",
        thrG: float = None,
        thrB: float = None,
        thrGMAX: float = None,
        constrains: list = [],
        maxstep : float = 0.2,
        fmax : float = 0.01, 
    ):
        self.number = number
        self.functional = functional.upper()
        self.basis = basis.upper()
        self.solvent = Solvent(solvent) if solvent else None
        self.opt = opt
        self.freq = freq
        self.add_input = add_input
        self.thrG = thrG
        self.thrB = thrB
        self.thrGMAX = thrGMAX
        self.get_thrs(self.load_threshold())
        self.calculator = calculator
        self.constrains = constrains
        self.maxstep = maxstep
        self.fmax = fmax

        self.freq_fact = freq_fact
        self.graph = graph

    @property
    def calculation_level(self):
        return LEVEL_DEFINITION[self.number_level].upper()

    @property
    def level(self):
        return f"{self.functional}/{self.basis}" + (
            str(self.solvent) if self.solvent else ""
        )

    @property
    def thr(self):
        return f"\tthrG    : {self.thrG} kcal/mol\n\tthrB    : {self.thrB} cm-1\n\tthrGMAX : {self.thrGMAX} kcal/mol\n"

    @property
    def number_level(self):
        c = 0
        if self.opt:
            c += 1
        if self.freq:
            c += 2
        return c

    def load_threshold(self) -> dict:
        """
        Load default thresholds

        return | dict : thresholds
        """

        default = os.path.join(
            os.path.abspath(os.path.dirname(__file__)),
            "parameters_file",
            "default_threshold.json",
        )
        return json.load(open(default))

    def get_calculator(self, cpu, charge: int, mult: int, mode: str):
        """
        Get the calculator from the user selector

        cpu | int : allocated CPU
        charge | int : charge of the molecule
        mult | int : multiplicity of the molecule
        mode | str {opt, freq, energy}: type of calculation required
        """

        calc = {"orca": {
            "opt": self.orca_opt,
            "freq": self.orca_freq,
            "energy": self.calc_orca_std,
            },
        }

        return calc[self.calculator][mode](cpu, charge, mult)

    def get_thrs(self, thr_json):
        """
        Get default thrs if not defined by user

        thr_json | dict : JSON default thresholds
        """
        c = LEVEL_DEFINITION[self.number_level]
        if not self.thrG:
            self.thrG = thr_json[c]["thrG"]
        if not self.thrB:
            self.thrB = thr_json[c]["thrB"]
        if not self.thrGMAX:
            self.thrGMAX = thr_json[c]["thrGMAX"]

    def __str__(self):  # pragma: no cover
        if self.solvent:
            return f"{self.functional}/{self.basis} - {self.solvent}"
        return f"{self.functional}/{self.basis}"

    def __repr__(self):  # pragma: no cover
        if self.solvent:
            return f"{self.functional}/{self.basis} - {self.solvent}"
        return f"{self.functional}/{self.basis}"






    # ORCA CALCULATOR

    def orca_common_str(self, cpu):
        if self.solvent:
            if "xtb" in self.functional.lower():
                solv = f"ALPB({self.solvent.solvent})"
            else:
                solv = f" cpcm({self.solvent.solvent})"
        else:
            solv = ""

        si = f'{self.functional} {self.basis} {solv} nopop'

        smd = ""
        if self.solvent and "xtb" not in self.functional.lower():
            smd = self.solvent.orca_input_smd()

        ob = (f"%pal nprocs {cpu} end "
            + smd
            + self.add_input
            + (" %maxcore 4000" if "maxcore" not in self.add_input else ""))

        return si, ob

    def calc_orca_std(self, cpu: int, charge: int, mult: int):

        simple_input, ob = self.orca_common_str(cpu)
        label = "ORCA"
        calculator = ORCA(
            label=label,
            orcasimpleinput=simple_input,
            orcablocks=ob,
            charge=charge,
            mult=mult,
            task="energy"
        )

        return calculator, label
    
    def orca_opt(self, cpu: int, charge: int, mult: int):
        calculator, label = self.calc_orca_std(cpu, charge, mult)
        calculator.parameters["orcasimpleinput"] += " engrad"

        return calculator, label
    
    def orca_freq(self, cpu: int, charge: int, mult: int):
        calculator, label = self.calc_orca_std(cpu, charge, mult)
        calculator.parameters["orcasimpleinput"] += " freq"

        return calculator, label





    @staticmethod
    def load_raw(json):
        return Protocol(
            number=json["number"],
            functional=json["functional"],
            basis=json["basis"],
            solvent=Solvent(**json["solvent"]) if json["solvent"] else None,
            opt=json["opt"],
            freq=json["freq"],
            add_input=json["add_input"],
            calculator=json["calculator"],
            thrB=json["thrB"],
            thrG=json["thrG"],
            thrGMAX=json["thrGMAX"],
            freq_fact=json["freq_fact"],
            graph=json["graph"],
            constrains=json["constrains"], 
            maxstep = json["maxstep"],
            fmax=json["fmax"]
        )


