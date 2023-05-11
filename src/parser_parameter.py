import re
import os
import numpy as np
import cclib

try:
    from src.regex_parsing import regex_parsing
    from src.rrho import free_gibbs_energy
except ImportError as e:
    print(e)
    from regex_parsing import regex_parsing
    from rrho import free_gibbs_energy

EH_TO_KCAL = 627.5096080305927


def get_param(x, calculator, param):
    """
    Parsing for Rotational Constant
    """
    if re.search(regex_parsing[calculator][param], x):
        return x


def get_freq(fl: str, calc: str) -> np.ndarray:
    """
    Parsing for frequencies

    fl | str : piece of log file
    calc | str : calculator type
    return | np.array : frequency [cm-1]
    """

    fl = (
        "\n".join(
            "".join(fl)
            .split(regex_parsing[calc]["s_freq"])[-1]
            .strip()
            .splitlines()[4:]
        )
        .split(regex_parsing[calc]["e_freq"])[0]
        .strip()
        .splitlines()
    )
    freq = np.array(
        [
            float(i.split()[regex_parsing[calc]["idx_freq"]])
            for i in fl
            if float(i.split()[regex_parsing[calc]["idx_freq"]]) != 0.0
        ]
    )
    return freq


def get_conf_parameters(conf, number: int, p, time, temp: float, log) -> bool:
    """
    Obtain the parameters for a conformer: E, G, B, m

    conf | Conformer : conformer
    number | int : protocol number
    p | Protocol : protocol executed
    time | datetime : elapsed time requested for the calculation
    temp | float : temperature [K]
    log : logger instance

    return | bool : calculation ended correctly and not crashed due to server error
    """

    data = cclib.io.ccread(f'ORCA.out')
    conf.last_geometry = data.atomcoords[-1]
    
    with open(os.path.join(conf.folder, f"protocol_{number}.out")) as f:
        fl = f.readlines()

    try:
        e = float(
            list(filter(lambda x: get_param(x, p.calculator, "E"), fl))[-1]
            .strip()
            .split()[-1]
        )
    except Exception as e:
        log.error(e)
        return False

    freq = np.array([])
    if p.freq:
        freq = get_freq(fl, p.calculator) * p.freq_fact
        log.info(f"{conf.number} has {freq[freq<0].size} imaginary frequency(s)")
        if freq.size == 0:
            log.error(("\n".join(fl[-6:])).strip())
            log.critical(
                f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nNo frequency present in the calculation output.\n{'='*20}\nExiting\n{'='*20}\n"
            )
            raise IOError("No frequency in the output file")

    B = np.array(
        list(filter(lambda x: get_param(x, p.calculator, "B"), fl))[-1]
        .strip()
        .split(":")[-1]
        .split(),
        dtype=float,
    )
    b = np.linalg.norm(B)

    M = np.linalg.norm(
        np.array(
            list(filter(lambda x: get_param(x, p.calculator, "m"), fl))[-1]
            .strip()
            .split(":")[-1]
            .split(),
            dtype=float,
        )
    )

    g = ""
    if freq.size > 0:
        g = free_gibbs_energy(
            SCF=e, T=temp, freq=freq, mw=conf.weight_mass, B=B, m=conf.mult
        )

    conf.energies[str(number)] = {
        "E": e * EH_TO_KCAL if e else e,  # Electronic Energy [kcal/mol]
        "G": g * EH_TO_KCAL if g else None,  # Free Gibbs Energy [kcal/mol]
        "B": b if b else None,  # Rotatory Constant [cm-1]
        "m": M if M else None,  # dipole momenti [Debye]
        "time": time,  # elapsed time [sec]
    }

    return True


if __name__ == "__main__":
    # from logger import create_log

    # class Conf:
    #     def __init__(self, number, mult, folder):
    #         self.number = number
    #         self.mult = mult
    #         self.folder = folder
    #         self.energies = {}

    # c = Conf('1', 1, 'conf_1')
    # log = create_log('test.out')

    # get_conf_parameters(c, 0, 1, 298.15, log)
    # print(c.energies)

    with open("files/opt.out") as f:
        fl = f.readlines()
    get_freq(fl)
