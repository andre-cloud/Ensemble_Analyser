import ase, os, re
from ase import Atoms
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.vibrations import Vibrations, Infrared

from src.IOsystem import tail
from src.regex_parsing import regex_parsing
from src.parser_parameter import get_opt_geometry

MAX_TRY = 5

def check_output(label, calc):
    file = f"{label}.{regex_parsing[calc]['ext']}"
    if os.path.exists(file):
        return regex_parsing[calc]['finish'] in tail(file, 5)
    
    return False
    

def optimize(conf, protocol, cpu : int, log, try_ = 0):
    calc, label = protocol.get_calculator(
            cpu=cpu, charge=conf.charge, mult=conf.mult, mode="opt"
        )
    
    atoms = conf.get_ase_atoms(calc)

    if protocol.constrains: 
        c = FixAtoms(protocol.constrains)
        atoms.set_constraint(c)
    
    opt = BFGS(atoms, maxstep=protocol.maxstep, logfile=f"opt_p{protocol.number}_{conf.number}.log", trajectory=f"p{protocol.number}_{label}.trj")
    opt.run(protocol.fmax)

    if not check_output(label, protocol.calculator):
        if try_ < MAX_TRY:
            optimize(conf, protocol, cpu, log, try_+1)
        else:
            raise 

    set_last_geometry(conf, atoms.get_positions())

    os.rename(f"opt_p{protocol.number}_{conf.number}.log", f"{conf.folder}/opt_p{protocol.number}_{conf.number}.log")
    os.rename(f"p{protocol.number}_{label}.trj", f"{conf.folder}/p{protocol.number}_{label}.trj")


    return atoms, label

def calc_freq(conf, protocol, cpu : int, log, try_ = 0):
    calc, label = protocol.get_calculator(
            cpu=cpu, charge=conf.charge, mult=conf.mult, mode="freq"
        )
    
    atoms = conf.get_ase_atoms(calc)

    vib = Infrared(atoms)
    try:
        vib.run()
    except ase.calculators.calculator.PropertyNotImplementedError: 
        pass

    return atoms, label

def single_point(conf, protocol, cpu : int, log, try_=0):


    calc, label = protocol.get_calculator(
            cpu=cpu, charge=conf.charge, mult=conf.mult, mode="energy"
        )
    
    atoms = conf.get_ase_atoms(calc)

    if "opt" in protocol.functional.lower():
        with open(f"{label}.{regex_parsing[calc]['ext']}") as f:
            fl = f.readlines()

        geom = get_opt_geometry(fl, protocol.calculator, log)
        set_last_geometry(conf, geom)

    atoms.get_potential_energy()
    
    return atoms, label




def set_last_geometry(conf, geometry):

    conf.last_geometry = geometry[:]
    return None