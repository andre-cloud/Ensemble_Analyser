#!/data/bin/python_env/bin/python3

import datetime
import json, os
import sys
from ase.calculators.orca import ORCA

from ensemble_analyser.launch import launch
from ensemble_analyser.logger import save_snapshot
from ensemble_analyser.pruning import calculate_rel_energies, check_ensemble
from tabulate import tabulate

try:
    from .logger import ordinal
except ModuleNotFoundError:
    pass


DEBUG = os.getenv('DEBUG')


def load_protocol(file:str):
    default = 'ensemble_analyser/parameters_file/default_protocol.json'
    return json.load(open(default if not file else file ))

def load_threshold(file:str):
    default = 'ensemble_analyser/parameters_file/default_threshold.json'
    return json.load(open(default if not file else file ))




LEVEL_DEFINITION = {
    0 : 'SP'.lower()          , # mere energy calculation
    1 : 'OPT'.lower()         , # optimisation step
    2 : 'FREQ'.lower()        , # single point and frequency analysis
    3 : 'OPT+FREQ'.lower()      # optimisation and frequency analysis
}


class Solvent: 
    def __init__(self, solv:dict):
        self.solvent = solv['solvent']
        self.smd = solv['smd']

    def __str__(self):
        return f' - SMD({self.solvent})' if self.smd else  f' - CPCM({self.solvent})'
    def __repr__(self):
        return f' - SMD({self.solvent})' if self.smd else  f' - CPCM({self.solvent})'

    def orca_input_smd(self):
        if self.smd:
            return f'%cpcm smd true smdsolvent "{self.solvent}" end'
        return ''


class Protocol: 

    def __init__(self, number:int, functional:str, basis:str, solvent:Solvent, opt:bool, freq:bool, add_input:str , thrs_json, calculator='orca', thrG: float = None, thrB: float = None, thrGMAX: float = None):

        self.number = number
        self.functional = functional.upper()
        self.basis = basis.upper()
        self.solvent = solvent
        self.opt = opt
        self.freq = freq
        self.add_input = add_input
        self.thrG = thrG
        self.thrB = thrB
        self.thrGMAX = thrGMAX
        self.get_thrs(thrs_json)
        self.calculator = calculator

    @property
    def calculation_level(self):
        return LEVEL_DEFINITION[self.number_level].upper()

    @property
    def level(self):
        return f'{self.functional}/{self.basis}' + (str(self.solvent) if self.solvent else '')
    
    @property
    def thr(self):
        return f'\tthrG    : {self.thrG} kcal/mol\n\tthrB    : {self.thrB} cm-1\n\tthrGMAX : {self.thrGMAX} kcal/mol\n'
    
    @property
    def number_level(self):
        c = 0
        if self.opt: c += 1
        if self.freq: c += 2
        return c

    def get_calculator(self, cpu, charge:int, mult:int):
        calc = {
            'orca' : self.get_orca_calculator(cpu, charge, mult)
        }

        return calc[self.calculator]
        

    def get_thrs(self, thr_json):
        c = LEVEL_DEFINITION[self.number_level]
        if not self.thrG: self.thrG = thr_json[c]['thrG']
        if not self.thrB: self.thrB = thr_json[c]['thrB']
        if not self.thrGMAX: self.thrGMAX = thr_json[c]['thrGMAX']


    def __str__(self):
        return f'{self.functional}/{self.basis}' + (str(self.solvent) if self.solvent.solvent else '')
    def __repr__(self):
        return f'{self.functional}/{self.basis}' + (str(self.solvent) if self.solvent.solvent else '')
    

    def get_orca_calculator(self, cpu:int, charge:int, mult:int):
        # possibilities for solvent definitions
        if self.solvent.solvent:
            if 'xtb' in self.functional.lower():
                solv = f"ALPB({self.solvent.solvent})"
            else:
                solv = f" cpcm({self.solvent.solvent})"
        else:
            solv = ''

        # ! B3LYP def2-SVP FREQ CPCM(solvent)
        simple_input = f'{self.functional} {self.basis} {"freq" if self.freq else ""} {"opt" if self.opt else ""} {solv} nopop'


        # %cpcm
        #     smd True
        #     smdSolvent solvent
        # end
        smd = ''
        if self.solvent and 'xtb' not in self.functional.lower(): smd = self.solvent.orca_input_smd()

        label = 'ORCA'
        calculator = ORCA(
            label = label,
            orcasimpleinput = simple_input,
            orcablocks=f'%pal nprocs {cpu} end ' + smd + self.add_input + ('%maxcore 4000' if 'maxcore' not in self.add_input else ''),
            charge = charge, 
            mult = mult, 
            task='energy'
        )

        return calculator, label
    
    @staticmethod
    def load_raw(json):
        return Protocol(
            number = json['number'],
            functional = json['functional'],
            basis = json['basis'],
            solvent = Solvent(json['solvent']),
            opt = json['opt'],
            freq = json['freq'], 
            add_input = json['add_input'],
            calculator=json['calculator'],
            thrs_json = None, 
            thrB=json['thrB'], 
            thrG=json['thrG'], 
            thrGMAX=json['thrGMAX']
        )



def create_protocol(p, thrs, log):
    protocol = []

    log.info('Loading Protocol\n')
    for idx, d in p.items():
        func    = d.get('func', None)
        basis   = d.get('basis', 'def2-svp')
        opt     = d.get('opt', False)
        freq    = d.get('freq', False)

        add_input = d.get('add_input', '')

        solv    = d.get('solv', None)
        if solv or solv.get('solvent', None): solv = Solvent(solv)

        thrG    = d.get('thrG', None)
        thrB    = d.get('thrB', None)
        thrGMAX = d.get('thrGMAX', None)

        if not func:
            log.critical(f"{'='*20}\nCRITICAL ERROR\n{'='*20}\nFUNC key must be passed in order to calculate energy. DFT functional or HF for Hartree-Fock calculation or semi-empirical methods (XTB1/XTB2/PM3/AM1 or similar supported by the calculator) (Probelm at {ordinal(int(idx))} protocol definition)\n{'='*20}\n")
            raise IOError('There is an error in the input file with the definition of the functional. See the output file.')

        protocol.append(Protocol(number =idx, functional=func, basis=basis, solvent= solv, opt= opt, freq= freq, add_input = add_input, thrs_json= thrs, thrG = thrG, thrB = thrB, thrGMAX = thrGMAX))

    log.info('\n'.join((f"{i.number}: {str(i)} - {i.calculation_level}\n {i.thr}" for i in protocol)) + '\n')
    return protocol


def run_protocol(conformers, p, temperature, cpu, log):
    log.info(f'STARTING PROTOCOL {p.number}')
    log.info(f'\nActive conformers for this phase: {len([i for i in conformers if i.active])}\n')
    for i in conformers:
        if not i.active: continue
        if i.energies.get(str(p.number)): continue
        launch(i, p, cpu, log, conformers)

    conformers_tmp = sorted(conformers)

    calculate_rel_energies(conformers, temperature)

    log.info('')
    log.info('Ended Calculations')
    log.info('')
    log.info('Summary')
    log.info('')
    log.info(tabulate([i.create_log() for i in conformers_tmp if i.active], headers=['conformers', 'E[Eh]' ,'G[Eh]', 'B[cm-1]', 'E. Rel [kcal/mol]', 'Pop [%]', 'Elap. time [sec]'], floatfmt=".6f"))
    log.info('')
    log.info('Total elapsed time: ' + str(datetime.timedelta(seconds = sum([i._last_energy['time'] for i in conformers_tmp if i.active]))))

    log.info('Start Pruning')
    conformers = check_ensemble(conformers, p, log)
    save_snapshot(f'ensemble_after_{p.number}.xyz', conformers, log)

    log.info(f'{"="*15}\nEND PROTOCOL {p.number}\n{"="*15}\n\n')
    #{'charge': 0, 'mult': 1, 'task': 'gradient', 'orcasimpleinput': 'B3LYP def2-svp ', 'orcablocks': ''}



