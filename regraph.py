#!/data/bin/python_env/bin/python3

from src.launch import restart, create_protocol
from src.logger import create_log
from src.grapher import Graph
import json, logging
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('idx', nargs='+', help="Protocol's number to (re-)generate the graphs")


args = parser.parse_args()

conformers, protocol, start_from = restart()

settings = json.load(open("settings.json"))
output = 'regraph.log'
temperature = settings.get("temperature", 298.15)
final_lambda = settings.get("final_lambda", 800)
definition = settings.get("definition", 4)
fwhm = settings.get("fwhm", None)
shift = settings.get("shift", None)
invert = settings.get("invert", False)

# initiate the log
log = create_log(output)
# deactivate the log of matplotlib
logging.getLogger("matplotlib").disabled = False

protocol = create_protocol('protocol_dump.json', log)



for p in args.idx: 
    Graph(
        confs=conformers,
        protocol=p,
        log=log,
        T=temperature,
        final_lambda=final_lambda,
        definition=definition,
        FWHM=fwhm,
        shift=shift,
        invert=invert,
        regraph = True
    )
