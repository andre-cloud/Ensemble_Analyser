# Conformer Ensemble Pruning Software

This software is a tool designed for molecular modeling to simplify and optimize their pruning and optimizing workflow of an ensemble of conformers. 

With this software, users can quickly and easily identify the most important conformers and reduce the computational burden of simulations.

---

## Installation

To use the software, follow these steps:

- Clone the repository onto your local machine using Git.
```bash
git clone https://github.com/username/ensemble-calculation.git
```

- Install the required dependencies listed in the requirements.txt file. This can be done by running the following command in your terminal:
```bash
pip install -r requirements.txt
```

- Install [ORCA](https://orcaforum.kofo.mpg.de/app.php/portal) from the ORCA Forum

- Export the ORCA command in order to start the ASE calculation
```bash
export ASE_ORCA_COMMAND="/complete/path/to/orca/folder/orca PREFIX.inp > PREFIX.out"
```

## Usage

Usage

- Prepare the ensemble:
    - Create a file containing the conformer structures in XYZ format.
    - Adjust the charge and multiplicity values as needed.

- Define the protocol:
    - Create a JSON file specifying the protocol steps, calculation levels, and thresholds.
    - Specify the functional, frequency calculation, and additional inputs for graph generation if required.

- Run the calculation:
```bash
python ensemble_calculation.py --ensemble <path/to/ensemble.xyz> --protocol <path/to/protocol.json> --cpu <#cpus>
```
   - Optional arguments:
        - ```--output <output_file>```: Specify the output file name (default: output.out).
        - ```--cpu <cpu_count>```: Set the number of CPUs to allocate (default: maximum available CPUs).
        - ```--temperature <temperature>```: Set the temperature in Kelvin (default: 298.15 K).
        - ```--final_lambda <final_lambda>```: Set the final lambda value for graph generation (default: 800.0).
        - ```--definition <definition>```: Set the definition value for graph generation (default: 4).
        - ```--fwhm <fwhm>```: Set the full width at half maximum (FWHM) for graph generation (default: None).
        - ```--shift <shift>```: Set the shift value for graph generation (default: None).
        - ```--invert```: Invert the energy axis for graph generation (default: False).
        - ```--restart```: Restart the calculation from a previous checkpoint.

- View the results:
    - The calculation output and log file will be saved in the current directory.
    - Summary tables and graphs will be generated for each protocol step.
    -  Final results and a summary will be displayed in the log.








<!-- To use the software, first, create two JSON files specifying the protocol to use and conformers' pruning parameters. 
The file format is described in detail
```bash
python ensemble_analyser.py -h-p # to obtain an example of the protocol.json file
python ensemble_analyser.py -h-t # to obtain an example of the treshold.json file
```

Next, run the `ensemble_analyser.py` script with the path to the JSON file as a command-line argument. The script will read the input file, perform the conformer pruning, and output a file with the pruned conformers.

```bash
python ensemble_analyser.py ensemble.xyz
````

## Parameters

The software uses the following parameters to determine the most relevant conformers:

- thrG: Refers to the energy (G or E) threshold to consider two conformers equivalent together with thrB.
- thrB: Refers to the rotary constant threshold to consider two conformers equivalent together with thrG.
- thrGMAX: Refers to the maximum energy window considered. Conformers lying above it will be sorted out immediately. -->



---

## Contributing

Contributions to this software are always welcome! If you have any ideas or suggestions, please feel free to submit a pull request or open an issue.

## License

This software is licensed under the MIT-3.0 License. See the LICENSE file for details.

## Contact

For any questions or support, please contact [by email](mailto:andrea.pellegrini15@unibo.it).
