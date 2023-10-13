.. modules:

Code
====

.. contents:: Table of Contents

IOsystem Module
---------------

This module provides versatile functionalities. It includes a function to parse the `xyz geometry descriptor` data format, yielding lists of atoms and their XYZ positions.

The module also offers a directory creation function, raising an exception if the directory already exists. Additionally, it contains a custom JSON encoder class, `SerialiseEncoder`, which facilitates the serialization of objects like NumPy arrays. Lastly, there's a function for displaying the last lines of an output file, useful for quick inspection.

.. automodule:: src.IOsystem
   :members:
   :undoc-members:
   :show-inheritance:

Calculations Module
--------------------

This module contains functions for performing various computational chemistry calculations, including *geometry optimization*, *vibrational frequency* calculations, and *single-point energy* calculations. It interfaces with the Atomic Simulation Environment (ASE).

The module includes functions for checking calculation output, handling calculation errors, and setting the last calculated geometry.

.. automodule:: src.calculations
   :members:
   :undoc-members:
   :show-inheritance:

Clustering Module
------------------

This module leverages the power of *Machine Learning* and *Data Analysis* tools from the scikit-learn library to perform **Principal Component Analysis** (PCA) and **cluster analysis** on molecular conformations. It enables users to discover relationships between conformers based on their Cartesian coordinates.

The module provides functions to find the optimal number of clusters using silhouette scores, perform PCA, and create visual snapshots of the analysis. These tools are valuable for gaining insights into the structural diversity of molecular ensembles.

.. automodule:: src.clustering
   :members:
   :undoc-members:
   :show-inheritance:

Conformer Module
------------------

This module contains the `Conformer` class, which is responsible for storing and managing information about conformers used during the computational protocol. It provides functionalities to create and manage conformers, including their geometries, energies, and other properties.

The `Conformer` class is equipped with methods for data storage, file I/O, and serialization, making it a fundamental component for handling molecular structures in computational chemistry applications.

.. automodule:: src.conformer
   :members:
   :undoc-members:
   :show-inheritance:

Grapher Module
----------------

This `Graph` module provides functionality for generating electronic spectra graphs, including **UV** and **ECD** (Electronic Circular Dichroism) spectra.

It can automatically align calculated spectra with reference data, perform Gaussian convolution, and offer various customization options.

.. automodule:: src.grapher
   :members:
   :undoc-members:
   :show-inheritance:

ioFile Module
----------------

This module provides a set of functions for working with molecular conformer ensembles. It includes functionalities for converting input files into multi-geometry XYZ files, reading initial ensembles in various formats, and saving snapshots of conformer ensembles to XYZ files.

The module also utilizes Open Babel for file conversion.

.. automodule:: src.ioFile
   :members:
   :undoc-members:
   :show-inheritance:

Launch Module
----------------

The module includes various functions and procedures related to running calculations on the molecular conformers using a user-defined protocol.

.. automodule:: src.launch
   :members:
   :undoc-members:
   :show-inheritance:

Logger Module
---------------

The module includes the main logger function.

.. automodule:: src.logger
   :members:
   :undoc-members:
   :show-inheritance:

Parser Arguments Module
--------------------------

The module includes the ``argparse`` parsing method.

.. automodule:: src.parser_arguments
   :members:
   :undoc-members:
   :show-inheritance:

Parser Parameter Module
-------------------------

The module provides essential functionalities for extracting critical parameters and properties from output files. These functionalities enable streamlined data extraction and analysis.

.. automodule:: src.parser_parameter
   :members:
   :undoc-members:
   :show-inheritance:

Protocol Module
-----------------

This module provides a structured framework for configuring and running quantum chemistry calculations, streamlining the setup and execution of calculations in a computational chemistry workflow. For now, it is particularly focused for users working with the ORCA quantum chemistry software, but it will be adapted for other quantum chemistry software packages.

.. automodule:: src.protocol
   :members:
   :undoc-members:
   :show-inheritance:

Pruning Module
----------------

This module offers essential functions for managing molecular conformers and performing various analyses in computational chemistry workflows. Key functionalities include **energy-based** pruning, **structural comparison** using **RMSD**, and conformer deactivation.

.. automodule:: src.pruning
   :members:
   :undoc-members:
   :show-inheritance:

Rrho Module
-------------

This module provides various functions for calculating thermodynamic properties of molecular systems, and it's based upon the thermodynamic formalism of `Grimme`_ with the dumping of the frequency value.

.. _Grimme: https://doi.org/10.1039/D1SC00621E
.. automodule:: src.rrho
   :members:
   :undoc-members:
   :show-inheritance:
