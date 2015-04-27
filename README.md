# Pulseq

This project is an open source framework for the development, 
representation and execution of magnetic resonance (MR) sequences. A central contribution 
of this project is an **open file format** to compactly describe MR sequences 
suitable for execution on an MR scanner. 
MATLAB and C++ source code is provided for reading and writing sequence files. The main homepage for Pulseq is
http://github.com/. This is the GitHub repository of the source code.



This project is open source under the MIT License. See [LICENSE](LICENSE) for details.

The directories are organized as follows:

* `doc/` - Contains the file specification and HTML source code documentation
* `examples/` - Contains example sequence files (`*.seq`)
* `src/` - C++ class for reading sequence files
* `matlab/` - MATLAB code for reading, writing, modifying and visualizing sequence files

## System requirements

System requirements vary depending on which features one intends to use or modify. In short,

- doxygen is required to generate HTML source code documentation
- latex is required to build the file specification PDF
- python is required to run the build tests (i.e. `make check`)

These are optional and not essential to start using Pulseq.


