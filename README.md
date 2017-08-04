<div align="center">
<a href="http://pulseq.github.io/" target="_blank">
<img src="http://pulseq.github.io/logo_hires.png" width="250" alt="Pulseq logo"></img>
</a>
</div>

# Introduction

*Pulseq* is an open source framework for the development, 
representation and execution of magnetic resonance (MR) sequences. A central contribution 
of this project is an **open file format** to compactly describe MR sequences 
suitable for execution on an MRI scanner or NMR spectrometer. 
MATLAB and C++ source code is provided for reading and writing sequence files. The main homepage for Pulseq is
http://pulseq.github.io/. This readme is part of the GitHub repository of the source code.



This project is open source under the MIT License. See [LICENSE](LICENSE) for details.

The directories are organized as follows:

* `doc/` - Contains the file specification and HTML source code documentation
* `examples/` - Contains example sequence files (`*.seq`)
* `src/` - C++ class for reading sequence files
* `matlab/` - MATLAB code for reading, writing, modifying and visualizing sequence files

## Performance
Improved performance of the MATLAB toolbox can be obtained by compiling key functions into mex functions. This is **not** necessary to use the toolbox although your code may run substantially faster. In later versions of MATLAB with the *MATLAB Coder* toolbox, this can be done with the following MATLAB script:

	>> compile_mex

## System requirements

System requirements vary depending on which features one intends to use or modify. In short,

- doxygen is required to generate HTML source code documentation
- latex is required to build the file specification PDF
- python is required to run the build tests (i.e. `make check`)

These are optional and not essential to start using Pulseq.

