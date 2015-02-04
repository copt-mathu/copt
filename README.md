
<h1>
COPT: A C++ optimization open library
<br/>
<img src = "https://cloud.githubusercontent.com/assets/6701845/5139813/cee4c87a-719d-11e4-9b90-5d1826f2bc39.png", width = "100px"/>
</h1>
<h3>Copyright (C) MathU <img src = "https://cloud.githubusercontent.com/assets/6701845/5139913/2ac461ae-719f-11e4-8d72-3cdea63566b3.png", width = "60px" />
</h3> 

<p>
COPT is a C++ based optimization open library. 
The very simple idea is to provide a complete computing library for basic optimization algorithms. 
It is powered by MathU organization in University of Science and Technology of China (USTC). 
The leaders are Zhouwang Yang, associated professor in USTC and Ruimin Wang, PhD candidate in USTC.
</p>

<p>
The library contains basic mathematical types like vector, matrix, function and solver. 
Matlab-like element access and assignment are allowed. 
It is quite simple to claim mathematical types after fixing basic traits.
</p>

-----------------------------------------------------------------
<h2>Development Guide</h2>
<p>This is the introduction and development for C++ based open source library COPT focusing on optimization problems. Current version of COPT is 0.0.1 and the document is only for developers right now. This document tells  you how to install COPT on your personal computer and further develop optimization algorithms based on COPT. </p>

<h3>Installation Guide</h3>

<p>COPT currently depends on three third party libraries as BLAS, LAPACK and SuiteSparse. COPT is developed on OS X system and it currently only support OS X and Linux operation system. MinGW on windows is another choice but not a good one. Now installation of COPT on OS X is introduced. BLAS and LAPACK are built-in libraries for OS X as default. Thus the only thing you have to do is to install SuiteSparse. At first, command line tool must be installed which can be easily done by downloading it in AppStore. Then you download the source code of SuiteSparse on http://faculty.cse.tamu.edu/davis/suitesparse.html. After that, remove the  file “SuiteSparse_config.mk” in the sub-dictionary ‘SuiteSparse_config’ and rename  the file “SuiteSparse_config_mac.mk” to “SuiteSparse_config.mk”. Then just type ‘make’ and ‘sudo make install’ in command line tool. The installation on Linux is not complex either. Just type the following commands in terminal: “sudo apt-get libblas-dev”, “sudo apt-get liblapack-dev” and “sudo apt-get libsuitesparse-dev”. </p>

<p>Now the dependency has been installed and you can clone or fork our codes on https://github.com/fromradio/COptdev. Then you can use the library now and you can type ‘make all’ to generate all examples written by us. </p>

<h2>Development Guide</h2>

<p>COPT is a C++ based library which uses a lot of C++ features. The library mainly contains the following parts: the basic types for describing mathematical objects, the optimization problems and the optimization algorithms. To use COPT, one only needs to add the header file “Header” under folder include for optimization part and “IO” for input/output.</p>

<h3>1. Mathematical types</h3>

<p>The mathematical types of COPT are all contained in ‘KernelTrait’ which takes the type of scalar and integer as template. The default templates are set as ‘double’ and ‘int’. However, one can modify the kernel he or she wants to use for instance one can easily define a complex kernel via the following code:</p>

->typedef COPT::KernelTrait<std::complex<double>,int> kernel;<-

