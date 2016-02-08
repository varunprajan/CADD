# CADD-Code

## Overview

Performs concurrent multi-scale mechanics simulations in 2D, by implementing the Coupled Atomistic and Discrete Dislocation (CADD) algorithm (see Shilkrot, Miller, and Curtin, JMPS, 2004). Written in modern Fortran (Fortran 2008+), with pre- and post- processing capabilities offered by various Python modules. Is also written in a modular fashion so that other types of simulations can be performed: to wit, 1) molecular dynamics ('atomistic'), 2) finite elements ('fe'); 3) discrete dislocation dynamics ('dd'); and 4) finite elements coupled to atomistics with no discrete dislocations ('cadd_nodisl').

## Layout

1) Fortran. The Fortran source code is in the folder "Fortran". THe files with prefix "mod_" are module files I have written myself. The files without prefix "mod_" are program files I have written myself, or routines from the Harwell Subroutine Library for linear algebra. I have used MA57, which implements the multifrontal method for solving Ax = b where A is sparse, symmetric, and not necessarily positive definite. This is used to decompose/solve Kx = f in finite elements, where K is linear, and not necessarily positive definite because constraints are imposed using Lagrange multipliers. If this code is used to generate results for a scientific paper, HSL should be given appropriate credit.

2) Python. The Python module files are in the submodule "PythonModules". Their use is illustrated in the IPYthon notebook files in the folder "Python", where each notebook is for a different simulation. (I would suggest starting with "Atomistic_2.ipynb", which is for an atomistic crack simulation, and "Source_with_Obstacles.ipynb", which is for a simple DD dipole simulation.) In short, the modules "cadddatadump" and "cadd_plot" are for visualizing the dump files, and the module "cadd_io" is for various input/output tasks, and the module "cadd_main" has the main routines (in an object with all the simulation information) that drives everything. I have also used the modules "compile_cadd" and "run_cadd" to automate the compilation and running of simulations on the EPFL Linux clusters.

3) Tests. Contains documentation for unit tests and tests of broader functionality of the code. Some of this (particularly the unit tests) is quite messy.

4) References. Contains useful (and not-so-useful) scientific papers for CADD, DD, etc., as well as guidelines for best practices for Fortran code.

## Notes on code philosophy and goals

The previous CADD code suffered from numerous flaws, including 1) memory leaks; 2) an almost complete lack of documentation; 3) spaghetti-like code; and 4) algorithmic deficiencies, particularly with dislocation detection; and (most importantly) 5) not working properly. In my view, these deficiencies stemmed from several sources: A) the multiple authors with inconsistent coding styles; B) the inherent limitations of Fortran 77; C) the lack of forethought in organizing the code in a modular fashion (which reflects the layer-by-layer and haphazard nature of its construction); D) the lack of documentation or use of version control to track changes; and E) the lack of an explicit statement of the actual algorithm for CADD (a description in a scientific paper does not suffice).

I have attempted to rectify these flaws by writing the code from scratch, using a uniform coding style and modern Fortran (Fortran 2008). All subroutines and functions have been unit tested (to the greatest extent possible); public and private variables and functions/subroutines ensure that namespace pollution does not occur; every function is carefully documented; code reuse is maximized by using short, single-responsibility subroutines; and the best practices summarized in the folder "Fortran Best Practices" are followed for the most part. (Two exceptions are that I use fixed format Fortran, which seems more visually appealing, and I haven't yet added intent declaration for variables, although this should be fairly easy.)

(Related note: When adding functionality to this code, please follow these style conventions!)

Furthermore, and equally importantly, I have taken away responsibility from the Fortran code for tasks that are difficult to implement in Fortran: particularly, sophisticated input/output. As a consequence, the Fortran input files are not easily readable by humans, and should not be generated directly. Instead, the Python pre-processing module should be used to generate them.

I have also attempted to describe the CADD algorithm in the text file "algorithm.txt."

There are three possible sacrifices made in this effort to redo CADD from scratch: 1) functionality (the code as currently constituted is minimal, working only for 2D lattices, single-material continuum domains, etc.); 2) speed (since it is unclear which sections of the code are actually slow, I have opted against "premature optimization"); and 3) correctness (the code has not survived the rigors of 10+ years of testing, although since the previous code didn't work, perhaps this testing wasn't so rigorous after all). I have included notes where possible indicating where and how the code's functionality needs to be extended. One modification that may be very important is implementing the multipole method, which is significantly faster than the direct approach for calculating dislocation fields for moderate numbers of dislocations (~a few hundred).

The code would probably benefit (eventually) from being rewritten in a truly object-oriented language like C++ (this was the same evolution that LAMMPS underwent --- from Fortran to C++). I have left this task to a future programmer who knows C++ better than I do.
