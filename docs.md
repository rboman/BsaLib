---
project: BsaLib
Summary: BsaLib, a Modern Fortran Library for the Bispectral Stochastic Analysis.
author: Michele Esposito Marzino
email: michele.espositomarzino@gmail.com
source: false
src_dir: ./src
output_dir: ./doc
exclude_dir: src/*
preprocess: true
preprocessor: ifort /E /I:.\src\CONSTANTS\
predocmark: >
docmark_alt: #
predocmark_alt: <
display: public
         protected
graph: false
coloured_edges: true
search: true
sort: alpha
license: gfdl
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
dbg: true
extra_mods: iso_fortran_env:https://fortranwiki.org/fortran/show/iso_fortran_env
md_extensions: markdown.extensions.toc
---

--------------------

[TOC]

Brief Description
-----------------

`BsaLib`, a Modern Fortran Library for the Bispectral Stochastic Analysis 
of structures under non-Gaussian stationary random actions.

License
-------

`BsaLib` is release under the **GNU General Public License**.
Visit the [GPL official website](https://www.gnu.org/licenses/gpl-3.0.html) for more information.


Related Scientific Publications
-------------------------------

1. [Non-Gaussian buffeting analysis of large structures by means of a Proper Orthogonal Decomposition](https://www.sciencedirect.com/science/article/abs/pii/S0167610523002799?via%3Dihub)
2. [A multiple timescale approach of bispectral correlation](https://www.sciencedirect.com/science/article/abs/pii/S0167610522003786?via%3Dihub)
3. [On the background and biresonant components of the random response of single degree-of-freedom systems under non-Gaussian random loading](https://www.sciencedirect.com/science/article/abs/pii/S0141029611001507)
