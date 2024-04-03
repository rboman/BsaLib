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
license: LGPL v3
print_creation_date: false
creation_date: %Y-%m-%d %H:%M %z
dbg: true
extra_mods: iso_fortran_env:https://fortranwiki.org/fortran/show/iso_fortran_env
md_extensions: markdown.extensions.toc
---

--------------------

[TOC]

{!README.md!}
