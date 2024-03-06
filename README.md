<!-- ![BSA logo](./resources/images/BSA_logo_extended.PNG "BSA logo") -->

# BsaLib - Bispectral Stochastic Analysis Library

**BsaLib** is library/framework for the Bispectral Stochastic Analysis (BSA) of linear systems, 
under non-Gaussian stationary random excitations.

**NOTE**: currently, only wind action is included within, but other phenomena (waves for instance) 
can be easily integrated. See [further developments](#what's-missing?-further-developments).



# Code structure

There are two parts in this repository: (i) `BsaLib` (library) and (ii) `BSA` (executable).

## `BsaLib` - core library

`BsaLib` is the main core of this repository (in `./src/`). 
It consists of the main library and its API to which anyone could link to and interact with.
To use `BsaLib`, simply import the [main API module](./src/BsaLib.f90):

```Fortran
program test
    use, non_intrinsic :: BsaLib
    implicit none (type, external)
    ! your declarations here

    ! your logic here

    ! initialise BsaLib
    call bsa_Init()

    ! set BsaLib internal state through its API procedures

    ! once done, run BsaLib
    call bsa_Run( ... args ... )

    ! finally, release BSA memory
    call bsa_Finalise()
end program
```

Being designed as a *plug-in* library, it needs the hosting program/library 
to provide some data needed by `BsaLib` in order to function properly.
For example, currently, all structural modal data is required to be given, 
such as Mass, Damping and Stiffness matrices ($M^*, C^*, K^*$), as well as modal damping ratios ($\xi$).
However, this is usually data that almost all Finite Element Software are capable of computing.
For details about all data required by `BsaLib`, please refer to [next section](#`bsa`---executable-program) 
and [`bsa.extdata` section](./readme_files/extdata.md).



## `BSA` - executable program

As a side project package, a single-source executable file is provided under `./app/`. 
It emulates what one would normally do when using `BsaLib` as a *plug-in* for its own 
library/program. 
On the other hand, this program is thought and provided for all those interested in using 
`BsaLib` but not having any hosting library/program. 
Nonetheless, even if this program is provided, the user would yet need to provide the data 
for `BsaLib` to function properly. If any of this data is not provided, the `BsaLib` runtime check 
would detect it and abort correct logic flow.

For that, the provided executable program relies on reading two input files:

1. `BsaLib` related settings (formatted file, named `bsa.bsadata`). 
For details, read the [dedicated section](./readme_files/bsadata.md).
2. External data file (named `bsa.extdata`), in binary format, containing 8-byte floating point records 
(`real64` of the `iso_fortran_env` compiler intrinsic module). 
For full details, read the [external data section](./readme_files/extdata.md).




# Related published works

This work includes some novel work, on a mathematical and numerical level.
Please refer to following paper(s) for detailed information:

- [Non-Gaussian buffeting analysis of large structures by means of a Proper Orthogonal Decomposition](https://doi.org/10.1016/j.jweia.2023.105576)



# Known Issues

There is one main known issue in the current version. 

> Using `OpenMP` parallelisation in the second Post-meshing phase, 
> execution time is higher compared to serialised version. 
> This is due to the necessary synchronisation between threads when accessing 
> shared file I/O when reading each zone's dumped data, causing the `critical` 
> section to be the main bottleneck in this part.
> For a proper use of `OpenMP` parallelisation in Post-meshing phase, a fundamental
> rethinking of the algorithmic structure needs to be done.
> A first, temporary, possible solution, could be the usage of thread-level private
> I/O with dedicated units, avoiding the need to synchronise when reading back 
> information in post-meshing phase. However, this might soon become inelegant solution 
> when the number of threads would start increasing considerably.
> For this, as a temporary solution, a conditional compilation flag 
> (`BSA_USE_POST_MESH_OMP`) can be used to control effective use of 
> `OpenMP` parallelisation in Post-meshing phase, **disabled by default**.



# What's missing? Further developments

- [ ] Integrate models for other non-Gaussian actions (waves, for instance)
- [ ] Complete full support for spatial (in-plane) symmetries of a real-valued bispectrum
- [ ] Compute nodal correlation internally (don't require it as user data)
- [ ] Add support for $\mathtt{Mesher}$ zones' interest modes
- [ ] Add functionality to generate spectra from time series
- [ ] Add a local caching system
- [ ] Add `MPI` support (for running `BsaLib` on multi-node clusters)
- [ ] Add ability to export using formats of most common Visualisation (Paraview, for instance). 
Better, in a more sustainable way, add the possibility to let the user provide its own desired 
exporting function, so that `BsaLib` is not tight to any specific 
exporting format.
- [ ] Integrate a built-in ad-hoc bispectrum post-processing Visualiser (using [Vulkan](https://www.vulkan.org/), for optimal performances)

