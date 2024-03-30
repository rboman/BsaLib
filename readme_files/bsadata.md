# `bsa.bsadata` input file structure

This file is required by the `BSA` standalone executable to gather information 
related to `BsaLib` specific settings. 

There are 5 sections (**WARNING**: to be provided in the given order!):

1. [-GENERAL](#general)
2. [-CLASSIC](#classic)
3. [-MESHER](#mesher)
4. [-DIRECTIONS](#directions)
5. [-TURBULENCE](#turbulence)



## General

Currently supported general settings (in order):

#### Analysis type ID (`int`)

Valid values: `1` (classic type); `2` (mesher type); `3` (both). 

Default: `1`.

#### Version ID (`int`) [UNUSED]


#### PSD (Power Spectral Density) scaling convetion ID (`int`)

Valid values: `1` (circular frequencies convention, $\sigma^2 = \int_{-\infty}^{\infty} S(\omega) d\omega$); 
`2` (frequency convention, $\sigma^2 = \int_{0}^{\infty} S(f) df$). 

Default: `1`.

#### PSD computation switch (`int`)

Valid values: `0` (OFF, no PSDs computed); `1` (ON). 

Default: `1`.

#### Bispectra computation switch (`int`)

Valid values: `0` (OFF, no bispectra computed); `1` (ON). 

Default: `1`.

#### Only diagonal elements computation switch (`int`)

Controls whether all tensors elements are computed, 
or computation limited to tensors elements whose indices are equal in all dimensions (i.e. $a_{iii}$). 

Valid values: `0` (OFF, all elements computed); `1` (ON, only diagonal elements computed). 

Default: `0`. 

**NOTE**: indeed, `1` results in a much faster computation, specially for Bispectra computation. 
It is however less precise due to the neglecting of outer-diagonal elements, which might play a crucial role 
in the final statistical moments estimates. This is much better explained in the [companion paper](https://doi.org/10.1016/j.jweia.2023.105576).

### [EXPERIMENTAL] Bispectra spatial symmetry value (`int`) 

Controls how much information is implicitly assumed 
leveraging spatial (in-plane) symmetries of real part of bispectra. 

Valid values: `0` (FULL, all information computed); `2` (HALF); `4` (FOURTH). 

Default: `0`. 

**NOTE**: at one point, default should become `2`, since this would be the best solution.

#### [EXPERIMENTAL] Tensor symmetry switch (`int`) 

Allows to control whether symmetrical elements of a spectra tensor 
(PSD or Bispectra) are computed (i.e. $a_{ijk}$ computed only once for all possible permutations of 
indices $(i,j,k)$.) or not. 

Valid values: `0` (OFF, all tensor elements are computed); `1` (ON, only unique 
symmetrical tensor elements are computed). 

Default: `0`. 

**NOTE**: at one points, default should become `1`. 

**NOTE2**: this flag has only sense if diagonal elements only flag is OFF.

#### Test flag (`int`)

If ON, disables some internal checks, specially at the level of discretisation refinement 
compared to suggested values. 

Valid values: `0` (OFF, check enabled); `1` (ON, checks disabled). 

Default: `0`.



## Classic

The classic group controls settings regarding the "classic" approach, which includes 
"conventional" spectral and bispectral analysis implementations.

There are only two main values to be set:

#### Number of discretisation frequencies $\tt N_{freqs}$ (`int`)

To be considered from `0` to $f_{max} \text{ [Hz]}$ (or $\omega_{max}\ [\frac{rad}{s}]$ if 
circular frequency convetion used).

Valid values: `> 0`.

Default: NONE, a value must be provided.


#### Delta frequency (refinement) $\Delta f \text{ [Hz]}$ (`double`)

After n. of frequencies is give, $\Delta f$ specifies the (regular) spacing between each frequency.

Valid values: `> 0.`.

Default: NONE, a value must be provided.



## Mesher


#### POD (Proper Orthogonal Decomposition) Switch (`int`)

Controls whether POD is used for decomposing base wind flow field. See [reference](https://doi.org/10.1016/j.jweia.2023.105576) 
for more details.

Valid values: `0` (OFF, no POD); `1` (ON).

Default: `1`.

#### Background zone base refinement (`int`)

Specifies the base number of meshing points (per side) to disretise the rectangular zone 
covering the background peak at the origin $(0, 0)$.

Valid values: `> 0`.

Default: NONE, must be provided.

#### Background zone area extension (`int`)

Specifies the integer multiplier of the base background zone extension. 

Valid values: `> 0`.

Default: `1` (NO extension).

**NOTE**: at some point, this parameter must be changed to be a floating point number.

#### Peak zone area extension (`int`)

Specifies the integer multiplier of a peak zone extension, where a peak zone is a meshing zone 
covering a part of the 2D space where a (bi-)resonant peak is located.

Valid values: `> 0`.

Default: `1` (NO extension).

**NOTE**: at some point, this parameter must be changed to be a floating point number.

#### Max area extension (`int`)

Specifies the integer multiplier of the max area extension. 
Internaly, a maximum area (of interest) limit is determined. This factor allows to extend this limit 
up to the desired value.

Valid values: `> 0`.

Default: `1` (NO extension).

**NOTE**: at some point, this parameter must be changed to be a floating point number.


#### Full coverage switch [DEPRECATED]


#### Dump modal info switch (`int`)

Controls whether to include modal info in dump file (used in post-processing tasks), or not.

Valid values: `0` (NO, do not include); `1` include.

Default: `1`.


## Directions

Ignore. See [example](#example-of-working-`bsa.bsadata`-file).



## Turbulence (`int - int[]`)

Specifies the (3) spatial turbulent component to be accounted in the determination of the 
stochastic wind loading. 

First entry is an `int` specifying the effective number of turbulent components.

Following, a list (one entry per line) of component indices, where each index should be in the set 
${1, 2, 3}$ ($(x, y, z)$ components in an Euclidean 3D space).

**NOTE**: in the [example below](#example-of-working-`bsa.bsadata`-file), we declare to consider 
2 turbulent components, component `1` (x) and `3` (z).




### Example of working `bsa.bsadata` file

```text
-GENERAL:
   2
   0
   0
   1
   1
   0
   0
   0
   1
-CLASSIC:
   500
   0.010
-MESHER:
   1
   50
   2
   2
   2
   1
   1
-DIRECTIONS:
   1
   1
-TURBULENCE:
   2
   1
   3
```






