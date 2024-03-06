# `bsa.extdata` input file structure

The other input file read by `BsaLib` is named `bsa.extdata`, and as the name suggests 
should provide all the data that is "external" to `BsaLib` itself (i.e. do not affect its 
internal settings directly). 
It is however data needed for the proper functioning of `BsaLib`.

This is a binary (Fortran `unformatted`) file, where:

- Each integer entry must be compatible to the `int32` kind of the `iso_fortran_env` Fortran 
intrinsic module;
- Each floating point entry must be `real64` compatible (from the same ISO module).


This data mostly concerns:

1. Data related to the [structure/system](#structural-data) under investigation
2. Data related to the characterisation of the [wind action](#wind-data)


**NOTE**: in the following, whenever an array of size $\tt N$ should be read from 
this file, it is intrinsically meant to be read like this:

C:
```c
type *data = (type *)malloc(N * sizeof(type));
// check allocation OK

// loop version
for (unsigned int i = 0; i < N; ++i) {
    fread( (data + i), sizeof(type), 1, fhandle );
}

// all-at-once version
fread( data, sizeof(type), N, fhandle );
```

Fortran:
```Fortran
allocate(data(N), stat=istat)
! check allocation OK

! (implicit) loop version
read(fhandle) (data(i), i = 1, N)

! all-at-once version
read(fhandle) data
```




## Structural data

Structural data must be ordered as such.

#### Number of nodes (total) $\tt NN$ (`int`)

Specifies the total number of structural nodes with which the model is discretised.

#### Number of degrees-of-freedom per node $\tt NLIB$ (`int`)

Specifies the total number of degrees of freedom available at each node.

#### Number of nodes loaded $\tt NNL$ (`int`)

For optimal performances, it is good to separate the total number of nodes and 
the actual number of nodes that are under the considered action (wind, wave, etc.) 
for which a Bispectral analysis is desired.

This separation allows for both better performances and less memory footprint.

#### Number of degrees-of-freedom (per node) loaded (`int`)

Same concept as for number of loaded nodes.

**DEPRECATED**.

#### List of loaded nodes indexes $\tt NL$ (`int[]`)

Specifies the array of $\tt NNL$ indexes (between $1$ and $\tt NN$) of nodes 
that are effectively loaded.

#### List of loaded nodal degrees-of-freedom (`int[]`)

Specifies the array of indexes of nodal degrees-of-freedom that are effectively loaded.

**DEPRECATED**.

#### Node coordinates (all) (`double[][]`)

2D-array of nodal coordinates.

Coordinates should be provided as a series of $(x, y, z)$ coordinates, for each node, 
from first to last (in order).





## Wind data


#### Vertical wind profile ID (`int`)

#### Wind PSD type ID (`int`)

#### Global vertical axis ID (`int`)

**DEPRECATED**

#### Degree of wind model transformation (`int`)

#### Air density (`double`)

#### Wind-to-Global rotation matrix (`double[3][3]`)

Rotation matrix that translates WRS (Wind Reference System) euclidean coordinates into 
GRS (Global Reference System) ones.

#### Number of wind zones 



