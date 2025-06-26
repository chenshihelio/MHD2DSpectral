This is the 2D pseudo-spectral MHD code for the course UCLA M263A - Solar System Magnetohydrodynamics by Prof. Vassilis Angelopoulos.

## Setup

### Dependencies

First, install the Conda Package Manager. Miniconda is enough: (https://www.anaconda.com/docs/getting-started/miniconda/install).

You can install the required packages by running the following command:

```sh
conda env create --file environment.yml
```

or with pixi:

```sh
pixi install
```

To activate the environment, run:

```sh
conda activate laps
```

or with pixi:

```sh
eval "$(pixi shell-hook -e dev)"
```

### Installation

To install the package, run the following command:

```sh
cd src && make install && make clean
```

Note that `pylasp` (Python interface to LAPS) is also installed as a package when using `pixi install`.

## Run

In order to run a new simulation:

1. create a new directory, where the simulation will be run

2. make sure the `lasp` executable is either copied into this directory or in your PATH environment variable

3. add an inputs file `mhd.input`

4. run

Run the executable, e.g. with MPI:

```sh
mpirun -np 4 laps
```

or using the `run` command in `justfile`:

```sh
just run
```

This will also create a new directory `diags` where the output files will be stored.

<hr>

Copyright 2024 Chen Shi

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

If you use the code to generate any publications, you should cite the following articles 

@article{shi11laps, title={LAPS: An MPI-parallelized 3D pseudo-spectral Hall-MHD simulation code incorporating the expanding box model}, author={Shi, Chen and Tenerani, Anna and Rappazzo, Antonio Franco and Velli, Marco}, journal={Frontiers in Astronomy and Space Sciences}, volume={11}, pages={1412905}, publisher={Frontiers} }

@article{shi2020propagation, title={Propagation of Alfv{'e}n waves in the expanding solar wind with the fast--slow stream interaction}, author={Shi, Chen and Velli, Marco and Tenerani, Anna and Rappazzo, Franco and R{'e}ville, Victor}, journal={The Astrophysical Journal}, volume={888}, number={2}, pages={68}, year={2020}, publisher={IOP Publishing} }
