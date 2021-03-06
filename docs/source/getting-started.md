# Getting Started

## Dependencies

The IPE source code makes use of both HDF5 and NetCDF libraries. Additionally, parallelism is supported using MPI (via OpenMPI or MPICH). In order to build and run IPE, you will need to have HDF5 and NetCDF installed on your system. For parallel builds of IPE, you will need to have parallel HDF5 and parallel NetCDF installed.

## Obtaining Source Code

The official repository for WAM-IPE is hosted on GitHub: <https://github.com/CU-SWQU/GSMWAM-IPE>.

(input_decks)=

## Input Decks

IPE depends on a set of files that are needed for the accompanying empirical models to specify the model grid, model runtime parameters, and the initial conditions. These files, along with a description and where it is used in the source code is listed in the [Input Files](input_files) section.

The input decks can be downloaded using the provided python script under **IPELIB/scripts/download_ipe-data.py**
This script uses wget to download data over https.

```{code-block} bash
cd IPELIB/run
python ../scripts/download_ipe-data.py
```

````{admonition} CAUTION
:class: warning

If you are on theia, the python script will not work. Instead, you will need to copy files from **/scratch3/NCEPDEV/swpc/refactored_ipe_input_decks/**
```{code-block} bash
cp /scratch3/NCEPDEV/swpc/noscrub/refactored_ipe_input_decks/* IPELIB/run/
```
````

(input_files)=

## Input Files

```{list-table}
:header-rows: 1
:widths: 10 20 30

* - File Name
  - Subroutine
  - Description
* - gd2qd.dat
  - 
  - 
* - global_idea_coeff_hflux.dat
  - 
  - 
* - global_idea_wei96.cofcnts
  - 
  - 
* - wei96.cofcnts
  - 
  - 
* - hwm123114.bin 
  - msise00/hwm14.f90::initqwm
  - Binary file ( little-endian containing a mix of 4-byte integers and double precision ) containing data needed for the 2014 Horizontal Wind Model
* - dwm07b104i.dat
  - msise00/hwm14.f90::initdwm
  - Binary file ( little-endian containing a mix of 4-byte integers and double precision ) containing data needed for the 2014 Disturbance Wind Model
* - ionprof
  - IPE_Plasma_Class.F90::Auroral_Precipitation
  - ASCII file containing values used for calculating ionization rates.
* - tiros_spectra
  - IPE_Plasma_Class.F90::Auroral_Precipitation
  - ASCII file containing values used for calculating ionization rates.
* - IPE_Grid.nc
  - IPE_Grid_Class.F90::Read_NetCDF_IPE_Grid
  - NetCDF file containing IPE Grid information and mapping weights onto geographic grid
* - IPE.inp
  - IPE_Model_Parameters.F90
  - Namelist file used for specifying IPE parameters, including timestep, file IO frequency, and forcing parameters
```

(ipe_build)=

## Building IPE

IPE uses autoconf as it’s make system. To build the code, execute the following from the **IPELIB/** directory, replacing **/path/to/install** with the path where you would like to install IPE binaries and libraries.

```{code-block} bash
INSTALL_DIR=/path/to/install
./configure --prefix=${INSTALL_DIR}
make install
```

This process will create the following directory structure.

```{code-block} bash
${INSTALL_DIR} 
│
└───bin/
│   │   ipe
│   │   ipe_postprocess
│   │   eregrid
│   
└───include/
│   
└───lib/
```

````{note}
If you are working on Theia, you can use the provided **build_ipe_from_scratch.csh** script to build the code.
```{code-block} bash
cd IPELIB/
./build_ipe_from_scratch.csh
```
````

**bin/** contains binaries (programs) for running IPE (standalone), a post-processor to generate diagnostics on a geographic grid, and a program for regridding electric field data from an external model (Geospace or OpenGGCM) onto the IPE Grid.

**lib/** contains library and object files from the compilation process. This directory should be added to the linker when coupling IPE with another code or in developing additional tools using the IPE API.
`-L/${INSTALL_DIR}/lib -lipe`

**include/** contains the **.mod** files, generated by the Fortran source code compiliation. This directory should be passed as an includes directory at compile time for coupling IPE with another code or in developing additional tools using the IPE API.
`-I/${INSTALL_DIR}/include`

## Setting the model input parameters

Modify **IPE.inp**

```{list-table}
:header-rows: 1

* - Parameter Name
  - Namelist
* - NetCDF_Grid_File
  - SpaceManagement
* - NLP
  - SpaceManagement
* - NMP
  - SpaceManagement
* - NPTS2D
  - SpaceManagement
* - NFluxTube
  - SpaceManagement
* - Time_Step
  - TimeStepping
* - Start_Time
  - TimeStepping
* - End_Time
  - TimeStepping
* - Initial_TimeStamp
  - TimeStepping
* - Solar_Forcing_Time_Step
  - Forcing
* - F107_kp_size
  - Forcing
* - F107_kp_interval
  - Forcing
* - F107_kp_skip_size
  - Forcing
```

## Running IPE

To run IPE, it's important that you've [downloaded the input decks
necessary to run IPE](input_decks). The input
decks should be placed in a directory where you wish to run the IPE
source code ( e.g. GSMWAM-IPE/IPELIB/run )

For a serial run, all you need to do is invoke the executable

```{code-block} bash
${INSTALL_DIR}/bin/ipe
```

in the directory where all of the input decks are located.

To run in parallel, make sure that the appropriate MPI tools are in your
search path and call the executable with

```{code-block} bash
mpirun -np 20 ${INSTALL_DIR}/bin/ipe
```

to run the model with 20 ranks. Note that, you can change the number of
ranks to any value between 1 and 40.

### On Theia

Underneath the **IPELIB/run/** directory, batch submission scripts are
provided for running the serial and parallel versions through a PBS job
scheduler. These scripts can be modified to change the submission queue,
the wall time, and the number of ranks.

If you are on theia, you can use the provided job submission scripts to
run the code. For a parallel run

```{code-block} bash
cd ipelib/run
qsub ipe_parallel.pbs
```

For serial,

```{code-block} bash
cd ipelib/run
qsub ipe_serial.pbs
```

For testing both serial and parallel, and comparing the results

```{code-block} bash
cd ipelib/run
qsub ipe_test.pbs
```
