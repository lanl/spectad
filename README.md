# Speculatively-Parallel Temperature-Accelerated Dynamics (`SpecTAD`)

**Principal SpecTAD Author:**

Richard J. Zamora (2015-2017)

rjzamora@lanl.gov

Theoretical Division,Los Alamos National Laboratory, Los Alamos, NM

**Original TAD Authors:**

Arthur F. Voter and Blas P. Uberuaga

***
This software is open source software available under the BSD-3 license.
 
Copyright (c) 2017, Los Alamos National Security, LLC
All rights reserved.
 
Copyright 2017. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.
 
Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1.       Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2.       Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3.       Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
***


**OFFICIAL IDENTIFYING INFORMATION:**

**Identifying Number**: C17142    
**Software Name**:  SpecTAD       

**DESCRIPTION:**

This SpecTAD code is the state-wise parallel version of the Temperature Accelerated Dynamics (TAD) code (`tad2` & `tad3`) written mostly by Art Voter and Blas Uberuaga (LA-CC-02-05).

SpecTAD improves upon TAD by allowing for *speculative* parallelization, *replica* parallelization, spatial *force-call* parallelization (through an optional `LAMMPS` binding), and a file IO-free synthetic TAD approach (making it **faster**, but more memory hungry than serial TAD).

Beware that this code is still primarily a research code. Therefore, not all features have been tested and debugged and it is **not** user friendly (*yet*).

## Basic Structure of the Code

SpecTAD is a state-wise parallel version of the original TAD code. It currently requires a minimum of **THREE** MPI ranks to run:

* The first rank is the **Manager** rank. The Manager rank constructs the official state-to-state trajectory. It does this by piecing together serial (single-state) TAD runs that were executed on the **Slave** ranks. The Manager is also responsible for assigning work to the slaves (work being a specific state to perform *serial* TAD in).

* The last rank is the **State Manager** (previously called the *File* Manager). This rank keeps a history of all visited states if `irecognize>0` is set in the input file.

* All other ranks are **Slave** processes (although this term is not really used in the source code). The slaves wait for the Manger to assign them a state to perform the TAD algorithm in. When they finish with their assigned state, they just wait (in a loop) for another state to be assigned.

## Getting Started

### Unpack the code:

```
gunzip spectad.tar.gz
tar –xvf spctad.tar
cd spectad
```

Alternatively, using **github.com**:

```
git clone https://github.com/lanl/spectad.git
cd spectad
```

## Building the Code

The non-makefile script to build the code is called **`maktad`**. Right now the SpecTAD code has only been tested with OpenMPI. Beware that using other mpi implementations may result in issues.

This default version of the maktad script will not allow for LAMMPS to be used for force calls (**`maktad-lammps`** is needed to enable LAMMPS calls).

*See below for building the code with LAMMPS.*

### Building the code (Without a LAMMPS binding):

If you have OpenMPI installed and ready to go (with the `mpif90` executible working), just type the following commands (otherwise, you will need to modify the ‘maktad’ file a bit):

```
cd ~/spectad/src/
bash maktad
```

If the `maktad` file does not work for your machine, you probably need to change `compiler="mpif90 "` and `linker="mpif90 "` lines to refer to the correct compiler for your system. To do this, you can either modify one of the existing `machinetype` options, or you can just create a new option (ex: `machinetype="cray"`).

### Setup the `pots/` directory:

When LAMMPS is not used, SpecTAD requires the *clsman* potential format for the native inter-atomic potential implementations. Many eam potentials are included in the `pots/` directory, but they need to be compiled before use. To compile these potentials...

```
cd ~/spectad/pots/.
```

Note: The next line requires that you use the specific fortran compiler setup for your machine. For example, <Fortran compiler> might be *`f77`* or *`gfortran`* on your machine:

```
<Fortran-compiler> ../src/potcong.f -o potcong
./potcong <potcong.dat
```

If you don't want to copy the potential files to the currnt working directory whenever you run SpecTAD, you need to set a `POTSDIR` environment variable. For example:

```
export POTSDIR='/Users/username/spectad/pots/'
```


## Running the Code

To run the code, start with one of the provided examples in the `examples/` directory: 

```
cd ~/spectad/examples/
```

This directory contains some input files to simulate the dynamics of the simple Ag surface trimer, as discussed in [DOI: 10.1146/annurev-chembioeng-080615-033608](http://www.annualreviews.org/doi/full/10.1146/annurev-chembioeng-080615-033608).

The `examples/` directory contains the `ag-cls/` and `ag-lmp/` directories. The first of these directories uses the native *clsman*-based EAM implementation in the code. The second uses LAMMPS to perform force-calls, so the example must provide a *setfl*-formatted EAM potential for Ag. Both versions of the example will use *Synthetic* SpecTAD (with direct vineyard calculations to determine pre-exponential factors) to simulate 1 millisecond at 200K (using T_High=750K). Note that SpecTAD will always use synthetic mode when `irecognize=1` is set in the input file (`trimer.input` for both versions of the example). The code will only use the vineyard prefactor in synthetic mode if `ivineyard=1` (otherwise the prefactor will be estimated using the observed transition rate). 

Choose whether to use the native force-call implementations (`ag-cls/`), or to use LAMMPS (`ag-lmp/`), and go to the corresponding directory...

#### The `ag-lmp` Example:

This example uses `ietype=199` in the input file - This means that the code will **NOT** use your POTSDIR environment variable to find the potential files. Instead, you must have a `script.lmp` file and the `Ag.eam` potential file in the current working directory. (For this example, these files are included in the `ag-lmp/inputs/` directory).  Notice that `script.lmp` includes the potential definition that would normally be written in a lammps input file. When ietype is set to 198 or 199, this file must be used to tell the LAMMPS engine what potential you are using. The two alternatives correspond to **charge** and **atomic** LAMMPS systems, respectively (*warning* ietype 198 and 199 are very new, and may need to be tweaked for potentials that have yet to be tested).

To run this example, copy the contents of this `inputs/` directory into your current wording directory (which should be `~/spectad/ag-lmp/` for this example). Then, the code can be executed using the following commend (with the appropriate mpirun command for your system):

```
mpirun -np <2+NS> ../../src/spectad.exe <trimer.input >trimer.lis &
```
*Where <2+NS> is the number of mpi ranks to use, and NS (>=1) is the nuber of speculative processes you want to use to parallelize your simulation. Note that NS=1 is not really using any speculations, because the single speculation will be exploring the same offical state that Serial TAD would be exploring.*

#### The `ag-cls` Example:

This example uses `ietype=0` in the input file - This means that the code **WILL** use your POTSDIR environment variable to find the EAM Ag potential files. If the code does not run, there is a good chance that your POTSDIR variable is NOT correctly defined…

To run this example, copy the contents of this `inputs/` directory into your current wording directory (which should be `~/spectad/ag-cls/` for this example).  Then, the code can be executed using the following commend (with the appropriate mpirun command for your system):

```
mpirun -np <2+NS> ../../src/spectad.exe <trimer.input >trimer.lis &
```

*Again - Where <2+NS> is the number of mpi ranks to use, and NS (>=1) is the nuber of speculative processes you want to use to parallelize your simulation. Note that NS=1 is not really using any speculations, because the single speculation will be exploring the same offical state that Serial TAD would be exploring.*

### Other Considerations

The TAD/SpecTAD options can be modified in the `trimer.input` file for both veersions of the example. 

The SpecTAD input file should be defined as **`<casename>.input`**, where `<casename>` is determined by the user. By convension, the starting geometry file is consistantly defiled as **`<casename>.start`**, while the output listing is usually chosen to be **`<casename>.lis`**.

For example, if you wanted to manually setup a SpecTAD run, you would do something similar to this:

```
cd ~/new-run-directory/
    *** Populate your case1.input and case1.start files ***
mpirun <path to spectad.exe> -np 3 <case1.input >case1.lis &
```

Note that we are asking SpecTAD to use the minimum number of MPI ranks here (3). Only **one** of these MPI ranks will be used to perform TAD work.

Note that the **`itad`** value in the input file should be set to 100. Also, note that **`ivineyard=1`** is usually efficient only for simple systems (`ivineyard=0` is a more robust setting).

If you want to take advantge of parallel force calls or replica-based parallelism of each speculative process, you can change the first two lines in the input file. For example,

```
2   nprocs_force - Number of cores for force calls
10  nprocs_rep - Number of ParRep Slaves for each SpecTAD process
...
```

will assign 10 replica processes to EACH speculative process, and each replica process will use 2 MPI ranks to perform a parallel force call. This means that every speculative process will need 20 MPI ranks. In this case, the minumum number of MPI ranks needed to run the code with 1 speculative process (`nprocs_spec=1`) is: **`nprocs_force x nprocs_rep x nprocs_spec + 2 = 22`**

## The SpecTAD Output

Much of the SpecTAD output is similar to the original TAD code. One of the main differences is that SpecTAD summarizes the official trajectory in the **`spawnlist.dat`** file. Every minute or so, during runtime, `spawnlist.dat` will be updated with the current *official* trajectory. This file lists information for every assigned state visit (**spawn**) explored by slave process that is on the official trajectory (has not been *pruned* by an ancestor spawn). The first five lines of the output would look something like this:

```
1  1  2.990892E-07  8.018179E+01  1.590439E+03  2  -8.156072E+02  2.050352E-01
2  2  3.254350E-07  1.352604E+02  1.025854E+03  3  -8.156072E+02  2.050190E-01
3  3  3.270469E-07  1.584866E+02  8.798512E+02  2  -8.156072E+02  2.050352E-01
4  4  4.342972E-07  2.303036E+02  8.040409E+02  3  -8.156072E+02  2.050230E-01
5  5  4.498858E-07  2.819405E+02  6.803568E+02  4  -8.156072E+02  2.050225E-01
```

Here, the Columns define the following information:

1. Total number of state visits (or total transitions) on the official low-temperature trajectory at the end of this spawn
2. The **Spawn ID** for this specific spawn
3. Global low-temperature time for the (final) winning transition of this spawn. Not that *final* refers to the fact that a single spawn can skip through states (using kMC-based *Synthetic* TAD priciples) if it recognizes that the assigned state (or subsequently kMC-chosen states) does not require any additional high-temperature MD. This is also the reason that many state visits can be *skipped* in **column 1**
4. The human wall clock time (WCT), with respect to the start of the global simulation, when this spawn finished
5. The estimated **Boost** achieved, over MD, for the entire simulation up to this spawn
6. The official **State ID** for the last state visited by this spawn
7. The total energy of the last state visited by this spawn
8. The final escape barrier accepted by this spawn

Note that the official states (as referenced in `spawnlist.dat`) are written out as `state.<State ID>.dat` during run time.


## Compiling with LAMMPS
Before you can compile SpecTAD with LAMMPS, you need to build LAMMPS as a shared library in a separate directory (see section 2.5.2 [here](http://lammps.sandia.gov/doc/Section_start.html#start-5)).

This creates the LAMMPS library, but does not allow you to interface with a Fortran code without a bit more work. Since SpecTAD is written in Fortran and LAMMPS is written in C++, you also need to follow the instructions to build the Fortran wrapper in the **`/examples/COUPLE/fortran2`** directory within LAMMPS.

Note that, as an example, the contents of the **`fortran2/`** directory used to test the code is included in the **`lammps-fortran-example/`** directory.

After the fortran wrapper is built, you should set **`LAMMPS_ROOT`**, **`LAMMPS_SRC`** and **`LAMMPS_FORT`** correctly in the **`maktad-lammps`** script.

You should also set your **`LD_LIBRARY_PATH`** (**`DYLD_LIBRARY_PATH`** on a Mac) to include the **`LAMMPS_SRC`** and **`LAMMPS_FORT`** paths.

Once your **`maktad-lammps`** script and environment variables are set correctly, you should be able to build the code by typing:

```
cd ~/spectad/src/
bash maktad-lammps
```

Once the code is built with LAMMPS, you can use a LAMMPS-based potential to perform a run as long as the specific potential type is taken into account in the **mod_lammps.f** source code file in SpecTAD.

### Example SpecTAD+LAMMPS build steps on Mac

In some designated top directory, clone lammps:

``` 
git clone https://github.com/lammps/lammps.git
cd lammps
git checkout stable
```
 
Install necessary packages:

```
cd src
make yes-manybody
...
```
 
Modify the makefile at `~/lammps/src/MAKE/MACHINES/Makefile.mac_mpi` to agree with your system and compile the code as a shared library:

```
make mac_mpi
make mac_mpi mode=shlib
```
 
Copy the library files to a directory called `lib/` that should live in your home directory. If you don’t have a lib directory there – then make one first. Then:

```
cp liblammps.so ~/lib/.
cp liblammps_mac_mpi.so ~/lib/.
```
 
Go back to the designated top directory, and unpack spectad:

```
git clone https://github.com/lanl/spectad.git
cd spectad
```
 
Go back to `/example/COUPLE/fortran2/` directory in lammps, and copy over the makefile distributed in SpecTAD:

```
cd lammps/example/COUPLE/fortran2
cp ~/SpecTAD-modes/lammps-fortran-example/makefile .
```
 
Change the makefile to have the correct compilers for your system, and build the fortran wrapper for the lammps interface:

```
make
```
 
Copy the `liblammps_fortran.so` file to the `~/lib/` directory:

```
cp liblammps_fortran.so ~/lib/.
```
 
Modify the maktad-lammps script to agree with your system, and build:
 
* Make sure the `LAMMPS_ROOT` variable is set correctly for your lammps location
* Make sure the compiler and linker are the correct mpif90 compiler for your system
* Type: `bash maktad-lammps`
 


*More documentation is definitely needed/missing (Detailed output listings, adding LAMMPS capabilities, other input file considerations, atom deposition, etc..)*
