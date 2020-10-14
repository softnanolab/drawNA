# LAMMPD installation

## Summary

LAMMPS is our preffered simulation software due its rigourous usability, and track record as a molecular dynamics simulation software. For simulating DNA, we use the oxDNA2 potential, on which the oxDNA software is based.

## Install

To install, we clone the repository from the lammps group repository. It is important to note that there is a wealth of information included in the documentation for LAMMPS which can be found [here](https://lammps.sandia.gov/doc/Manual.html)

```bash
$ git clone https://github.com/lammps/lammps.git
$ cd lammps
$ cd src
$ make yes-user-cgdna yes-molecule yes-asphere yes-kspace
$ make mpi
```

## Test

Once we have correctly installed LAMMPS, we can run one of the many examples to ensure that our oxDNA2 simulations will run correctly. Starting in the source directory, here is what we can do:

```
$ cd ../examples/USER/cgdna/duplex3
$ ../../../../src/lmp_mpi in.duplex3
```

## Run

To run LAMMPS we have the following shell command:

```
$ /path/to/lammps/src/lmp_mpi [command line options]
```

A typical file structure would look like this:

```
|-  workspace/
|--|-   lammps.in.main
|--|-   lammps.structure.conf
```
There are 4 ways of calling the oxDNA command, depending on a user's preference:

1. Just call the complete path e.g. 
    ```
    $ /path/to/lammps/src/lmp_mpi
    ```

2. Create an `alias` or `export` an environment variable to use as a shortcut e.g.
    ```
    alias lammps=/path/to/lammps/src/lmp_mpi
    $ lammps

    export lammps=/path/to/lammps/src/lmp_mpi
    $ $lammps 
    ```

3. Append the path to the oxDNA binary in your `$PATH` environment variable e.g.
    ```
    export PATH=$PATH:/path/to/lammps/src/
    ```

4. Create a linker using `ln` and create the link in a folder that already exists in your path (may require root priveleges) e.g.
    ```
    sudo ln /path/to/lammps/src/lmp_mpi /usr/local/bin/lammps
    ```

## View

Trajectories can be easily viewed using Ovito, which is our preferred visualisation platform.
