# oxDNA installation

## Summary

oxDNA is both a potential used in molecular simulations and a software that can be used for molecular dynamics and Monte Carlo simulations.
In this set of instructions, we will describe how to download and install oxDNA, and run it.

## Install

To install, we clone the repository of the softnanolab group, where a live version-controlled update of the latest oxDNA version found from soruceforge is located.

```bash
$ git clone https://github.com/debeshmandal/oxDNA.git
$ cd oxDNA
$ mkdir build
$ cd build
$ cmake ..
$ make -j8
```

## Run

To run oxDNA we have the following shell command:

```
$ /path/to/oxDNA/build/bin/oxDNA input_script
```

A typical file structure would look like this:

```
|-  workspace/
|--|-   oxdna.in.main
|--|-   oxdna.structure.conf
|--|-   oxdna.structure.top
```
There are 4 ways of calling the oxDNA command, depending on a user's preference:

1. Just call the complete path e.g. 
    ```
    $ /path/to/oxDNA/build/bin/oxDNA
    ```
2. Create an `alias` or `export` an environment variable to use as a shortcut e.g.
    ```
    alias oxDNA=/path/to/oxDNA/buid/bin/oxDNA
    $ oxDNA input_script

    export oxDNA=/path/to/oxDNA/buid/bin/oxDNA
    $ $oxDNA input_script
    ```

3. Append the path to the oxDNA binary in your `$PATH` environment variable e.g.
    ```
    export PATH=$PATH:/path/to/oxDNA/build/bin/
    ```

4. Create a linker using `ln` and create the link in a folder that already exists in your path (may require root priveleges) e.g.
    ```
    sudo ln /path/to/oxDNA/build/bin/oxDNA /usr/local/bin/oxDNA
    ```

## View

Trajectories can be easily viewed using Ovito, which is our preferred visualisation platform.
