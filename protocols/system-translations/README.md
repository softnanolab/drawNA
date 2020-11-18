# System Translations

## Summary

Demonstration of how to create multiple replicas of oxDNA systems that are translated a certain distance apart.

## Usage

```
python main.py [-n NUMBER] [-o OXDNA] [-i INPUT_FILE]

    NUMBER: Number of replicas
    OXDNA: Path to oxDNA binary
    INPUT_FILE: Path to oxDNA input file
```



## Description

Runs a selection of simulations to gradually create replicas of the same system inside its own box. 
This folder can be used as a template for more complicated setups by changing the following functions:

- `calculate_box_size` - there are many ways of deciding the box size, it can be customised here
- `generate_system` - the main manual generation function, create a system in any starting configuration
- `run_simulation` - if not using oxDNA, then use this function to do file conversions and simulation running
