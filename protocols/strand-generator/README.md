# Strand Generator

## Summary

In this folder, we manually create DNA strands using the `oxdna.Nucleotide` and `oxdna.Strand` classes.

A function that takes care of this for us is stored as `oxdna.strand.generate_helix`.

## Usage

```python main.py -n <strand_length> [--double]```

## Description

In `main.py` there are a number of functions relating to the generation of new `Strand` instances via creating lists of nucleotides that can be passed to the `Strand` initialisation phase. Note that some of the constant names are incorrectly labelled, but this is noted in the actual implementation. Here are the steps of the main functions:

### `get_5p` (now `Nucleotide.make_5p`)

- Create a rotation matrix to describe a rotation of 35.9&deg; aroung the `a3` axis of the previous `Nucleotide`
- Get a new `a1` vector that has now been rotated appropriately
- Shift the new position of the base in the `a3` direction by a spacing of 0.39
- Shift again by `-(POS_BASE=0.5)` in the new `a1` direction
- Take the cross product of the `a3` and `a1` direction to give the new `a2` direction and shift by -0.105 to give the correct position

### `get_across` (now `Nucleotide.make_across`)

- Reverse the vector directions
- Create a new COM by shifting 2x the sum of `POS_BASE` and `BASE_BASE`
