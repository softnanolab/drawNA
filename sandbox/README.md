# Sandbox

The sandbox exists to experiment with code, in particular standalone scripts

## Usage

When adding a script to the sandbox, try and do the following:

- Add a docstring to the top of the file (using triple quotes)
- if the script will be callable from the command-line (which it probably will be) put an if __name__ == '__main__': statement at the end of the script
- Add a summary of the script, including filename and brief description below

## Scripts

- `empty.py` - contains nothing, acts as a template/placeholder
- `file2vertices.py` - contains an algorithm to extract vertices from different file types and make them cyclical, currently does not work
- `point_in_polygon.py` - experimental 2D point in polygon code
- `old_polygons.py`- old functions from `origamiUROP.polygons`
- `find_crossovers.py` - experimental algorithms which could be used to determine row sizes for DNA lattices which have multiple columns as opposed to one in DNA Snake. Also contains updated find_crossover_locations function.
- `staples_from_route` - new set of `Staple` classes based off of `DNANode` and `DNAEdge` for the generation of staples using the scaffold strand route `LatticeRoute`
