# Contributing

## Initial Notes

This package is currently not deployed and installation requires only `pip` and other Python packages. That being said, most simulation tools are more complicated to build especially on other operating systems so this should be considered when making additions.

- In the future there may be some C++ libraries and bindings but this does not exist yet.

- Only Python 3 is supported and there are no plans for using Python 2.

- Function annotations are used extensively but are not checked - these are often used instead of docstrings for simple functions:
```python
def function(arg: type) -> return-type:
    return something

def complex(arg_1: int, arg_2: float) -> pd.DataFrame:
    """
    Description of complex function 
    because it is complex
    """
    return something_complex
```
- No linters are currently used - make an issue and correct this if you think it is necessary

- We are looking for someone to configure an automatic documentation generator so please raise an issue if you can do this

## Directory Structure

### `img`

- Contains media used for examples and documentation
- Likely to be moved in the future

### `drawNA`

- Main Python package
- Contains all Python code that will be installed
- Structure:
  - `oxdna` - all DNA related code
    - `oxdna.nucleotide`
    - `oxdna.strand`
    - `oxdna.system`
    - `oxdna.utils`
  - `lattice` - all lattice and routing code
    - `lattice._lattice`
    - `lattice.edge`
    - `lattice.node`
    - `lattice.route`
    - `lattice.utils`
  - `polygons` - basic code for storing geometry
  - `readers` - simulation file readers
  - `tools` - general computing tools

### `polygons`

- Contains defunct examples of 2D polygons

### `protocols`

- This folder is important and where all successful examples should be kept
- There is a script called `create_protocol.py` which is used like this: 

```python create_protocol.py protocol-name```
- It will properly create a protocols folder with the minimum required files to ensure consistency between protocols
- Please provide detailed descriptions of usage so that others can learn from this folder

### `sandbox`

- This folder is also important but less formal than the `protocols` folder
- It contains random experiments by any contributors and is less maintained
- Start here when cloning to experiment with the code and see what others have done
- If copying and playing with another entry than make a new file where possible
- Code in this folder is not guaranteed to work and may be removed in the future
- Code in this folder is not kept up-to-date and requests to update this code will probably not be granted

### `software`

- Here are general instructions on using some of the external software that this repository is used for

### `test`

- `pytest` code is in here
- All code should be tested in here using standard `pytest` conventions

### `.gitattribute`

- Remove Jupyter Notebooks from the language line count


### `.gitignore`

- Standard usage
- Note that all generated files should be on here except those for testing
- Make sure no `.DS_Store` files are added ever!

### `.travis.yml`

- Currently CI is done using this file alone
- CI is only done on `xenial`
- CI is likely to migrate to GitHub workflows so keep an eye out for `.github/workflows` in the main directory

### `CONTRIBUTING.md`

- This file

### `LICENSE`

- Auto-generated when repository was created

### `README.md`

- Main README and long description for `pip`

### `setup.py`

- Usage:
```
pip install .
```
## Adding new functionality

Say we wanted to add a new set of functions to the repository. For example, a new data format is released and we want the `oxdna.system.System` class to output to this file format. Here are the steps that should be followed:

1. Create new branch
```
git checkout -b new-feature
```
2. Add our changes and commit
```
git commit -a -m 'new feature finished'
```
3. Add tests, run `pytest` and commit when finished
```
touch test/test_new_feature.py
...
# Add testing etc.
...
pytest
git commit -a -m 'added testing and confirmed pytest works'
```
4. Make new upstream branch and push
```
git push --set-upstream origin new-feature
```
5. Make a pull request
6. (Administrator) Do review and merge when no errors are raised

## Tracking issues

Issues are tracked in the standard way, by using GitHub issues. There is no system for doing this at the time of writing.
