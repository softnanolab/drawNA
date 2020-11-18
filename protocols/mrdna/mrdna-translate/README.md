#mrDNA Translate

## Summary

In this directory is an example of how to run mrDNA using Python whilst creating replicas that represent different stages during the simulation of the same DNA nanostructure.

## Prerequisites

- Python packages
    - `numpy`
    - `scipy`
    - `appdirs`
    - `mdanalysis`
    - `cadnano`
    - `matplotlib`
    - `shapely`
    - `pandas`

- `CUDA`
This is a GPU interface for NVIDIA GPUs, if you do not have a CUDA GPU, then you cannot use mrDNA and you  cannot use this example.


- `cadnano2.5`
The classic DNA design program can be installed like this:
```
git clone https://github.com/cadnano/cadnano2.5
cd cadnano2.5
python setup.py install
```

- `arbd`
This program cannot be downloaded from the command line but instead must be downloaded using a web browser from [here](http://bionano.physics.illinois.edu/arbd). Once downloaded do the following:
```sh
pushd /tmp
cp /path/to/arbd .
pushd arbd/src
make -j8
export PATH=$PATH:/tmp/arbd/src 
popd
popd
```

- `mrdna`
```sh
pushd /tmp
git clone https://gitlab.engr.illinois.edu/tbgl/tools/mrdna
pushd mrdna
python setup.py install
popd
```

- `softnanotools`
```sh
pip install softnanotools
```

- `oxDNA`
```sh
pushd /tmp
git clone https://github.com/lorenzo-rovigatti/oxDNA
pushd oxDNA
mkdir build
pushd build
cmake ..
make -j8
export PATH=$PATH:/tmp/oxDNA/build/bin 
popd 
popd 
popd
```

## Usage

```python main.py```

## Description

Create an extended description of the program here.
