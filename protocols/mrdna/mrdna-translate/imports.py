import subprocess

try:
    subprocess.check_output(['nvidia-smi'])
except FileNotFoundError:
    raise FileNotFoundError(
        'nvidia-smi did not work so CUDA probably '
        'is not installed, so this example will not work!'
    )
    

try:
    import mrdna
except ImportError:
    raise ImportError(
        'mrDNA is not installed, '
        'see README.md for installation guide'
    )

try:
    import softnanotools
except ImportError:
    raise ImportError(
        'softnanotools is not installed, '
        'install using `pip install softnanotools`'
    )