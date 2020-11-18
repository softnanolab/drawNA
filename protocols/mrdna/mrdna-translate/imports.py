import subprocess
from softnanotools.logger import Logger
logger = Logger(__name__)

logger.info('Checking for nvidia-smi...')
try:
    subprocess.check_output(['nvidia-smi'])
    logger.info('Success!')
except FileNotFoundError:
    raise FileNotFoundError(
        'nvidia-smi did not work so CUDA probably '
        'is not installed, so this example will not work!'
    )
    
logger.info('Trying to import mrdna...')
try:
    import mrdna
    logger.info('Success!')
except ImportError:
    raise ImportError(
        'mrDNA is not installed, '
        'see README.md for installation guide'
    )

logger.info('Trying to import softnanotools...')
try:
    import softnanotools
    logger.info('Success!')
except ImportError:
    raise ImportError(
        'softnanotools is not installed, '
        'install using `pip install softnanotools`'
    )