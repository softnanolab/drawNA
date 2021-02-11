import subprocess
from softnanotools.logger import Logger
logger = Logger(__name__)

logger.info('Checking for nvcc...')
try:
    subprocess.check_output(['nvcc'], stderr=subprocess.DEVNULL)
except FileNotFoundError:
    raise FileNotFoundError(
        'nvcc did not work so CUDA probably '
        'is not installed, so this example will not work!'
    )
except subprocess.CalledProcessError:
    logger.info('Success!')
    
logger.info('Trying to import mrdna...')
try:
    import mrdna
    logger.info('Success!')
except ImportError:
    raise ImportError(
        'mrDNA is not installed, '
        'see README.md for installation guide'
    )

logger.info('Trying to run oxDNA...')
try:
    subprocess.check_output(['oxDNA'], stderr=subprocess.DEVNULL)
except FileNotFoundError:
    raise FileNotFoundError(
        'oxDNA was not found at oxDNA - please add '
        'the oxDNA executable to your PATH'
    )
except subprocess.CalledProcessError:
    logger.info('Success!')

logger.info('Trying to import softnanotools...')
try:
    import softnanotools
    logger.info('Success!')
except ImportError:
    raise ImportError(
        'softnanotools is not installed, '
        'install using `pip install softnanotools`'
    )

#logger.info('Trying to import vmd...')
#try:
#    import vmd
#    logger.info('Success!')
#except ImportError:
#    raise ImportError(
#        'vmd is not installed, '
#    )