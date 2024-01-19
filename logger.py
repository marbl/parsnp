import logging
import io
############################################# Logging ##############################################
BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
#These are the sequences need to get colored ouput
CSI=""#"\x1B["
BOLDME = ""#CSI+'\033[1m'
STATUS_BLUE = ""#CSI+'\033[94m'
OK_GREEN = ""#CSI+'\033[92m'#'32m'
SKIP_GRAY = ""#CSI+'\033[37m'
WARNING_YELLOW = ""#CSI+'\033[93m'
ERROR_RED = ""#CSI+'\033[91m'
ENDC = ""#CSI+'0m'
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"

MIN_TQDM_INTERVAL=30


# Logging redirect copied from https://stackoverflow.com/questions/14897756/python-progress-bar-through-logging-module
class TqdmToLogger(io.StringIO):
    """
        Output stream for TQDM which will output to logger module instead of
        the StdOut.
    """
    logger = None
    level = None
    buf = ''
    def __init__(self,logger,level=None):
        super(TqdmToLogger, self).__init__()
        self.logger = logger
        self.level = level or logging.INFO
    def write(self,buf):
        self.buf = buf.strip('\r\n\t ')
    def flush(self):
        self.logger.log(self.level, self.buf)


def formatter_message(message, use_color = True):
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message

COLORS = {
    'DEBUG': BLUE,
    'INFO': WHITE,
    'WARNING': YELLOW,
    'ERROR': RED,
    'CRITICAL': RED
}

class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, datefmt = None, use_color = True):
        logging.Formatter.__init__(self, fmt=msg, datefmt=datefmt)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)

# Instantiate logger
logger = logging.getLogger("Parsnp")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter and add it to the handlers
formatter = ColoredFormatter('%(asctime)s - %(levelname)s - %(message)s',
                             datefmt="%H:%M:%S")
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(ch)
####################################################################################################
