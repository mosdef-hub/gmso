import logging

logging.captureWarnings(True)

formatter = logging.Formatter('[%(asctime)s::%(filename)s::%(funcName)s::%(levelname)s]%(message)s',
        datefmt='%a %d %b %Y %H:%M:%S')

# This logger listens for all debug-level messages
logger = logging.getLogger("TopLog")
logger.setLevel(logging.DEBUG)

# This logger listesn for messages from py.warnings
warn_logger = logging.getLogger('py.warnings')


# This file will print warning messages
warn_handler = logging.FileHandler('TopWarning.out')
warn_handler.setLevel(logging.WARNING)
warn_handler.setFormatter(formatter)

# This file will print debug messages
debug_handler = logging.FileHandler('TopDebug.out')
debug_handler.setLevel(logging.DEBUG)
debug_handler.setFormatter(formatter)

#warn_handler file receives input from both loggers
warn_logger.addHandler(warn_handler)
logger.addHandler(warn_handler)

#debug_handler file receives input from both loggers
logger.addHandler(debug_handler)
warn_logger.addHandler(debug_handler)
