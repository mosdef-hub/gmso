import os
import logging
import logging.config
import json


def start_logging():
    """ Initialize the logger ("TopLog")

    Notes
    ---
    Logger (and py.warnings) messages will be printed to TopWarning.out and
    TopDebug.out in the user's local/exec directory
    Configuration is loaded from 'log_config.json'
    Any modifications to the configuration can be done by modifying the json file
    Or by modifying the logger in this function
    Two loggers are used to catch py.warnings messages and logging messages
    """
    logging.captureWarnings(True)
    path_to_log_config = os.path.join(os.path.dirname(__file__), 'log_config.json')
    logging.config.dictConfig(json.load(open(path_to_log_config, 'r')))

def end_logging():
    """ End logging

    Notes
    ----
    This is done by clearing the logger's associated handlers 
    
    Comments
    --------
    This seems to be cleaner/safer than just setting `logger.handlers = []`
    Creating another list of the handlers is important to avoid simultaneously
    looping and modifying the same list
    We need to reset captureWarnings in order to allow py.warnings to instead
    print to stream rather than the associated filehandler"""
    logger = logging.getLogger("TopLog")
    logging.captureWarnings(False)
    all_handlers = [h for h in logger.handlers]
    for h in all_handlers:
        logger.removeHandler(h)
