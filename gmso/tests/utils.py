
import os


def get_path(filename):
    """Given a test filename return its path"""
    _path = os.path.join(os.path.split(__file__)[0], 'files', filename)
    return _path
