__version__ = "1.0.0"


import os

_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data_path():
    return os.path.join(_ROOT, 'data')