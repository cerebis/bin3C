import numpy as np
import os
import sys


def make_random_seed():
    """
    Provide a random seed value between 1 and 10 million.
    :return: integer random seed
    """
    return np.random.randint(1000000, 10000000)


def make_dir(path, exist_ok=False):
    """
    Convenience method for making directories with a standard logic.
    An exception is raised when the specified path exists and is not a directory.
    :param path: target path to create
    :param exist_ok: if true, an existing directory is ok. Existing files will still cause an exception
    """
    if not os.path.exists(path):
        os.mkdir(path)
    elif not exist_ok:
        raise IOError('output directory already exists!')
    elif os.path.isfile(path):
        raise IOError('output path already exists and is a file!')


def app_path(subdir, filename):
    """
    Return path to named executable in a subdirectory of the running application

    :param subdir: subdirectory of application path
    :param filename: name of file
    :return: absolute path
    """
    return os.path.join(sys.path[0], subdir, filename)


