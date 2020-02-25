from __future__ import division, print_function, absolute_import

import os
import sys
from os.path import join, dirname
from distutils.sysconfig import get_python_inc
import subprocess
import numpy
from numpy.distutils.misc_util import get_numpy_include_dirs

try:
    from numpy.distutils.misc_util import get_info
except ImportError:
    raise ValueError("numpy >= 1.4 is required (detected %s from %s)" %
                     (numpy.__version__, numpy.__file__))


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from scipy._build_utils.system_info import get_info as get_system_info

    config = Configuration('cy_special', parent_package, top_path)

    # adding the cython port files
    config.add_extension('_cy_special', sources=['_cy_special.c'], depends=['*.pyx', '*.pxd', '*.pxi'])
    config.add_extension('round', sources=['round.c'])

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path=' ').todict())