# -*- coding: utf-8 -*-
"""This module is needed for the installation of the package."""

from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.build_ext import build_ext

import sys
import subprocess

def pip_install(package):
    proc_id = subprocess.check_call([sys.executable, "-m", "pip", "install", package])
pip_install('wheel')
pip_install('numpy')

from wheel.bdist_wheel import bdist_wheel
from numpy.distutils.core import setup
from numpy.distutils.core import Extension

import pkg_resources
import os


version = {}
with open("smuthi/version.py") as fp:
    exec(fp.read(), version)
__version__ = version['__version__']


class PrepareCommand(develop):
    def run(self):
        prepare_nfmds()
        

class CustomDevelopCommand(develop):
    def run(self):
        prepare_nfmds()
        develop.run(self)


class CustomInstallCommand(install):
    def run(self):
        prepare_nfmds()
        install.run(self)


class CustomBdistWheelCommand(bdist_wheel):
    def run(self):
        if sys.platform.startswith('win'):
            # Skip F2Py. Before "python setup.py bdist_wheel", you need to call in advance:
            # python setup.py prepare
            # python setup.py build_ext --inplace --compiler=mingw32 --fcompiler=gnu95 -f
            self.distribution.ext_modules = []
        bdist_wheel.run(self)


def prepare_nfmds():
    """If Windows: Call encoding converter to get Fortran sources with valid encoding"""
    if sys.platform.startswith('win'):
        currdir = os.getcwd()
        nfmds_sources_dirname = pkg_resources.resource_filename('smuthi.linearsystem.tmatrix.nfmds', 'NFM-DS')
        os.chdir(nfmds_sources_dirname + '/TMATSOURCES')
        with open("encoding_converter.py") as fp:
            exec(fp.read(), version)
        os.chdir(currdir)


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def get_requirements():
    """Return a list of requirements, depending on the operating system."""
    requirements = ['argparse',
                    'imageio',
                    'matplotlib',
                    'mpmath',
                    'numpy',
                    'numba',
                    'pyyaml',
                    'scipy',
                    'sympy',
                    'tqdm',
                    'h5py',
                    'pycparser',
                    'psutil']
    if sys.platform.startswith('win'):
        # this package offers windows binaries:
        requirements.append('pywigxjpf-win')
    else:
        requirements.append('pywigxjpf')
    return requirements


def get_extensions(static=True):
    """Depending on the platform, pick suitable Extension object. This is important such that if
    MinGW is used, the DLLs are statically linked (otherwise the user cannot use the binary if he 
    doesn't have MinGW on his computer, too."""
    f2py_options = ['only:', 'tlay','taxsym','tnonaxsym', ':']
    if os.environ.get('READTHEDOCS'):
        return []
    if sys.platform.startswith('win'):
        if static:
            extra_link_args = ["-static", "-static-libgfortran", "-static-libgcc"]            
        else:
            extra_link_args = []
        return [Extension('smuthi.linearsystem.tmatrix.nfmds.nfmds',
                          ['smuthi/linearsystem/tmatrix/nfmds/NFM-DS/TMATSOURCES/win/TAXSYM_SMUTHI.f90'],
                          extra_link_args=extra_link_args,f2py_options=f2py_options)]
    else:
        return [Extension('smuthi.linearsystem.tmatrix.nfmds.nfmds',
                          ['smuthi/linearsystem/tmatrix/nfmds/NFM-DS/TMATSOURCES/TAXSYM_SMUTHI.f90'],f2py_options=f2py_options)]


setup(
    name="SMUTHI",
    version=__version__,
    author="Amos Egel",
    author_email="amos.egel@gmail.com",
    url='https://gitlab.com/AmosEgel/smuthi',
    description="Light scattering by multiple particles in thin-film systems",
    long_description=read('README.rst'),
    packages=['smuthi',
              'smuthi.fields',
              'smuthi.linearsystem',
              'smuthi.linearsystem.tmatrix',
              'smuthi.linearsystem.tmatrix.nfmds',
              'smuthi.linearsystem.particlecoupling',
              'smuthi.postprocessing',
              'smuthi.utility'],
    ext_modules=get_extensions(True),
    cmdclass={'prepare': PrepareCommand,
              'compile': build_ext,
              'install': CustomInstallCommand,
              'develop': CustomDevelopCommand,
              'bdist_wheel': CustomBdistWheelCommand},
    package_data={'smuthi.linearsystem.tmatrix.nfmds': ['NFM-DS/*.txt', 'NFM-DS/TMATSOURCES/*.f90', 'NFM-DS/TMATFILES/*',
                                                        'NFM-DS/INPUTFILES/*.dat', 'NFM-DS/OUTPUTFILES/*','nfmds*'],
                  },
    include_package_data=True,                  
    install_requires=get_requirements(),
    setup_requires=['numpy', 'wheel'],
    extras_require={'cuda':  ['PyCuda']},
    license='MIT',
)
