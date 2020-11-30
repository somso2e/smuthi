# -*- coding: utf-8 -*-
"""This module is needed for the installation of the package."""

import os
from setuptools import setup
from setuptools.command.install import install
from setuptools.command.develop import develop
import pkg_resources
import os
import subprocess
import sys
import warnings
import shutil
import glob

version = {}
with open("smuthi/version.py") as fp:
    exec(fp.read(), version)
__version__ = version['__version__']


class CustomDevelopCommand(develop):
    def run(self):
        install_pywigxjpf()
        compile_nfmds()
        develop.run(self)


class CustomInstallCommand(install):
    def run(self):
        install_pywigxjpf()
        compile_nfmds()
        install.run(self)


def install_pywigxjpf():
    """If Windows: try to install pywigxjpf"""
    if sys.platform.startswith('win'):
        sys.stdout.write('Try to install pywigxjpf ... \n')
        sys.stdout.flush()
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", "pywigxjpf"], stdout=None, stderr=None)
        except Exception as e:
            warnings.warn('\n*****************************************************\n'
                          'pywigxjpf installation failed.\n'
                          'If you want to benefit from fast Wigner3j calculations:\n'
                          'Try to install manually by "pip install pywigxjpf".'
                          '\n*****************************************************\n',
                          UserWarning)


def compile_nfmds():
    """compile nfmds Fortran extension (unless built on readthedocs)"""
    currdir = os.getcwd()
    if not os.environ.get('READTHEDOCS'):
        nfmds_sources_dirname = pkg_resources.resource_filename('smuthi.linearsystem.tmatrix.nfmds', 'NFM-DS')
        sys.stdout.write('\nCompiling sources at ' + nfmds_sources_dirname + ' ...')
        sys.stdout.flush()
        os.chdir(nfmds_sources_dirname + '/TMATSOURCES')

        sys.stdout.flush()
        try:
            if sys.platform.startswith('win'):
                import encoding_converter
                os.chdir('win')
                subprocess.check_call(
                    ['f2py', '-c', '--compiler=mingw32', '--fcompiler=gnu95', 'TAXSYM_SMUTHI.f90', '-m', 'nfmds'])
                for file in glob.glob(r'nfmds*'):
                    shutil.move(file, '../../../' + file)
            else:
                subprocess.check_call(['f2py', '-c', 'TAXSYM_SMUTHI.f90', '-m', 'nfmds'])
                for file in glob.glob(r'nfmds*'):
                    shutil.move(file, '../../' + file)
            sys.stdout.write(' done.\n')
            sys.stdout.flush()
        except Exception as e:
            raise NameError('Compiling failed.')
    os.chdir(currdir)


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def get_requirements():
    """Return a list of requirements. If windows, don't include pywigxjpf, because we want to have it as an optional
    dependency in that case (will be taken care of in custom install command)"""
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
                    'pycparser']
    if sys.platform.startswith('win'):
        # skip pywigxjpf (to have it optional, if the user has no C compiler)
        sys.stdout.write('Compiling from Windows machine. Skipping pywigxjpf for the moment. '
                         'I will try to install it in post processing.\n')
        sys.stdout.flush()
    else:
        requirements.append('pywigxjpf')
    return requirements

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
    cmdclass={'install': CustomInstallCommand,
              'develop': CustomDevelopCommand},
    package_data={'smuthi.linearsystem.tmatrix.nfmds': ['NFM-DS/*.txt', 'NFM-DS/TMATSOURCES/*.f90', 'NFM-DS/TMATFILES/*',
                                                        'NFM-DS/INPUTFILES/*.dat', 'NFM-DS/OUTPUTFILES/*','nfmds*'],
                  'smuthi': ['_data/*']},
    include_package_data=True,                  
    install_requires=get_requirements(),
    setup_requires=['numpy'],
    extras_require={'cuda':  ['PyCuda']},
    entry_points={'console_scripts': ['smuthi = smuthi.__main__:main']},
    license='MIT',
)
