import os
import sys
import subprocess
from setuptools import setup
from setuptools.command.build_ext import build_ext

version = {}
with open("smuthi/version.py") as fp:
    exec(fp.read(), version)
__version__ = version['__version__']


class CMakeBuild(build_ext):
    def run(self):
        """Run CMake and Makefile"""
        build_dir = os.path.abspath(self.build_temp)
        src_dir = os.path.abspath("smuthi/linearsystem/tmatrix/nfmds/NFM-DS/TMATSOURCES")

        os.makedirs(build_dir, exist_ok=True)
        # Run Makefile inside TMATSOURCES
        subprocess.check_call(["make", "-C", src_dir])


setup(
    name="SMUTHI",
    version=__version__,
    author="Amos Egel",
    author_email="amos.egel@gmail.com",
    url='https://gitlab.com/AmosEgel/smuthi',
    description="Light scattering by multiple particles in thin-film systems",
    packages=['smuthi',
              'smuthi.fields',
              'smuthi.linearsystem',
              'smuthi.linearsystem.tmatrix',
              'smuthi.linearsystem.tmatrix.nfmds',
              'smuthi.linearsystem.particlecoupling',
              'smuthi.postprocessing',
              'smuthi.periodicboundaries',
              'smuthi.utility',
              'smuthi.utility.cython'],
    cmdclass={'build_ext': CMakeBuild},
    include_package_data=True,
    install_requires=[
        'numpy>=1.19.3',
        'argparse',
        'imageio',
        'matplotlib',
        'mpmath',
        'numba',
        'pyyaml',
        'scipy',
        'sympy',
        'tqdm',
        'h5py',
        'pycparser',
        'psutil'
    ],
    extras_require={'cuda': ['PyCuda']},
    license='MIT',
)

