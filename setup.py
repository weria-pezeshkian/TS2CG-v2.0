import os, sys
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from pathlib import Path
import multiprocessing


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])

class CMakeBuild(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)

    def build_cmake(self, ext):
        sourcedir = Path().absolute()
        extdir = Path(self.get_ext_fullpath(ext.name)).parent.absolute()

        cmake_args = [
            f'-DCMAKE_INSTALL_PREFIX={extdir}',
            f'-DPYTHON_EXECUTABLE={sys.executable}',
            '-DCMAKE_BUILD_TYPE=Release'
        ]

        build_temp = Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)

        num_cores = multiprocessing.cpu_count()
        build_args = ['--config', 'Release', '--', f'-j{num_cores}']

        os.chdir(build_temp)
        print(sourcedir, cmake_args)

        self.spawn(['cmake', sourcedir] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.'] + build_args)
            self.spawn(['cmake', '--install', '.'])
        os.chdir(sourcedir)

setup(
    name='TS2CG',
    version='2.0',
    packages=find_packages(include=['TS2CG']),
    ext_modules=[CMakeExtension('TS2CG')],
    cmdclass={'build_ext': CMakeBuild},
    package_data={
        'TS2CG': ['SOL', 'PLM', 'PCG', 'CMakeLists.txt',
                  'core/*', 'tools/*', 'cpp/*'],
    },
    include_package_data=True,
    python_requires='>=3.6',
    install_requires=[],
    entry_points={
        'console_scripts': [
            'TS2CG=TS2CG.run_modules:main',
        ],
    },
)
