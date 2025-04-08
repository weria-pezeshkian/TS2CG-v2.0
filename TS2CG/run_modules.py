import argparse
import subprocess
from pathlib import Path
from importlib.metadata import version

from TS2CG.tools.domain_placer import DOP
from TS2CG.tools.circular_domains import DAI
from TS2CG.tools.inclusion_updater import INU
from TS2CG.tools.dir_visualizer import VIS
from .cpp import PCG, PLM, SOL


from typing import List, Union, Optional
import logging
import os

logger = logging.getLogger(__name__)

def run_python_module(module_name, args):
    """
    run the specified python module with given arguments.
    """
    if module_name == 'DOP':
        DOP(args)
    elif module_name == 'DAI':
        DAI(args)
    elif module_name == 'INU':
        INU(args)
    elif module_name == 'VIS':
        VIS(args)
    else:
        print(f"Unknown Python module: {module_name}")

def run_cpp_module(module_name, args):
    """Run the specified CPP module with given arguments."""
    cpp_modules = {
        'PCG': PCG(),
        'PLM': PLM(),
        'SOL': SOL()
    }

    #if module_name == "PLM":
    #    try:
    #        from TS2CG.core.mesh_converter import adjust_mesh_file
    #        args=adjust_mesh_file(args)
    #    except ValueError:
    #        pass
 
    try:
        module = cpp_modules[module_name]
        module(*args)
    except KeyError:
        logger.error(f"Unknown CPP module: {module_name}")
        raise ValueError(f"Unknown CPP module: {module_name}")
    except Exception as e:
        logger.error(f"Failed to run {module_name}: {e}")
        raise

def main():
    """
    main entry point for the TS2CG command-line interface.
    """

    # submodules are handled differently based on their type
    # define the c++ submodules
    cpp_modules = ['SOL', 'PLM', 'PCG']

    # define the python based modules
    python_modules = ['DOP', 'DAI', 'INU','VIS']

    # parse arguments before a calling module
    parser = argparse.ArgumentParser(
        description='TS2CG: converts triangulated surfaces (ts) to coarse-grained membrane models',
        prog='TS2CG',
    )

    parser.add_argument(
        'module',
        choices=cpp_modules + python_modules,
        help='choice of which module to run'
    )

    parser.add_argument(
        'args',
        nargs=argparse.REMAINDER,
        help='arguments for the chosen module'
    )

    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f'%(prog)s {version("ts2cg")}'
    )

    args = parser.parse_args()

    # call the right subroutine based on the module type
    if args.module in cpp_modules:
        run_cpp_module(args.module, args.args)
    elif args.module in python_modules:
        run_python_module(args.module, args.args)

if __name__ == '__main__':
    main()
