from pathlib import Path
from typing import Union, List
from .base import CPPModule

class PCG(CPPModule):
    """Wrapper for PCG (Point Cloud Generator) C++ module"""
    def __init__(self):
        super().__init__("PCG")

    def __call__(self, *args: str):
        """Run PCG with given arguments"""
        return self.run(list(args))

class PLM(CPPModule):
    """Wrapper for PLM (Pointillism) C++ module"""
    def __init__(self):
        super().__init__("PLM")

    def __call__(self, *args: str):
        """Run PLM with given arguments"""
        return self.run(list(args))

class SOL(CPPModule):
    """Wrapper for SOL (Solvate) C++ module"""
    def __init__(self):
        super().__init__("SOL")

    def __call__(self, *args: str):
        """Run SOL with given arguments"""
        return self.run(list(args))
