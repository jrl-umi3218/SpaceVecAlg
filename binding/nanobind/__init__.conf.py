# Copyright 2025-2026 Inria

import os
import platform
import contextlib

class DllDirectoryManager(contextlib.AbstractContextManager):
    """Context manager for managing DLL directories on Windows.
    On Windows with Python 3.8+, Python doesn't search DLL in PATH anymore
    We must specify DLL search path manually with `os.add_dll_directory`
    On Linux and MacOS, we can use DYLD_LIBRARY_PATH and LD_LIBRARY_PATH
    environment variables, or simply the rpath search for shared libraries.
    ref: https://github.com/python/cpython/issues/87339#issuecomment-1093902060
         https://docs.python.org/3/library/os.html#os.add_dll_directory
    """

    def safe_add_dll_directory(self, dll_dir: str):
        """Add a DLL directory to the search path, if path exists."""
        if os.path.isdir(dll_dir):
            self.dll_dirs.append(os.add_dll_directory(dll_dir))

    def __enter__(self):
        self.dll_dirs = []
        return self

    def __exit__(self, *exc_details):
        for d in self.dll_dirs:
            d.close()

if platform.system() == 'Windows':
    with DllDirectoryManager() as mgr:
        module_path = os.path.dirname(__file__)
        mgr.safe_add_dll_directory(os.path.join(module_path, '..', '..', '..', 'bin'))
        dll_dirs = []
        for dll_dir in dll_dirs:
            mgr.safe_add_dll_directory(os.path.join(module_path, dll_dir))

# from .sva_pywrap import *  # noqa
# if hasattr(sva_pywrap, "__version__"):
#     __version__ = sva_pywrap.__version__ # noqa

# Re-export all public names from sva_pywrap
# This is done to get a flat package structure where types appear as
# sva.<classname> instead of sva.sva_pywrap.<symbol>
#
# Consider symbols prefixed with an underscore _<symbol> as private.
# They remain accessible as sva.sva_pywrap._<symbol>
from . import sva_pywrap
for name in dir(sva_pywrap):
    if not name.startswith("_"):
        globals()[name] = getattr(sva_pywrap, name)

__all__ = [name for name in dir(sva_pywrap) if not name.startswith("_")]
