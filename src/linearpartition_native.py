"""
Thin ctypes-based bindings to the native LinearPartition shared library.
"""

from __future__ import annotations

import ctypes
import os
import sys
from pathlib import Path
from typing import Optional, Tuple

import numpy as np


def _default_library_path() -> Path:
    root = Path(__file__).resolve().parent
    ext = "dylib" if sys.platform == "darwin" else "so"
    return root / f"liblinearpartition_v.{ext}"


class _LinearPartitionLib:
    """
    Internal helper that loads the shared library once and exposes ctypes
    signatures.
    """

    _instance: Optional["_LinearPartitionLib"] = None

    def __init__(self, library_path: Optional[Path] = None) -> None:
        lib_path = library_path or Path(
            os.environ.get("LINEARPARTITION_NATIVE_LIB", _default_library_path())
        )
        if not lib_path.exists():
            raise OSError(
                f"LinearPartition native library not found at {lib_path}. "
                "Run `make liblinearpartition` to build it or set "
                "LINEARPARTITION_NATIVE_LIB to the compiled library path."
            )

        self._lib = ctypes.CDLL(str(lib_path))
        self._lib.lp_compute_ensemble_energy.argtypes = [
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
        ]
        self._lib.lp_compute_ensemble_energy.restype = ctypes.c_int
        self._lib.lp_compute_bpp_matrix.argtypes = [
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_longlong,
        ]
        self._lib.lp_compute_bpp_matrix.restype = ctypes.c_int
        self._lib.lp_calculate_mea_structure.argtypes = [
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_double,
            ctypes.c_char_p,
            ctypes.c_longlong,
            ctypes.POINTER(ctypes.c_double),
        ]
        self._lib.lp_calculate_mea_structure.restype = ctypes.c_int
        self._lib.lp_get_last_error.argtypes = []
        self._lib.lp_get_last_error.restype = ctypes.c_char_p

    @classmethod
    def get(cls) -> "_LinearPartitionLib":
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    @property
    def lib(self) -> ctypes.CDLL:
        return self._lib


class LinearPartitionNativeClient:
    """
    Minimal client that exposes ensemble free energy calculation through the
    in-process LinearPartition API.
    """

    def __init__(
        self,
        beamsize: int = 100,
        no_sharp_turn: bool = True,
        verbose: bool = False,
        dangles: int = 2,
        library_path: Optional[Path] = None,
    ) -> None:
        self.beamsize = beamsize
        self.no_sharp_turn = no_sharp_turn
        self.verbose = verbose
        self.dangles = dangles
        self._lib = _LinearPartitionLib(library_path).lib

    def _call_native(self, sequence: str) -> float:
        seq = self._prepare_sequence(sequence)
        c_sequence = seq.encode("ascii")
        result = ctypes.c_double()
        status = self._lib.lp_compute_ensemble_energy(
            c_sequence,
            ctypes.c_int(self.beamsize),
            ctypes.c_int(1 if self.no_sharp_turn else 0),
            ctypes.c_int(1 if self.verbose else 0),
            ctypes.c_int(self.dangles),
            ctypes.byref(result),
        )
        if status != 0:
            error = self._lib.lp_get_last_error()
            message = error.decode("utf-8") if error else "Unknown error"
            raise RuntimeError(f"LinearPartition native call failed: {message}")
        return float(result.value)

    def calculate_partition_function(self, sequence: str) -> float:
        return self._call_native(sequence)

    def _prepare_sequence(self, sequence: str) -> str:
        seq = sequence.strip().upper().replace("T", "U")
        if not seq:
            raise ValueError("Sequence must be non-empty")
        return seq

    def calculate_bpp_matrix(self, sequence: str, cutoff: float = 0.0) -> np.ndarray:
        seq = self._prepare_sequence(sequence)
        n = len(seq)
        buffer = (ctypes.c_double * (n * n))()
        status = self._lib.lp_compute_bpp_matrix(
            seq.encode("ascii"),
            ctypes.c_int(self.beamsize),
            ctypes.c_int(1 if self.no_sharp_turn else 0),
            ctypes.c_int(1 if self.verbose else 0),
            ctypes.c_int(self.dangles),
            ctypes.c_double(cutoff),
            buffer,
            ctypes.c_longlong(n * n),
        )
        if status != 0:
            error = self._lib.lp_get_last_error()
            message = error.decode("utf-8") if error else "Unknown error"
            raise RuntimeError(f"LinearPartition native call failed: {message}")
        matrix = np.ctypeslib.as_array(buffer, shape=(n * n,))
        return matrix.reshape((n, n)).copy()

    def calculate_mea_structure(self, sequence: str, gamma: float = 3.0) -> Tuple[str, float]:
        seq = self._prepare_sequence(sequence)
        buf = ctypes.create_string_buffer(len(seq) + 1)
        energy = ctypes.c_double()
        status = self._lib.lp_calculate_mea_structure(
            seq.encode("ascii"),
            ctypes.c_int(self.beamsize),
            ctypes.c_int(1 if self.no_sharp_turn else 0),
            ctypes.c_int(1 if self.verbose else 0),
            ctypes.c_int(self.dangles),
            ctypes.c_double(gamma),
            buf,
            ctypes.c_longlong(len(seq) + 1),
            ctypes.byref(energy),
        )
        if status != 0:
            error = self._lib.lp_get_last_error()
            message = error.decode("utf-8") if error else "Unknown error"
            raise RuntimeError(f"LinearPartition native call failed: {message}")
        return buf.value.decode("ascii"), float(energy.value)

