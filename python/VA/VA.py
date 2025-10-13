# ppp.py — simple ctypes wrapper for an EXISTING libppp_ball.so.
# This does NOT recompile or modify the shared library.

import ctypes
import numpy as np

# Load the shared library (edit path if needed)
_LIB = ctypes.CDLL("libppp_ball.so")

# Bind function signatures ONCE
_LIB.ppp_generate.argtypes = [
    ctypes.c_double,             # lambda
    ctypes.c_ulonglong,          # M
    ctypes.c_int,                # d
    ctypes.c_ulonglong,          # seed
    ctypes.POINTER(ctypes.c_double),  # out points (M*d)
    ctypes.POINTER(ctypes.c_double)   # out R
]
_LIB.ppp_generate.restype = ctypes.c_int

# Optional convenience wrapper
def ppp_generate(lambda_val, M, d, seed):
    """
    Calls ppp_generate from libppp_ball.so and returns (points ndarray 
[M,d], R).
    """
    M = int(M); d = int(d)
    pts = np.empty((M, d), dtype=np.float64)
    R = ctypes.c_double(0.0)

    rc = _LIB.ppp_generate(
        float(lambda_val), M, d, seed,
        pts.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(R)
    )
    if rc != 0:
        raise RuntimeError("Library call failed")

    return pts, float(R.value)


# Load the shared library (edit path if needed)
_LIB2 = ctypes.CDLL("libinitial_proj.so")

# Bind function signatures ONCE
_LIB2.proj_select_from_points.argtypes = [
    ctypes.POINTER(ctypes.c_double),  # Points (N*d)
    ctypes.c_ulonglong,               # N
    ctypes.c_int,                     # d
    ctypes.POINTER(ctypes.c_int)      # out_indices (size d+1)
]
_LIB2.proj_select_from_points.restype = ctypes.c_int

# Optional convenience wrapper
def proj_select_from_points(Points):
    """
    Calls proj_select_from_points from libinitial_proj.so and returns
    the indices (array of length d+1) of the selected points.
    """
    Points = np.asarray(Points, dtype=np.float64)
    if Points.ndim != 2:
        raise ValueError("Points must be a 2D array of shape (N, d).")
    N, d = Points.shape
    out_idx = np.empty(d + 1, dtype=np.int32)

    rc = _LIB2.proj_select_from_points(
        Points.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_ulonglong(N),
        ctypes.c_int(d),
        out_idx.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    )
    if rc != 0:
        raise RuntimeError(f"Library call failed (rc={rc}).")

    return out_idx

_VA = ctypes.CDLL("libVAedge.so")

_VA.VA_RunLoop.argtypes = [
    ctypes.POINTER(ctypes.c_double),  # Points (N*d)
    ctypes.c_ulonglong,               # N
    ctypes.c_int,                     # d
    ctypes.POINTER(ctypes.c_int),     # Aa_init (d+1)
    ctypes.c_char_p,                  # order
    ctypes.c_int,                     # max_iter
    ctypes.c_char_p,                  # log_path or None
    ctypes.POINTER(ctypes.c_int),     # out_Aa_final (d+1)
    ctypes.POINTER(ctypes.c_double),  # out_p_last (d)
    ctypes.POINTER(ctypes.c_double),  # out_r_last (1)
    ctypes.POINTER(ctypes.c_int),     # out_iters (1)
    ctypes.POINTER(ctypes.c_double),  # out_path_length (1)
    ctypes.POINTER(ctypes.c_int),     # out_K_sum (1)
    ctypes.POINTER(ctypes.c_int),     # out_k_last (1)
    ctypes.POINTER(ctypes.c_int),     # out_hit_inside (1)
]
_VA.VA_RunLoop.restype = ctypes.c_int

def VA_RunLoop(Points, Aa_init, order="max", max_iter=100000, log_path=None):
    Points = np.asarray(Points, dtype=np.float64)
    N, d = Points.shape
    Aa_init = np.asarray(Aa_init, dtype=np.int32)
    assert Aa_init.size == d+1
    out_Aa = np.empty_like(Aa_init)
    p_last = np.empty(d, dtype=np.float64)
    r_last = ctypes.c_double(0.0)
    iters  = ctypes.c_int(0)
    L      = ctypes.c_double(0.0)
    Ksum   = ctypes.c_int(0)
    k_last = ctypes.c_int(0)
    hit    = ctypes.c_int(0)

    rc = _VA.VA_RunLoop(
        Points.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_ulonglong(N), ctypes.c_int(d),
        Aa_init.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        order.encode("ascii"),
        ctypes.c_int(max_iter),
        (log_path.encode("utf-8") if log_path else None),
        out_Aa.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        p_last.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(r_last),
        ctypes.byref(iters),
        ctypes.byref(L),
        ctypes.byref(Ksum),
        ctypes.byref(k_last),
        ctypes.byref(hit),
    )
    if rc != 0:
        raise RuntimeError(f"VA_RunLoop failed rc={rc}")
    return {
        "Aa_final": out_Aa,
        "p_last": p_last,
        "r_last": float(r_last.value),
        "iters": int(iters.value),
        "path_length": float(L.value),
        "K_sum": int(Ksum.value),
        "k_last": int(k_last.value),
        "hit_inside": bool(hit.value),
    }
