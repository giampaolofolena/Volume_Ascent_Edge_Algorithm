import os, sys, ctypes, numpy as np

# locate built libs in ../../build or alongside this file if installed
def _find_lib(name):
    here = os.path.abspath(os.path.dirname(__file__))
    candidates = []
    # local dev build dir
    candidates += [os.path.abspath(os.path.join(here, "..", "..", "build", f"lib{name}.so")),
                   os.path.abspath(os.path.join(here, "..", "..", "build", f"{name}.dll")),
                   os.path.abspath(os.path.join(here, "..", "..", "build", f"lib{name}.dylib"))]
    # installed next to package
    candidates += [os.path.join(here, f"lib{name}.so"),
                   os.path.join(here, f"{name}.dll"),
                   os.path.join(here, f"lib{name}.dylib")]
    for p in candidates:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"Cannot find {name} shared library. Looked in: {candidates}")

# Load libs
_LIB_PPP = ctypes.CDLL(_find_lib("ppp_ball"))
_LIB_IP  = ctypes.CDLL(_find_lib("initial_proj"))
_LIB_VA  = ctypes.CDLL(_find_lib("VAedge"))

# Signatures
_LIB_PPP.ppp_generate.argtypes = [
    ctypes.c_double, ctypes.c_ulonglong, ctypes.c_int, ctypes.c_ulonglong,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)
]
_LIB_PPP.ppp_generate.restype = ctypes.c_int

_LIB_IP.proj_select_from_points.argtypes = [
    ctypes.POINTER(ctypes.c_double), ctypes.c_ulonglong, ctypes.c_int,
    ctypes.POINTER(ctypes.c_int)
]
_LIB_IP.proj_select_from_points.restype = ctypes.c_int

_LIB_VA.VA_RunLoop.argtypes = [
    ctypes.POINTER(ctypes.c_double), ctypes.c_ulonglong, ctypes.c_int,
    ctypes.POINTER(ctypes.c_int), ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)
]
_LIB_VA.VA_RunLoop.restype = ctypes.c_int

# Optional: expose CircumPart for tests (rows-aware)
try:
    _LIB_VA.VA_CircumPart.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)
    ]
    _LIB_VA.VA_CircumPart.restype = ctypes.c_int
    HAS_CIRCUMPART = True
except Exception:
    HAS_CIRCUMPART = False

# Python convenience wrappers
def ppp_generate(lambda_val, M, d, seed):
    pts = np.empty((int(M), int(d)), np.float64)
    R   = ctypes.c_double(0.0)
    rc = _LIB_PPP.ppp_generate(
        float(lambda_val), int(M), int(d), int(seed),
        pts.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(R)
    )
    if rc != 0: raise RuntimeError(f"ppp_generate failed rc={rc}")
    return pts, float(R.value)

def proj_select_from_points(points: np.ndarray):
    points = np.asarray(points, dtype=np.float64)
    N, d = points.shape
    idx = np.empty(d+1, np.int32)
    rc = _LIB_IP.proj_select_from_points(
        points.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_ulonglong(N), ctypes.c_int(d),
        idx.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    )
    if rc != 0: raise RuntimeError(f"proj_select_from_points failed rc={rc}")
    return idx

def VA_RunLoop(points: np.ndarray, Aa_init: np.ndarray, order="max", max_iter=100000, log_path=None):
    points = np.asarray(points, dtype=np.float64)
    N, d = points.shape
    Aa_init = np.asarray(Aa_init, dtype=np.int32)
    if Aa_init.shape != (d+1,): raise ValueError("Aa_init must be shape (d+1,)")

    Aa_final = np.empty(d+1, np.int32)
    p_last   = np.empty(d, np.float64)
    r_last   = ctypes.c_double(0.0)
    iters    = ctypes.c_int(0)
    L_path   = ctypes.c_double(0.0)
    K_sum    = ctypes.c_int(0)
    k_last   = ctypes.c_int(0)
    hit      = ctypes.c_int(0)

    rc = _LIB_VA.VA_RunLoop(
        points.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_ulonglong(N), ctypes.c_int(d),
        Aa_init.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        order.encode("ascii"),
        ctypes.c_int(max_iter),
        (log_path.encode("utf-8") if log_path else None),
        Aa_final.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        p_last.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(r_last),
        ctypes.byref(iters),
        ctypes.byref(L_path),
        ctypes.byref(K_sum),
        ctypes.byref(k_last),
        ctypes.byref(hit)
    )
    if rc != 0: raise RuntimeError(f"VA_RunLoop failed rc={rc}")
    return dict(Aa_final=Aa_final, p_last=p_last, r_last=float(r_last.value),
                iters=int(iters.value), path_length=float(L_path.value),
                K_sum=int(K_sum.value), k_last=int(k_last.value),
                hit_inside=bool(hit.value))
