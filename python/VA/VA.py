# VAwrapper.py

import ctypes as ct
import numpy as np
import os
from pathlib import Path

# project root = parent of "python" folder
_PY_FILE = Path(__file__).resolve()
_PROJ_ROOT = _PY_FILE.parents[2]   # .../Volume_Ascent_Edge_Algorithm
_CWD = Path.cwd().resolve()

def _resolve_lib(env_var: str, base: str) -> str:
    """
    Resolve shared library path priority:
      1) ENV: env_var
      2) <PROJ_ROOT>/build/lib{base}.(so|dylib|dll)
      3) <PROJ_ROOT>/lib{base}.(so|dylib|dll)
      4) <PROJ_ROOT>/cpp/lib{base}.(so|dylib|dll)
      5) <CWD>/build/lib{base}.(so|dylib|dll)
      6) <CWD>/lib{base}.(so|dylib|dll)
      7) <CWD>/cpp/lib{base}.(so|dylib|dll)
    Returns absolute path or raises FileNotFoundError.
    """
    p = os.getenv(env_var)
    if p and Path(p).exists():
        return str(Path(p).resolve())

    roots = [
        _PROJ_ROOT / "build",
        _PROJ_ROOT,
        _PROJ_ROOT / "cpp",
        _CWD / "build",
        _CWD,
        _CWD / "cpp",
    ]
    exts = [".so", ".dylib", ".dll"]
    for root in roots:
        for ext in exts:
            cand = (root / f"lib{base}{ext}").resolve()
            if cand.exists():
                return str(cand)
    raise FileNotFoundError(
        f"{env_var} not set and no lib{base} found in search roots: "
        f"{', '.join(str(r) for r in roots)}"
    )

# ------------------------- ppp_ball.so -------------------------

# Load lib with fallback: libppp_ball.(so|dylib|dll)
_BALL_LIB_PATH   = _resolve_lib("PPP_BALL_LIB", "ppp_ball")
_BALL = ct.CDLL(_BALL_LIB_PATH)

# Signatures
# NOTE: Keeping your radius signature (d, lam, M) since your C may expect it.
_BALL.ppp_compute_radius.argtypes = [ct.c_int, ct.c_double, ct.c_int]
_BALL.ppp_compute_radius.restype = ct.c_double

# Canonical generate signature: (lambda, M, d, seed, out*)
_BALL.ppp_generate.argtypes = [
    ct.c_int, ct.c_double, ct.c_int, ct.c_uint, ct.POINTER(ct.c_double)
]
_BALL.ppp_generate.restype = ct.c_int


def give_ball_radius(d: int, lam: float, M: int) -> float:
    """Return radius (double) computed by the C library."""
    return float(_BALL.ppp_compute_radius(ct.c_int(d), ct.c_double(lam), ct.c_int(M)))


def generate(d: int, lam: float, M: int, seed: int) -> np.ndarray:
    """
    Generate M points in R^d (float64). Returns array with shape (M, d).
    """
    pts = np.empty((M, d), dtype=np.float64)
    rc = _BALL.ppp_generate(
        ct.c_int(d),
        ct.c_double(lam),           # lambda
        ct.c_int(M),
        ct.c_uint(seed),
        pts.ctypes.data_as(ct.POINTER(ct.c_double)),
    )
    if rc != 0:
        raise RuntimeError(f"generate failed (rc={rc})")
    return pts


# --------------------- initial_proj (libinitial_proj.so) ---------------------

# Load lib with fallback: initial_proj.(so|dylib|dll)
_PROJ_LIB_PATH   = _resolve_lib("PPP_PROJ_LIB", "initial_proj")
_PROJ = ct.CDLL(_PROJ_LIB_PATH)

# int proj_select_from_points(const double* Points, int N, int d, int* out_indices)
_PROJ.proj_select_from_points.argtypes = [
    ct.POINTER(ct.c_double), ct.c_int, ct.c_int, ct.POINTER(ct.c_int)
]
_PROJ.proj_select_from_points.restype = ct.c_int


def projectVV(points: np.ndarray) -> np.ndarray:
    """
    Call C driver to select indices.
    points: (N, d) float64, C-contiguous.
    returns: (d+1,) int32 indices.
    """
    points = np.asarray(points, dtype=np.float64, order="C")
    if points.ndim != 2:
        raise ValueError("points must be 2D (N, d)")
    N, d = map(int, points.shape)
    out = np.empty(d + 1, dtype=np.int32)

    rc = _PROJ.proj_select_from_points(
        points.ctypes.data_as(ct.POINTER(ct.c_double)),
        ct.c_int(N),
        ct.c_int(d),
        out.ctypes.data_as(ct.POINTER(ct.c_int)),
    )
    if rc != 0:
        raise RuntimeError(f"proj_select_from_points failed (rc={rc})")
    return out


# --------------------- VAedge wrapper (libVAedge.so) ---------------------
# ctypes-based Python wrapper for the VAedge C API.

import os
import ctypes as ct
import numpy as np

# Load lib with fallback: VA_edge.(so|dylib|dll)
_VA_LIB_PATH = _resolve_lib("PPP_VAEDGE_LIB",  "VA_edge")  
_VA = ct.CDLL(_VA_LIB_PATH)

# ------------------------------- signatures ------------------------------
# int VA_FindNewVertex(
#   const double* Points, unsigned long long N, int d,
#   const int* Aa, int Aa_len,
#   const char* order,
#   int* out_Aa_new,
#   double* out_NV,
#   double* out_NR,
#   int* out_k,
#   int* out_left,
#   int* out_new)
_VA.VA_FindNewVertex.argtypes = [
    ct.POINTER(ct.c_double), ct.c_ulonglong, ct.c_int,
    ct.POINTER(ct.c_int), ct.c_int,
    ct.c_char_p,
    ct.POINTER(ct.c_int),
    ct.POINTER(ct.c_double),
    ct.POINTER(ct.c_double),
    ct.POINTER(ct.c_int),
    ct.POINTER(ct.c_int),
    ct.POINTER(ct.c_int),
]
_VA.VA_FindNewVertex.restype = ct.c_int

# int VA_RunLoop(
#   const double* Points, unsigned long long N, int d,
#   const int* Aa_init,
#   const char* order,
#   int   max_iter,
#   const char* log_path,
#   int*  out_Aa_final,
#   double* out_p_last,
#   double* out_r_last,
#   int*  out_iters,
#   double* out_path_length,
#   int*  out_K_sum,
#   int*  out_k_last,
#   int*  out_hit_inside,
#   int*  out_num_pairs,
#   int*  out_pairs_flat)
_VA.VA_RunLoop.argtypes = [
    ct.POINTER(ct.c_double), ct.c_ulonglong, ct.c_int,
    ct.POINTER(ct.c_int),
    ct.c_char_p,
    ct.c_int,
    ct.c_char_p,
    ct.POINTER(ct.c_int),
    ct.POINTER(ct.c_double),
    ct.POINTER(ct.c_double),
    ct.POINTER(ct.c_int),
    ct.POINTER(ct.c_double),
    ct.POINTER(ct.c_int),
    ct.POINTER(ct.c_int),
    ct.POINTER(ct.c_int),
    ct.POINTER(ct.c_int),
    ct.POINTER(ct.c_int),
]
_VA.VA_RunLoop.restype = ct.c_int

# Optional utility signatures (available from the C API)
# Bound functions below for convenience.

# int VA_CircumAll(const double* Pt, int d, double* center_out, double* R_out)
_VA.VA_CircumAll.argtypes = [
    ct.POINTER(ct.c_double), ct.c_int,
    ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
]
_VA.VA_CircumAll.restype = ct.c_int

# int VA_IsInsideSimplex_Bary(const double* point, const double* simplex, int d)
_VA.VA_IsInsideSimplex_Bary.argtypes = [
    ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_int
]
_VA.VA_IsInsideSimplex_Bary.restype = ct.c_int

# ------------------------------- helpers ---------------------------------

def _as_c_double_ptr(arr: np.ndarray) -> ct.POINTER(ct.c_double):
    return arr.ctypes.data_as(ct.POINTER(ct.c_double))


def _as_c_int_ptr(arr: np.ndarray) -> ct.POINTER(ct.c_int):
    return arr.ctypes.data_as(ct.POINTER(ct.c_int))


def _ensure_points(points: np.ndarray) -> tuple[np.ndarray, int, int]:
    points = np.asarray(points, dtype=np.float64, order="C")
    if points.ndim != 2:
        raise ValueError("points must be 2D (N, d)")
    N, d = map(int, points.shape)
    return points, N, d


def _ensure_indices(a: np.ndarray, length: int) -> np.ndarray:
    a = np.asarray(a, dtype=np.int32, order="C")
    if a.ndim != 1 or int(a.size) != length:
        raise ValueError(f"indices must be 1D with length={length}")
    return a


# ------------------------------- API ------------------------------------

def VA_find_new_vertex(points: np.ndarray, Aa: np.ndarray, order: str = "max") -> dict:
    """One step of the VAedge algorithm.

    Parameters
    ----------
    points : (N, d) float64, C-contiguous
    Aa     : (d+1,) int32 simplex indices
    order  : "max" or "min"

    Returns
    -------
    dict with keys: Aa_new (int32[d+1]), NV (float64[d]), NR (float), k (int), left (int), new (int)
    """
    points, N, d = _ensure_points(points)
    Aa = _ensure_indices(Aa, d + 1)

    if order not in ("max", "min"):
        raise ValueError("order must be 'max' or 'min'")

    out_Aa_new = np.empty(d + 1, dtype=np.int32)
    out_NV = np.empty(d, dtype=np.float64)
    out_NR = ct.c_double()
    out_k = ct.c_int()
    out_left = ct.c_int()
    out_new = ct.c_int()

    rc = _VA.VA_FindNewVertex(
        _as_c_double_ptr(points), ct.c_ulonglong(N), ct.c_int(d),
        _as_c_int_ptr(Aa), ct.c_int(d + 1),
        order.encode("ascii"),
        _as_c_int_ptr(out_Aa_new),
        _as_c_double_ptr(out_NV),
        ct.byref(out_NR),
        ct.byref(out_k),
        ct.byref(out_left),
        ct.byref(out_new),
    )
    if rc != 0:
        raise RuntimeError(f"VA_FindNewVertex failed (rc={rc})")

    return {
        "Aa_new": out_Aa_new,
        "NV": out_NV,
        "NR": float(out_NR.value),
        "k": int(out_k.value),
        "left": int(out_left.value),
        "new": int(out_new.value),
    }


def VA_run(
    points: np.ndarray,
    Aa_init: np.ndarray,
    order: str = "max",
    max_iter: int = 100000,
    log_path: str | None = None,
) -> dict:
    """Run the VAedge loop until the center lies inside the simplex or limits are hit.

    Returns dict with keys:
      Aa_final (int32[d+1]), p_last (float64[d]), r_last (float),
      iters (int), path_length (float), K_sum (int), k_last (int), hit_inside (int),
      pairs (int32[num_pairs, 2])
    """
    points, N, d = _ensure_points(points)
    Aa_init = _ensure_indices(Aa_init, d + 1)

    if order not in ("max", "min"):
        raise ValueError("order must be 'max' or 'min'")
    if max_iter < 0:
        raise ValueError("max_iter must be >= 0")

    Aa_final = np.empty(d + 1, dtype=np.int32)
    p_last = np.empty(d, dtype=np.float64)
    r_last = ct.c_double()
    iters = ct.c_int()
    path_length = ct.c_double()
    K_sum = ct.c_int()
    k_last = ct.c_int()
    hit_inside = ct.c_int()
    num_pairs = ct.c_int()

    # Preallocate flat pairs buffer length 2*max_iter (each swap contributes 2 ints)
    pairs_flat = np.empty(2 * max_iter if max_iter > 0 else 2, dtype=np.int32)

    rc = _VA.VA_RunLoop(
        _as_c_double_ptr(points), ct.c_ulonglong(N), ct.c_int(d),
        _as_c_int_ptr(Aa_init),
        order.encode("ascii"),
        ct.c_int(max_iter),
        (ct.c_char_p(log_path.encode("utf-8")) if log_path else None),
        _as_c_int_ptr(Aa_final),
        _as_c_double_ptr(p_last),
        ct.byref(r_last),
        ct.byref(iters),
        ct.byref(path_length),
        ct.byref(K_sum),
        ct.byref(k_last),
        ct.byref(hit_inside),
        ct.byref(num_pairs),
        _as_c_int_ptr(pairs_flat),
    )
    if rc != 0:
        raise RuntimeError(f"VA_RunLoop failed (rc={rc})")

    npairs = int(num_pairs.value)
    pairs = pairs_flat[: 2 * npairs].reshape(npairs, 2).copy()

    return {
        "Aa_final": Aa_final,
        "p_last": p_last,
        "r_last": float(r_last.value),
        "iters": int(iters.value),
        "path_length": float(path_length.value),
        "K_sum": int(K_sum.value),
        "k_last": int(k_last.value),
        "hit_inside": int(hit_inside.value),
        "pairs": pairs,
    }


def VA_CircumAll(Pt: np.ndarray) -> tuple[np.ndarray, float]:
    """Compute circumcenter and radius of a d-simplex.

    Parameters
    ----------
    Pt : ((d+1), d) float64, C-contiguous — simplex vertices as rows

    Returns
    -------
    (center, R) where center is float64[d], R is float
    """
    Pt = np.asarray(Pt, dtype=np.float64, order="C")
    if Pt.ndim != 2:
        raise ValueError("Pt must be 2D ((d+1), d)")
    m, d = Pt.shape
    if m != d + 1:
        raise ValueError("Pt must have shape ((d+1), d)")
    center = np.empty(d, dtype=np.float64)
    R = ct.c_double()
    rc = _VA.VA_CircumAll(
        _as_c_double_ptr(Pt), ct.c_int(d),
        _as_c_double_ptr(center), ct.byref(R),
    )
    if rc != 0:
        raise RuntimeError(f"VA_CircumAll failed (rc={rc})")
    return center, float(R.value)


def VA_inside_simplex(point: np.ndarray, simplex: np.ndarray) -> int:
    """Test if a point lies inside a simplex using barycentric coordinates.

    Parameters
    ----------
    point   : (d,) float64
    simplex : ((d+1), d) float64

    Returns
    -------
    int: 1 (inside), 0 (outside), -1 (singular/invalid)
    """
    point = np.asarray(point, dtype=np.float64, order="C")
    simplex = np.asarray(simplex, dtype=np.float64, order="C")
    if point.ndim != 1:
        raise ValueError("point must be 1D (d,)")
    if simplex.ndim != 2:
        raise ValueError("simplex must be 2D ((d+1), d)")
    d = int(point.size)
    if simplex.shape != (d + 1, d):
        raise ValueError("simplex must have shape ((d+1), d) matching point dimension")
    rc = _VA.VA_IsInsideSimplex_Bary(
        _as_c_double_ptr(point), _as_c_double_ptr(simplex), ct.c_int(d)
    )
    return int(rc)
