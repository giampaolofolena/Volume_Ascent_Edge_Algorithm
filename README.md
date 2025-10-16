## **Volume Ascent Edge Algorithm**

Voronoi-edge walk to a stable Voronoi vertex (stable Delaunay circumcenter) via iterative point swaps. Steps increase the distance to the current set of closest points.

This repository contains three **three main algorithms**:

---

### **1. Sampling**

* Uniformly sample **`M` points** in a **spherical shell** centered at the origin in `d` dimensions.
* Each point represents an **obstacle center**.

---

### **2. Initial Projection**

* **Find the closest center** to the origin:

  $\mathbf{p}_0 =\arg\min_i|\mathbf{o}_i|$

* Set the **initial velocity**:

  $\hat{\mathbf{v}}_0 = -\frac{\mathbf{p}_0}{|\mathbf{p}_0|}$
  
* **Project** the motion from the origin through $d$ steps:

  * At each step, select the **next closest point** in the direction of motion.
  * **Update** position and **project** velocity orthogonally to previously visited points.
  * After $d$ steps, reach the first **Voronoi vertex (VV)**, the **circumcenter** of $d+1$ points.

---

### **3. Edge Walk (VA-edge)**

* If the current VV is **not inside** its simplex:

  * For each point $p_j$ in the simplex:

    * Compute the **circumcenter and radius** of the $(d-1)$-simplex excluding $p_j$.
    * Determine the **Voronoi edge vector** $e_j$.
    * Check if this edge is a valid **expansion direction**.
  * **Select the edge** with the **maximal (or minimal) expansion ratio**:
    
    $k_{\text{max}} = \arg\max_j \frac{|e_j|}{r_j}$
    
  * Walk along this edge until intersecting a new point.
  * **Update** the simplex by replacing $p_{k_{\text{max}}}$ with the new point.
  * **Repeat** until the VV lies inside its simplex → **stable configuration** (an IS).

---

### Examples
<img src="media/GradientDescentB.gif" width="300" alt="VA-max">

*Geometry of one VA-max step in `d=2`*

<img src="media/Basins.png" width="300" alt="Basin">

*Example of Delaunay Basins in `d=2`*

---

## Layout

- `cpp/ppp_ball.cpp` -- Sampling of Poisson point process in a d-ball
- `cpp/initial_proj.cpp` -- selection of initial `d+1` indices
- `cpp/VA_edge.cpp` -- VA-edge algorithm (uses Eigen and nanoflann)
- `python/VA/VA.py` -- ctypes loader and wrappers

Dependencies are header-only ([Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), [nanoflann](https://github.com/jlblancoc/nanoflann)) and included as submodules under `cpp/third_party/`.

## Build Outputs

### Shared Libraries

- `build/libppp_ball.dylib`
- `build/libinitial_proj.dylib`
- `build/libVA_edge.dylib`

### Executables

- `build/ppp_ball_cli`
- `build/initial_proj_cli`
- `build/VA_edge_cli`

On Linux, Windows and macOS, shared libraries are built as `.so`,`.dll` and `.dylib` respectively.

## Quick start

```bash
git clone https://github.com/giampaolofolena/Volume_Ascent_Edge_Algorithm
cd Volume_Ascent_Edge_Algorithm
git submodule update --init --recursive

cmake -S . -B build
cmake --build build -j
```

---

## Python example

Minimal end-to-end usage with python.

```bash
cd python; python
```

```python
from VA import *

# Generate N points in d=3, inner radius=1.0, seed=1
pts = generate(3, 1, 100, 1)  # (N, d) float64

# Select initial (d+1) indices of a Voronoi vertex via initial projection
idx = projectVV(pts)          # (d+1,) int32

# Run the Voronoi-edge (volume-ascent) walk until a stable vertex or limits reached
res = VA_run(pts, idx)
print(res)
```

### Inputs

* `generate(d, r_in, N, seed)`

  * `d`: dimension.
  * `r_in`: inner radius of spherical shell.
  * `N`: number of points.
  * `seed`: RNG seed.
* `projectVV(points)`

  * `points` `(N×d, float64)` C-contiguous.
  * Returns `(d+1,) int32` indices.
* `VA_run(points, idx)`

  * `points` `(N×d, float64)`, `idx` `((d+1,), int32)`.
  * Returns a dict:

    * `Aa_final`: `((d+1,), int32)` — final simplex indices.
    * `p_last`: `(d,) float64` — final circumcenter.
    * `r_last`: `float` — final radius.
    * `iters`: `int` — iterations used.
    * `path_length`: `float` — path length of center.
    * `K_sum`, `k_last`: search diagnostics.
    * `hit_inside`: `int` — 1 if final center lies inside simplex.
    * `pairs`: `(num_swaps×2, int32)` — vertex swaps `(left -> new)`.

### Example session output

```
>>> pts = generate(3,1,1000,1)
>>> pts
array([[-0.04788522, -0.00893964, -0.07050876],
       [ 0.10402446,  0.07928355,  0.21671649],
       [-0.75496705,  0.038947  ,  0.19303341],
       ...,
       [ 3.97784881,  3.83475923, -2.80246953],
       [ 2.84458129,  1.58805878, -5.2696978 ],
       [ 1.58798733,  5.42454504, -2.55311521]])
>>> idx = projectVV(pts)
CHECK CONSISTENCY: distances of x from Aa *points* (they should be equal up to numerical error)
0.877092 0.877092 0.877092 0.877092
>>> idx
array([ 0,  1,  8, 10], dtype=int32)
>>> VA_run(pts,idx)
{'Aa_final': array([21, 17, 16, 64], dtype=int32), 'p_last': array([ 0.86489715, -0.90778117,  0.27611627]), 'r_last': 1.094844111034065, 'iters': 4, 'path_length': 0.504963121674542, 'K_sum': 6, 'k_last': 1, 'hit_inside': 1, 'pairs': array([[ 0, 21],
       [ 8, 16],
       [ 1, 17],
       [10, 64]], dtype=int32)}
>>>
```

---

## Executable example (CLI)

```bash
# 1) Generate points into a text file
./ppp_ball_cli 3 1 100 1 > points.txt

# 2) Select initial Voronoi-vertex indices
./initial_proj_cli points.txt 3 > indexes.txt

# 3) Run VA-edge (order can be 'max' or 'min')
./VA_edge_cli points.txt indexes.txt 3 max
```

### CLI I/O

* `ppp_ball_cli d r_in N seed` → writes `N` lines with `d` floats to STDOUT.
* `initial_proj_cli points.txt d` → prints `(d+1)` indices to STDOUT (space-separated).
* `VA_edge_cli points.txt indexes.txt d order` → prints final simplex, center, radius, diagnostics, swap pairs.

---

**Note**: Full algorithmic details are in the paper [Jamming the Random Lorentz Gas](https://arxiv.org/abs/2410.05784).
