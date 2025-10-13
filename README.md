## **Volume Ascent Edge Algorithm**

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

- `cpp/ppp_ball.cpp` -- Poisson point process in a d-ball (target number of points ? M)
- `cpp/initial_proj.cpp` -- selection of initial `d+1` indices
- `cpp/VAedge.cpp` -- VA-edge loop (uses Eigen and nanoflann)
- `python/VA/VA.py` -- ctypes loader and wrappers

Dependencies are header-only (Eigen, nanoflann) and included as submodules under `cpp/third_party/`.

## Build Outputs

### Shared Libraries

- `build/libppp_ball.so`
- `build/libinitial_proj.so`
- `build/libVAedge.so`

### Executables

- `build/ppp_ball_cli`
- `build/initial_proj_cli`
- `build/VAedge_cli`

On Windows and macOS, shared libraries are built as `.dll` and `.dylib` respectively.

## Quick start

```bash
git clone https://github.com/giampaolofolena/Volume_Ascent_Edge_Algorithm
cd Volume_Ascent_Edge_Algorithm
git submodule update --init --recursive

cmake -S . -B build
cmake --build build -j
```

# Python example
PYTHONPATH=python:$PYTHONPATH python python/examples/generate_and_run.py


---

**Note**: Full algorithmic details are in the paper [Jamming the Random Lorentz Gas](https://arxiv.org/abs/2410.05784).
