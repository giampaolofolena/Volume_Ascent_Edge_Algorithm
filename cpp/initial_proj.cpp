// initial_proj.cpp : 
// Projection toward the first (closest along the direction) Voronoi vertex
//   (i.e., a Delaunay circumcenter) via iterative face hits.
//
// Pipeline
//   0) Read flat input, build Eigen row vectors, and sort rows by ||row|| (closest first).
//   1) Translate by the closest point (now row 0); set x0 = 0 and v0 = -p0/||p0||; build A0 (unit directions) and b0 (distances).
//   2) Run d InitialProjection steps; collect active face (points) indices in Aa.
//
// I/O
//   • Library: int proj_select_from_points(const double* Points, int N, int d, int* out_indices)
//       - Points: row-major (N×d)
//       - out_indices[0..d]: indices in the sorted/translated frame {0, 1+Aa[0], …, 1+Aa[d-1]}
//   • CLI: reads a space-delimited file and d; prints indices only (one line).
//          ./initial_proj <points.txt> <dim> 

// Conventions
//   - Types: int and double only. Errors are returned as integer codes (no exceptions).
//
// Build
//   Exec:   g++ -O3 -std=c++17 initial_proj.cpp -o initial_proj -DBUILD_EXECUTABLE -Ithird_party/eigen/
//   Shared: g++ -O3 -std=c++17 -shared -fPIC initial_proj.cpp -o libinitial_proj.so -DBUILD_LIB -Ithird_party/eigen


#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;


/**
 * Build std::vector<Eigen::VectorXd> directly from a flat buffer (row-major N x d)
 * and sort the vectors by increasing squared norm.
 */
std::vector<Eigen::VectorXd> sort_points_by_norm_from_origin(const double* Points, int N, int d)
{
    if (!Points || N <= 0 || d <= 0) return {};

    std::vector<Eigen::VectorXd> pts;
    pts.reserve(N);

    // Directly push rows as vectorXd without mapping as a matrix
    for (int i = 0; i < N; ++i) {
        Eigen::VectorXd v(d);
        const double* src = Points + i * d;
        for (int j = 0; j < d; ++j) {
            v(j) = src[j];
        }
        pts.push_back(std::move(v));
    }

    // Sort ONLY if not already sorted by increasing squared norm
    auto cmp = [](const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
        return a.squaredNorm() < b.squaredNorm();
    };
    const bool already_sorted = std::is_sorted(pts.begin(), pts.end(), cmp);
    if (!already_sorted) {
        std::sort(pts.begin(), pts.end(), cmp);
        std::cerr << "points have been sorted based on the distance from the origin (the returned indexes correspond to the sorted set) ";
        for (int i = 0; i < N; ++i) {
                std::cout << pts[i].squaredNorm() << ' ' ;
            }
        std::cout << ' ';
    }

    return pts;
}


// Step 1 (+ A0,b0): pts must be sorted by norm, first is origin.
// - Translates all points by subtracting pts[0] in-place.
// - Sets x0 = 0, v0 = -p0/||p0||.
// - Builds A0 (unit directions of translated points i>=1) and b0 (their norms).
// Returns 0 on success; nonzero on bad input.
int reset_origin_init_and_build(std::vector<Eigen::VectorXd>& pts,
                                Eigen::VectorXd& p0,
                                Eigen::VectorXd& x0,
                                Eigen::VectorXd& v0,
                                Eigen::MatrixXd& A0,
                                Eigen::VectorXd& b0)
{
    if (pts.empty()) return 1;
    int d = (int)pts[0].size();
    for (const auto& v : pts) if ((int)v.size() != d) { std::cerr << "points do not have expected dimension"; return 2; };

    // p0 = first point (closest), set x0 and v0
    p0 = pts[0];
    x0 = -p0;
    double n0 = x0.norm();
    if(n0 > 0.0) { v0 = x0 / n0; } else { std::cerr << "the norm of the closest vector is negative"; return 2; };

    // translate all points: pts[i] -= p0
    for (auto& v : pts) v.noalias() -= p0;

    // Build A0 and b0 from translated points (skip row 0 which is now zero)
    int N = (int)pts.size();
    int m = N - 1;
    A0.resize(m, d);
    b0.resize(m);
    for (int i = 0; i < m; ++i) {
        const Eigen::VectorXd& q = pts[i + 1];
        double n = q.norm();
        b0(i) = n/2.;
        if (n > 0.0) A0.row(i) = (q / n).transpose();
        else         A0.row(i).setZero();
    }
    return 0;
}

// Given A (m x d), b (m), x0 (d), dir (d), and active set Aa,
// pick md = argmin_i (b_i - A_i·x0) / (A_i·dir) over valid i.
// Valid means: i ∉ Aa, (A_i·dir) > 0, distance finite and ≥ 0.
int pick_closest_face(const Eigen::MatrixXd& A,
                      const Eigen::VectorXd& b,
                      const std::vector<int>& Aa,
                      const Eigen::VectorXd& x,
                      const Eigen::VectorXd& dir,
                      int& md, double& step)
{
    const int m = (int)A.rows();
    const double INF = std::numeric_limits<double>::infinity();

    // Fast active flags
    std::vector<char> active(m, 0);
    for (int i : Aa) if (0 <= i && i < m) active[i] = 1;

    md = -1;
    step = INF;

    for (int i = 0; i < m; ++i) {
        if (active[i]) continue;                 // already active
        const double den = A.row(i).dot(dir);
        if (den <= 0.0) continue;                // not forward-facing
        const double num = b(i) - A.row(i).dot(x);
        const double ti  = num / den;
        if (ti >= 0.0 && std::isfinite(ti) && ti < step) {
            step = ti;
            md = i;
        }
    }
    return (md >= 0) ? 1 : 0; // 1 = ok, 0 = no feasible face
}


// InitialProjection: single step using pick_closest_face
// Returns 0 on success; 1 dim mismatch; 2 no feasible face.
int InitialProjection(const Eigen::MatrixXd& A,            // m x d
                      const Eigen::VectorXd& b,            // m
                      std::vector<int>& Aa,                // active rows
                      const Eigen::VectorXd& x,            // d
                      const Eigen::VectorXd& dir,          // d
                      Eigen::VectorXd& x_out,              // d
                      Eigen::VectorXd& v_out)              // d
{
    const int m = (int)A.rows();
    const int d = (int)A.cols();
    if (b.size() != m || x.size() != d || dir.size() != d) return 1;

    // Pick closest feasible face and step size
    int md; double t;
    if (!pick_closest_face(A, b, Aa, x, dir, md, t)) { std::cerr << "could not pick the closest face"; return 2; };

    // Step
    x_out = x + t * dir;

    // Activate the chosen face
    Aa.push_back(md);

    // v = direction projected onto orthogonal complement of span(A[Aa]^T)
    Eigen::MatrixXd M(d, (int)Aa.size());
    for (int c = 0; c < (int)Aa.size(); ++c)
        M.col(c) = A.row(Aa[c]).transpose();

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(M);
    const int rnk = qr.rank();
    if (rnk > 0) {
        Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(d, rnk);
        v_out = dir - Q * (Q.transpose() * dir);
    } else {
        v_out = dir;
    }

    // Keep direction sign consistent with x
    if (v_out.dot(x) < 0.0) v_out = -v_out;
    return 0;
}

// Euclidean distances from x to the Aa *points* (not faces).
// Each point q_i in translated coords is: q_i = b0(i) * A0.row(i).transpose()
// Output: out[t] = || x - q_{Aa[t]} ||, length = Aa.size()
// Returns 0 on success; 1 on dim mismatch.
int distances_to_Aa_points(const std::vector<Eigen::VectorXd>& pts,   // points
                           const std::vector<int>& out_idx,     // active indices
                           const Eigen::VectorXd& x)            // d   (current point, translated frame)
{
    std::cerr << "CHECK CONSISTENCY: distances of x from Aa *points* (they should be equal up to numerical error)" << std::endl;
    for (int t = 0; t < (int)out_idx.size(); ++t) {
        int i = out_idx[t];
        std::cerr << (x - pts[i]).norm() << " ";
    }
    std::cerr << '\n';
    
    return 0;
}


// =========================== Driver (d steps) ================================
// High-level flow for selection:
//   0) Convert flat buffer -> Eigen vectors; sort by ||row|| (ascending).
//   1) Translate by the closest point; set x0 and v0; build A0, b0.
//   2) Run d InitialProjection steps; Aa collects chosen face indices.
//   Output indices are in the sorted/translated frame:
//   out_idx[0] = 0 (closest point), out_idx[i+1] = 1 + Aa[i].

static int select_indices_from_points(const double* Points, int N, int d,
                                      std::vector<int>& out_idx)
{
    if (!Points || N < 2 || d <= 0) return 1;

    // 0) sort points by norm (closest first) and convert to Eigen vectors
    auto pts = sort_points_by_norm_from_origin(Points, N, d);
    if ((int)pts.size() != N) return 1;

    // 1) reset origin, initial velocity, and build A0, b0
    Eigen::VectorXd p0, x0, v0;
    Eigen::MatrixXd A0;
    Eigen::VectorXd b0;
    if (reset_origin_init_and_build(pts, p0, x0, v0, A0, b0) != 0) { std::cerr << "the norm of the closest vector is negative"; return 2; };

    // 2) run d Initial Projection steps
    std::vector<int> Aa; Aa.reserve(d);
    Eigen::VectorXd x = x0;
    Eigen::VectorXd v = v0;

    for (int s = 0; s < d; ++s) {
        Eigen::VectorXd xn, vn;
        int rc = InitialProjection(A0, b0, Aa, x, v, xn, vn);
        if (rc != 0) return rc;
        x = std::move(xn);
        v = std::move(vn);
    }

    // 3) indices in the (sorted & translated) frame:
    //     0 is the closest point (origin after translation),
    //     A0 rows correspond to pts[1..N-1], so add 1.
    out_idx.resize(d + 1);
    out_idx[0] = 0;
    for (int i = 0; i < d; ++i) out_idx[i + 1] = 1 + Aa[i];

    // >>> CHECK CONSISTENCY: compute distances from x to the Aa *points* and print them (they should be equal up to numerical error)
    distances_to_Aa_points(pts,out_idx,x);
    // <<< END CHECK

    return 0;
}




// ============================ LIB API ========================================
#if defined(BUILD_LIB)
extern "C" {

// Write (d+1) selected indices into out_indices[0..d].
int proj_select_from_points(const double* Points, int N, int d, int* out_indices)
{
    if (!Points || !out_indices) return 3;
    std::vector<int> idx;
    int rc = select_indices_from_points(Points, N, d, idx);
    if (rc != 0) return rc;
    for (int i = 0; i < d + 1; ++i) out_indices[i] = idx[i];
    return 0;
}

} // extern "C"
#endif // BUILD_LIB

// =========================== EXECUTABLE MAIN =================================
#if defined(BUILD_EXECUTABLE)
int main(int argc, char** argv)
{
    // Minimal CLI: read "points.txt" and dimension d; print only indices.
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <points.txt> <d>";
        return 1;
    }
    const char* path = argv[1];
    int d = std::stoi(argv[2]);

    std::ifstream fin(path);
    if (!fin) { std::cerr << "Cannot open: " << path << ""; return 1; }

    // Read data into a flat buffer Praw (row-major). Validate row width.
    std::vector<double> Praw;
    std::string line;
    int cols = -1;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::vector<double> row; double x;
        while (iss >> x) row.push_back(x);
        if (cols < 0) cols = (int)row.size();
        if ((int)row.size() != cols) { std::cerr << "Inconsistent row width"; return 1; }
        Praw.insert(Praw.end(), row.begin(), row.end());
    }
    if (cols != d) { std::cerr << "File has d=" << cols << " but provided d=" << d << ""; return 1; }
    if (Praw.empty()) { std::cerr << "Empty input"; return 1; }

    int N = (int)(Praw.size() / d);
    if (N < d + 1) { std::cerr << "Need at least d+1 points"; return 1; }

    std::vector<int> idx;
    int rc = select_indices_from_points(Praw.data(), N, d, idx);
    if (rc != 0) { std::cerr << "Projection failed (rc=" << rc << ")"; return 1; }

    // Print only indices, single line.
    for (int i = 0; i < d + 1; ++i) {
        if (i) std::cout << ' ';
        std::cout << idx[i];
    }
    std::cout << ' ';
    return 0;
}
#endif // BUILD_EXECUTABLE
