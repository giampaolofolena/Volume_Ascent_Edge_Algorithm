// VAedge.cpp :
// Volume-Ascent Voronoi edge walk to a stable Voronoi vertex (stable Delaunay circumcenter)
// via iterative point swaps (i.e., moves along Voronoi edges). Steps are constrained to
// monotonically increase the distance to the current set of closest points.
// In the 'max' case, the chosen face maximizes r / ||V - c||, giving the largest radius increase.

// Algorithm (high level)
// 0) Build a KD-tree over the row-major point set (N×d).
// 1) For the current simplex Aa (d+1 indices), compute its circumcenter V and radius R.
// 2) For each face of Aa, compute facet circumcenter c and radius r; score faces by r/||V−c|| and pick a face by `order` (min|max).
// 3) Move along v = (V−c)/||V−c|| and search a new point entering through that face by expanding radius queries.
// 4) Swap the dropped vertex with the found point; update V, R; repeat until V lies inside the simplex or limits are reached.
//
// I/O
//  • Library (C API):
//      - int VA_FindNewVertex(...): one face-swap step.
//      - int VA_RunLoop(...): iterate steps until termination; optional log.
//      - int VA_EqOfPlanes(...), int VA_Find_Direction(...), int VA_CircumAll(...),
//      - int VA_CircumPart(...), int VA_IsInsideSimplex_Bary(...): geometric utilities.
//      Conventions: row-major doubles, ints for sizes and indices. Error codes as ints.
//  • CLI (when BUILD_EXECUTABLE):
//      ./VAedge <points.txt> <indexes.txt> <dim> <order:min|max> [max_iter] [log_path]
//
// Conventions
// - No exceptions across the C API; return ints.
// - Points are immutable input; indices refer to original rows.
// - Uses Eigen for linear algebra and nanoflann for radius queries.
//
// Build
// Shared lib: g++ -O3 -std=c++17 -DNDEBUG -fPIC -shared VA_edge.cpp -o libVA_edge.so -DBUILD_LIB -ffast-math -Ithird_party/eigen -Ithird_party/nanoflann/include
// Executable: g++ -O3 -std=c++17 -DNDEBUG VA_edge.cpp -o VA_edge -DBUILD_EXECUTABLE -ffast-math -Ithird_party/eigen -Ithird_party/nanoflann/include

#include <Eigen/Dense>
#include <nanoflann.hpp>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

using Eigen::ColPivHouseholderQR;
using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::VectorXd;

// --------------------------- nanoflann adaptor ---------------------------

// Simple row-major point cloud adaptor around a raw double* (N x d)
struct RowMajorPoints {
    const double* pts = nullptr; // pointer to first element (row-major)
    size_t N = 0;                // number of points
    int d = 0;                   // dimensionality

    inline size_t kdtree_get_point_count() const { return N; }
    inline double kdtree_get_pt(size_t idx, size_t dim) const {
        return pts[idx * (size_t)d + (size_t)dim];
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

using KDAdaptor = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, RowMajorPoints>, RowMajorPoints, -1>;

// nanoflann version differences: index/dist accessors
template <typename T>
inline auto nlfi_idx(const T& m) -> decltype(m.index) { return m.index; }
template <typename T>
inline auto nlfi_idx(const T& m) -> decltype(m.first) { return m.first; }

template <typename T>
inline auto nlfi_d2(const T& m) -> decltype(m.dist_sq) { return m.dist_sq; }
template <typename T>
inline auto nlfi_d2(const T& m) -> decltype(m.second) { return m.second; }

// Pick the ResultItem type that matches the nanoflann version
template <typename KD>
struct NlfiResultItem {
#if defined(NANOFLANN_VERSION) && (NANOFLANN_VERSION >= 0x0150)
    using type = nanoflann::ResultItem<typename KD::IndexType, typename KD::DistanceType>;
#else
    using type = nanoflann::ResultItem<typename KD::IndexType>;
#endif
};

// --------------------------- small utilities ----------------------------

static inline double Vball(int d) { return std::pow(M_PI, 0.5 * d) / std::tgamma(0.5 * d + 1.0); }
static inline double phiJ(int d, double x) { return (d > 0) ? Vball(d) * std::pow(x, d) / double(d) : 0.0; }

static std::vector<int> stable_argsort_ascending(const std::vector<double>& a) {
    std::vector<int> idx(a.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&](int i, int j) { return a[i] < a[j]; });
    return idx;
}

// Read a simple indices file (space/newline separated ints). Optional header line starting with '#'.
static bool read_indices_file(const std::string& path, int expected_len, std::vector<int>& idx) {
    std::ifstream fin(path);
    if (!fin) return false;
    idx.clear();

    // Skip optional header line starting with '#'
    std::string tok; std::streampos pos = fin.tellg();
    if (fin >> tok) {
        if (tok.rfind("#", 0) != 0) { fin.clear(); fin.seekg(pos); }
        else { std::string rest; std::getline(fin, rest); }
    } else return false;

    int v;
    while (fin >> v) idx.push_back(v);
    return (expected_len <= 0) || ((int)idx.size() == expected_len);
}

// Assemble a MatrixXd by gathering rows Points[indices]
static inline MatrixXd gather_rows(const double* Points, int d, const int* indices, int len) {
    MatrixXd M(len, d);
    for (int i = 0; i < len; ++i) {
        const double* src = Points + (size_t)indices[i] * (size_t)d;
        for (int j = 0; j < d; ++j) M(i, j) = src[j];
    }
    return M;
}

// ---------------------- geometric core (circumcenters) ------------------

// Given m points X (m x d) and a reference X0 (1 x d),
// build the hyperplane equations of perpendicular bisectors: A x = b
static inline void EqOfPlanes(const MatrixXd& X, const RowVectorXd& X0, MatrixXd& A, VectorXd& b) {
    const int m = (int)X.rows();
    A = 0.5 * (X.rowwise() - X0);                  // directions
    const MatrixXd xc = 0.5 * (X.rowwise() + X0);  // midpoints
    b.resize(m);
    for (int i = 0; i < m; ++i) b(i) = A.row(i).dot(xc.row(i));
}

// Circumcenter of a d-simplex given by (d+1) x d matrix Pt
static inline VectorXd Find_Direction(const MatrixXd& Pt) {
    const int d = (int)Pt.cols();
    const MatrixXd X = Pt.bottomRows(d);  // Pt[1:]
    const RowVectorXd X0 = Pt.row(0);     // Pt[0]
    MatrixXd A; VectorXd b;
    EqOfPlanes(X, X0, A, b);
    ColPivHouseholderQR<MatrixXd> qr(A);
    return qr.solve(b);
}

static inline void CircumAll(const MatrixXd& Pt, VectorXd& center, double& R) {
    center = Find_Direction(Pt);
    VectorXd diff = center - Pt.row(0).transpose();
    R = diff.norm();
}

// Circumcenter of a facet (d points in R^d); returns center c and radius r
static inline void CircumPart(const MatrixXd& Pt, VectorXd& c, double& r) {
    const int d = (int)Pt.cols();
    const int rows = (int)Pt.rows(); // expect rows == d

    const VectorXd x = Find_Direction(Pt);           // circumcenter of full (with removed vertex previously)

    // Base = Pt[1:] - Pt[0] ( (d-1) x d )
    MatrixXd Base(rows - 1, d);
    for (int i = 0; i < rows - 1; ++i) Base.row(i) = Pt.row(i + 1) - Pt.row(0);

    // Thin QR on Base^T to get an orthonormal basis of the facet subspace
    MatrixXd BaseT = Base.transpose();               // d x (d-1)
    Eigen::HouseholderQR<MatrixXd> qr(BaseT);
    MatrixXd Q = qr.householderQ() * MatrixXd::Identity(d, rows - 1); // d x (d-1)

    // Project vector (x - Pt[0]) onto facet subspace
    const VectorXd w = x - Pt.row(0).transpose();
    const VectorXd v = Q * (Q.transpose() * w);

    c = v + Pt.row(0).transpose();
    r = v.norm();
}

// ---------------------------- inside-simplex ----------------------------
// Test if point is inside simplex (barycentric). Returns: 1 inside, 0 outside, -1 singular.
static inline int is_inside_simplex_bary(const VectorXd& point, const MatrixXd& simplex) {
    const int d = (int)point.size();
    if (simplex.rows() != d + 1 || simplex.cols() != d) return -1;

    // Build augmented system to solve for barycentric coordinates
    MatrixXd A(d + 1, d + 1);
    A.topRows(d) = simplex.transpose();
    A.row(d).setOnes();

    VectorXd rhs(d + 1); rhs.head(d) = point; rhs(d) = 1.0;

    Eigen::FullPivLU<MatrixXd> lu(A);
    if (!lu.isInvertible()) return -1;
    const VectorXd bary = lu.solve(rhs);

    for (int j = 0; j < d + 1; ++j) {
        const double y = bary(j);
        if (y < 0.0 || y > 1.0) return 0;
    }
    return 1;
}

// -------------------------- FindNewVertex core ---------------------------

// Internal worker that reuses an existing KD-tree (faster inside loops)
static inline int FindNewVertex_core(
    const double* Points, unsigned long long N, int d,
    const int* Aa, int Aa_len, const std::string& order,
    KDAdaptor& kdt,
    std::vector<int>& out_Aa_new,
    VectorXd& NV, double& NR, int& k_out,
    int* out_left, int* out_new) {

    // Current simplex and circumcenter
    MatrixXd PtSel = gather_rows(Points, d, Aa, Aa_len);
    VectorXd V; double R = 0.0; // V: circumcenter of current simplex; R: radius
    CircumAll(PtSel, V, R);
    const int DIM = Aa_len - 1; // ambient simplex dimension

    k_out = 0;
    std::vector<int> I;                // candidate face indices
    std::vector<double> Rat;           // their ratios

    // Evaluate all faces of current simplex
    for (int i = 0; i < Aa_len; ++i) {
        std::vector<int> Am; Am.reserve(Aa_len - 1);
        for (int t = 0; t < Aa_len; ++t) if (t != i) Am.push_back(Aa[t]);

        MatrixXd PtAm = gather_rows(Points, d, Am.data(), DIM);
        VectorXd c; double r = 0.0; // facet circumcenter and radius
        CircumPart(PtAm, c, r);

        const double cV = (V - c).norm();
        if (cV == 0.0) continue;
        const double rat = r / cV;

        const VectorXd VmC = V - c;
        const VectorXd Pi  = Eigen::Map<const VectorXd>(Points + (size_t)Aa[i] * (size_t)d, d) - V;
        const VectorXd Pj0 = Eigen::Map<const VectorXd>(Points + (size_t)Am[0] * (size_t)d, d) - V;

        if (VmC.dot(Pi) < VmC.dot(Pj0)) { Rat.push_back(rat); I.push_back(i); ++k_out; }
    }
    if (I.empty()) return 2; // no valid face

    // Stable argsort, pick min or max ratio
    const auto ord = stable_argsort_ascending(Rat);
    const int chosen_pos = (order == "min") ? ord.back() : ord.front();
    const int i_max = I[chosen_pos];

    // Face for i_max
    std::vector<int> Amx; Amx.reserve(Aa_len - 1);
    for (int t = 0; t < Aa_len; ++t) if (t != i_max) Amx.push_back(Aa[t]);

    MatrixXd PtAmx = gather_rows(Points, d, Amx.data(), DIM);
    VectorXd cmax; double rmax = 0.0;
    CircumPart(PtAmx, cmax, rmax);

    const double Vc = (V - cmax).norm();
    if (Vc == 0.0) return 3;
    const VectorXd vdir = (V - cmax) / Vc; // direction along which the new center moves

    // Candidate search radius expansion (matches Python logic)
    double rr = (DIM != 2) ? std::pow(4.0 * DIM, 1.0 / double(DIM)) : 10.0;
    double a_min = 1000.0; // sentinel
    int    pt_min = -1;

    while (a_min == 1000.0 && rr < 101.0) {
        const double rad = R * rr;
        const double rad2 = rad * rad;

        using RI = typename NlfiResultItem<KDAdaptor>::type;
        std::vector<RI> matches;

#if defined(NANOFLANN_HAS_SEARCH_PARAMS)
        nanoflann::SearchParameters params; params.sorted = false;
        kdt.radiusSearch(V.data(), (KDAdaptor::DistanceType)rad2, matches, params);
#else
        kdt.radiusSearch(V.data(), (KDAdaptor::DistanceType)rad2, matches);
#endif

        // Sort by index to mimic Python's order
        std::sort(matches.begin(), matches.end(), [](const RI& a, const RI& b) { return nlfi_idx(a) < nlfi_idx(b); });

        // Evaluate candidates within the radius
        for (const auto& m : matches) {
            const int pi = (int)nlfi_idx(m);

            // skip current simplex points
            bool inAa = false;
            for (int t = 0; t < Aa_len; ++t) if (Aa[t] == pi) { inAa = true; break; }
            if (inAa) continue;

            const VectorXd Pv = Eigen::Map<const VectorXd>(Points + (size_t)pi * (size_t)d, d) - V;
            if (Pv.squaredNorm() > rad2) continue; // inclusive radius check

            const double Pj = vdir.dot(Pv);
            if (Pj == 0.0) continue;

            const double rd = (Pv - Pj * vdir).norm();
            const double num = rd * rd + Pj * Pj - rmax * rmax - Vc * Vc;
            const double den = 2.0 * (Pj + Vc);
            if (den == 0.0) continue;

            const double a = num / den;
            if (a > 0.0 && a < a_min) { a_min = a; pt_min = pi; }
        }
        rr *= 1.1;
    }

    if (rr > 100.0) k_out = 9999;
    if (pt_min < 0) return 4; // no entering vertex found

    // Build new simplex by replacing the chosen vertex
    out_Aa_new.assign(Aa, Aa + Aa_len);
    out_Aa_new[i_max] = pt_min;

    if (out_left) *out_left = Aa[i_max]; // removed
    if (out_new)  *out_new  = pt_min;    // added

    // Compute circumcenter for new simplex
    MatrixXd PtNew = gather_rows(Points, d, out_Aa_new.data(), Aa_len);
    CircumAll(PtNew, NV, NR);
    return 0;
}

// Convenience shim that builds a temporary KD-tree
static inline int FindNewVertex_core_shim(
    const double* Points, unsigned long long N, int d,
    const int* Aa, int Aa_len, const std::string& order,
    std::vector<int>& out_Aa_new,
    VectorXd& NV, double& NR, int& k_out,
    int* out_left, int* out_new) {

    RowMajorPoints cloud{Points, static_cast<size_t>(N), d};
    KDAdaptor kdt(d, cloud, {10});
    kdt.buildIndex();

    return FindNewVertex_core(Points, N, d, Aa, Aa_len, order,
                              kdt, out_Aa_new, NV, NR, k_out,
                              out_left, out_new);
}

// ------------------------------ run loop --------------------------------

static inline int VA_RunLoop_impl(
    const double* Points, unsigned long long N, int d,
    const int* Aa_init,
    const char* order,
    int   max_iter,
    const char* log_path,
    int*  out_Aa_final,
    double* out_p_last,
    double* out_r_last,
    int*  out_iters,
    double* out_path_length,
    int*  out_K_sum,
    int*  out_k_last,
    int*  out_hit_inside,
    int* out_num_pairs,
    int* out_pairs_flat) {

    if (!Points || !Aa_init || !out_Aa_final || !out_p_last || !out_r_last ||
        !out_iters || !out_path_length || !out_K_sum || !out_k_last || !out_hit_inside)
        return 5;

    const int Aa_len = d + 1;
    std::vector<int> Aa(Aa_len);
    for (int i = 0; i < Aa_len; ++i) Aa[i] = Aa_init[i];

    // Initial circumcenter
    MatrixXd PtSel = gather_rows(Points, d, Aa.data(), Aa_len);
    VectorXd p_curr(d); double r_curr = 0.0;
    CircumAll(PtSel, p_curr, r_curr);

    // Build KD-tree once
    RowMajorPoints cloud{Points, (size_t)N, d};
    KDAdaptor kdt(d, cloud, {10});
    kdt.buildIndex();

    int i = 0, Ksum = 0, k_last = 0;
    double L = 0.0;
    const int DIM = Aa_len - 1;

    std::ofstream flog;
    if (log_path && std::strlen(log_path) > 0) flog.open(log_path);

    MatrixXd simplex(Aa_len, d);
    *out_hit_inside = 0;

    std::vector<int> pairs;
    pairs.reserve(std::max(0, max_iter) * 2);

    while (true) {
        // build simplex and test if current center is inside
        for (int t = 0; t < Aa_len; ++t)
            for (int j = 0; j < d; ++j)
                simplex(t, j) = Points[(size_t)Aa[t] * (size_t)d + (size_t)j];

        const int inside = is_inside_simplex_bary(p_curr, simplex);
        if (inside == 1) { *out_hit_inside = 1; break; }
        if (i >= max_iter) break;
        if (k_last >= 9999) break;

        // advance one step
        int left = -1, neu = -1; std::vector<int> Aa_new(Aa_len);
        VectorXd NV(d); double NR = 0.0; int k = 0;

        const int rc = FindNewVertex_core(
            Points, N, d,
            Aa.data(), Aa_len,
            order ? order : std::string("max"),
            kdt,
            Aa_new, NV, NR, k,
            &left, &neu);

        if (rc != 0) { k_last = k; break; }

        pairs.push_back(left); pairs.push_back(neu);

        if (k > 1 && flog.good()) {
            const double dist = (NV - p_curr).norm();
            flog << i << ' ' << dist << ' ' << NR << ' ' << phiJ(DIM, NR) << ' '
                 << L << ' ' << k << ' ' << p_curr.norm() << '\n';
        }

        L += (NV - p_curr).norm();
        ++i; Ksum += k; k_last = k;
        Aa.swap(Aa_new);
        p_curr = NV; r_curr = NR;
    }

    for (int t = 0; t < Aa_len; ++t) out_Aa_final[t] = Aa[t];
    for (int j = 0; j < d; ++j) out_p_last[j] = p_curr(j);
    *out_r_last = r_curr;
    *out_iters = i; *out_path_length = L;
    *out_K_sum = Ksum; *out_k_last = k_last;

    if (out_num_pairs) *out_num_pairs = static_cast<int>(pairs.size() / 2);
    if (out_pairs_flat && out_num_pairs && *out_num_pairs > 0) {
        const int need = (*out_num_pairs) * 2;
        for (int k = 0; k < need; ++k) out_pairs_flat[k] = pairs[k];
    }
    return 0;
}

// ------------------------------ C API -----------------------------------

#if defined(BUILD_LIB) && !defined(BUILD_EXECUTABLE)
extern "C" {

int VA_FindNewVertex(
    const double* Points, unsigned long long N, int d,
    const int* Aa, int Aa_len,
    const char* order,
    int* out_Aa_new,
    double* out_NV,
    double* out_NR,
    int* out_k,
    int* out_left,
    int* out_new) {

    if (!Points || !Aa || !out_Aa_new || !out_NV || !out_NR || !out_k) return 5;

    std::vector<int> newAa; VectorXd NV(d); double NR = 0.0; int k = 0;

    const int rc = FindNewVertex_core_shim(
        Points, N, d, Aa, Aa_len,
        order ? std::string(order) : std::string("max"),
        newAa, NV, NR, k, out_left, out_new);
    if (rc != 0) return rc;

    for (int i = 0; i < Aa_len; ++i) out_Aa_new[i] = newAa[i];
    for (int j = 0; j < d; ++j) out_NV[j] = NV(j);
    *out_NR = NR; *out_k = k;
    return 0;
}

int VA_RunLoop(
    const double* Points, unsigned long long N, int d,
    const int* Aa_init,
    const char* order,
    int   max_iter,
    const char* log_path,
    int*  out_Aa_final,
    double* out_p_last,
    double* out_r_last,
    int*  out_iters,
    double* out_path_length,
    int*  out_K_sum,
    int*  out_k_last,
    int*  out_hit_inside,
    int*  out_num_pairs,
    int*  out_pairs_flat) {

    return VA_RunLoop_impl(Points, N, d, Aa_init, order, max_iter, log_path,
                           out_Aa_final, out_p_last, out_r_last, out_iters,
                           out_path_length, out_K_sum, out_k_last,
                           out_hit_inside, out_num_pairs, out_pairs_flat);
}

// EqOfPlanes: X (m x d), X0 (d) -> A (m x d), b (m)
int VA_EqOfPlanes(const double* X, int m, int d, const double* X0, double* A_out, double* b_out) {
    if (!X || !X0 || !A_out || !b_out || m <= 0 || d <= 0) return 1;
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> Xmap(X, m, d);
    Eigen::Map<const Eigen::RowVectorXd> X0map(X0, d);

    MatrixXd A; VectorXd b; EqOfPlanes(Xmap, X0map, A, b);

    // row-major write-out
    for (int i = 0; i < m; i++) for (int j = 0; j < d; j++) A_out[(size_t)i * (size_t)d + (size_t)j] = A(i, j);
    for (int i = 0; i < m; i++) b_out[i] = b(i);
    return 0;
}

// Find_Direction: Pt ((d+1) x d) -> x (d)
int VA_Find_Direction(const double* Pt, int d, double* x_out) {
    if (!Pt || !x_out || d <= 0) return 1;
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> P(Pt, d + 1, d);
    const VectorXd x = Find_Direction(P);
    for (int j = 0; j < d; ++j) x_out[j] = x(j);
    return 0;
}

// CircumAll: Pt ((d+1) x d) -> center (d), R
int VA_CircumAll(const double* Pt, int d, double* center_out, double* R_out) {
    if (!Pt || !center_out || !R_out || d <= 0) return 1;
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> P(Pt, d + 1, d);
    VectorXd c; double R = 0.0; CircumAll(P, c, R);
    for (int j = 0; j < d; ++j) center_out[j] = c(j);
    *R_out = R; return 0;
}

// CircumPart: Pt (d x d) -> c (d), r
int VA_CircumPart(const double* Pt, int d, double* c_out, double* r_out) {
    if (!Pt || !c_out || !r_out || d <= 0) return 1;
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> P(Pt, d, d);
    VectorXd c; double r = 0.0; CircumPart(P, c, r);
    for (int j = 0; j < d; ++j) c_out[j] = c(j);
    *r_out = r; return 0;
}

// is_inside_simplex_bary: point (d), simplex ((d+1) x d) -> 1/0/-1
int VA_IsInsideSimplex_Bary(const double* point, const double* simplex, int d) {
    if (!point || !simplex || d <= 0) return -1;
    Eigen::Map<const VectorXd> p(point, d);
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> S(simplex, d + 1, d);
    return is_inside_simplex_bary(p, S);
}

} // extern "C"
#endif // BUILD_LIB && !BUILD_EXECUTABLE

// ------------------------------ CLI main --------------------------------
#if defined(BUILD_EXECUTABLE)

// Load points from a whitespace-separated text file (one point per line).
// If d>=0, verify each row has exactly d columns.
static bool load_points_txt(const char* path, int d, std::vector<double>& P, int& rows_out) {
    std::ifstream fin(path);
    if (!fin) return false;
    P.clear(); rows_out = 0; int cols = -1;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::vector<double> row; double x;
        while (iss >> x) row.push_back(x);
        if (cols < 0) cols = (int)row.size();
        if ((int)row.size() != cols) { std::cerr << "Inconsistent row width\n"; return false; }
        P.insert(P.end(), row.begin(), row.end());
        ++rows_out;
    }
    if (d >= 0 && cols != d) {
        std::cerr << "File has d=" << cols << " but provided d=" << d << "\n"; return false;
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <points.txt> <indexes.txt> <d> <order:min|max> [max_iter] [log_path]\n";
        return 1;
    }
    const char* points_path = argv[1];
    const char* idx_path    = argv[2];
    const int d = std::stoi(argv[3]);
    const std::string order = argv[4];
    if (order != "min" && order != "max") { std::cerr << "order must be 'min' or 'max'\n"; return 1; }

    // Load points
    std::vector<double> P; int rows = 0;
    if (!load_points_txt(points_path, d, P, rows)) { std::cerr << "Cannot open/parse " << points_path << "\n"; return 1; }

    // Load initial indices (d+1 integers)
    std::vector<int> Aa(d + 1), idx_loaded;
    if (!read_indices_file(idx_path, d + 1, idx_loaded) || (int)idx_loaded.size() != d + 1) {
        std::cerr << "Failed to read indices file or wrong count: " << idx_path << "\n"; return 1;
    }
    for (int i = 0; i < d + 1; ++i) Aa[i] = idx_loaded[i];

    int max_iter = 100000; const char* log_path = nullptr;
    if (argc >= 6) max_iter = std::stoi(argv[5]);
    if (argc >= 7) log_path = argv[6];

    std::vector<int> Aa_final(d + 1);
    std::vector<double> p_last(d, 0.0);
    double r_last = 0.0, L = 0.0; int iters = 0, Ksum = 0, k_last = 0, hit = 0;

    std::vector<int> pairs(2 * std::max(0, max_iter), -1);
    int num_pairs = 0;

    const int rc = VA_RunLoop_impl(
        P.data(), (unsigned long long)rows, d,
        Aa.data(), order.c_str(), max_iter, log_path,
        Aa_final.data(), p_last.data(), &r_last, &iters, &L, &Ksum, &k_last, &hit,
        &num_pairs, pairs.data());

    if (rc != 0) { std::cerr << "VA_RunLoop failed rc=" << rc << "\n"; return 1; }

    std::cout << "Aa_init:"; for (int v : Aa) std::cout << ' ' << v; std::cout << "\n";
    std::cout << "Aa_final:"; for (int v : Aa_final) std::cout << ' ' << v; std::cout << "\n";
    std::cout << "p_last:"; for (double v : p_last) std::cout << ' ' << v; std::cout << "\n";
    std::cout << "r_last: " << r_last
              << "\nphi_last: " << phiJ(d, r_last)
              << "\niters: " << iters
              << "\npath_length: " << L
              << "\nK_sum: " << Ksum
              << "\nk_last: " << k_last
              << "\nhit_inside: " << hit << "\n";

    std::cout << "pairs:";
    for (int s = 0; s < num_pairs; ++s) {
        const int left = pairs[2 * s + 0];
        const int neu  = pairs[2 * s + 1];
        std::cout << ' ' << left << "->" << neu;
    }
    std::cout << "\n";
    return 0;
}

#endif // BUILD_EXECUTABLE
