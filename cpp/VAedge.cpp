// VAedge.cpp — Eigen-based rewrite of EqOfPlanes / Find_Direction / CircumAll / CircumPart / FindNewVertex
// and loop driver with inside-simplex test via barycentric coordinates.
//
// Build:
//   Shared lib: g++ -O3 -std=c++17 -march=native -DNDEBUG -fPIC -shared VAedge.cpp -o libVAedge.so -DBUILD_LIB -ffast-math -I/path/to/eigen3 -I/path/to/nanoflann
//   Executable: g++ -O3 -std=c++17 -march=native -DNDEBUG VAedge.cpp -o VAedge -DBUILD_EXECUTABLE -ffast-math -I/path/to/eigen3 -I/path/to/nanoflann

#include <Eigen/Dense>
#include <nanoflann.hpp>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::ColPivHouseholderQR;

#include <vector>
#include <string>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iostream>
#include <sstream>

// ---- nanoflann compatibility helpers (global scope) ----
template <typename T>
inline auto nlfi_idx(const T& m) -> decltype(m.index) { return m.index; }
template <typename T>
inline auto nlfi_idx(const T& m) -> decltype(m.first) { return m.first; }

template <typename T>
inline auto nlfi_d2(const T& m) -> decltype(m.dist_sq) { return m.dist_sq; }
template <typename T>
inline auto nlfi_d2(const T& m) -> decltype(m.second) { return m.second; }

// Pick the ResultItem type that matches your nanoflann version
template <typename KD>
struct NlfiResultItem {
#if defined(NANOFLANN_VERSION) && (NANOFLANN_VERSION >= 0x0150)
    using type = nanoflann::ResultItem<typename KD::IndexType, typename KD::DistanceType>;
#else
    using type = nanoflann::ResultItem<typename KD::IndexType>;
#endif
};


// ---------- Small utilities ----------
#include <nanoflann.hpp>

// Row-major adaptor
struct RowMajorPoints {
    const double* pts = nullptr; size_t N = 0; int d = 0;
    inline size_t kdtree_get_point_count() const { return N; }
    inline double kdtree_get_pt(size_t idx, size_t dim) const {
        return pts[idx * (size_t)d + (size_t)dim];
    }
    template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
};

// KD-tree typedef
using KDAdaptor = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, RowMajorPoints>,
    RowMajorPoints, -1>;

static std::vector<int> stable_argsort_ascending(const std::vector<double>& a)
{
    std::vector<int> idx(a.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(),
        [&](int i, int j){ return a[i] < a[j]; });
    return idx;
}

static bool read_indices_file(const std::string& path, int expected_len, std::vector<int>& idx) {
    std::ifstream fin(path);
    if (!fin) return false;
    idx.clear();
    std::string tok;
    // optional header
    std::streampos pos = fin.tellg();
    if (fin >> tok) {
        if (tok.rfind("#", 0) != 0) {
            // first token is data, rewind to start of token
            fin.clear();
            fin.seekg(pos);
        } else {
            // consume rest of header line
            std::string rest; std::getline(fin, rest);
        }
    } else {
        return false;
    }
    int v;
    while (fin >> v) idx.push_back(v);
    return (expected_len <= 0) || ((int)idx.size() == expected_len);
}

static inline double sqr(double x){ return x*x; }
static inline double Vball(int d){ return std::pow(M_PI, 0.5*d) / std::tgamma(0.5*d + 1.0); }
static inline double phiJ(int d, double x){ return (d>0) ? Vball(d)*std::pow(x,d)/double(d) : 0.0; }
static inline double hmax(int d, double M, double rin=0.0){ return std::pow(M/Vball(d) + std::pow(rin,d), 1.0/d); }
static inline double Rin(int d, double Phi){ return std::pow(double(d)*Phi / Vball(d), 1.0/d); }

static inline bool in_set(int val, const int* arr, int n){
    for(int i=0;i<n;++i) if(arr[i]==val) return true;
    return false;
}

// Gather rows Points[Aa] into (len x d) matrix
static inline void gather_points(const double* Points, int d, const int* Aa, int len, MatrixXd& M){
    M.resize(len, d);
    for(int i=0;i<len;++i){
        const double* src = Points + (size_t)Aa[i]*(size_t)d;
        for(int j=0;j<d;++j) M(i,j) = src[j];
    }
}

// ---------- EqOfPlanes / Find_Direction / CircumAll / CircumPart ----------
static inline void EqOfPlanes(const MatrixXd& X, const RowVectorXd& X0, MatrixXd& A, VectorXd& b){
    // Python:
    // v = (X - X0)/2; xc = (X + X0)/2; b_i = sum(v_i * xc_i)
    const int m = (int)X.rows();
    const int d = (int)X.cols();
    A = 0.5*(X.rowwise() - X0);
    MatrixXd xc = 0.5*(X.rowwise() + X0);
    b.resize(m);
    for(int i=0;i<m;++i) b(i) = A.row(i).dot(xc.row(i));
}

static inline VectorXd Find_Direction(const MatrixXd& Pt){
    // Pt: (d+1) x d
    const int d = (int)Pt.cols();
    MatrixXd X = Pt.bottomRows(d);       // Pt[1:]
    RowVectorXd X0 = Pt.row(0);          // Pt[0]
    MatrixXd A; VectorXd b;
    EqOfPlanes(X, X0, A, b);
    // Least squares A x ≈ b
    ColPivHouseholderQR<MatrixXd> qr(A);
    return qr.solve(b);
}

static inline void CircumAll(const MatrixXd& Pt, VectorXd& center, double& R){
    center = Find_Direction(Pt);
    VectorXd diff = center - Pt.row(0).transpose();
    R = diff.norm();
}

// Pt: (d x d) here (because one vertex was removed before calling)
// Returns: c (center) and r (norm) via references.
static inline void CircumPart(const Eigen::MatrixXd& Pt, Eigen::VectorXd& c, double& r)
{
    const int d    = (int)Pt.cols();
    const int rows = (int)Pt.rows();      // expected: rows == d

    // x = Find_Direction(Pt)
    Eigen::VectorXd x = Find_Direction(Pt);

    // Base = Pt[1:] - Pt[0]          // shape: (d-1) x d
    Eigen::MatrixXd Base(rows - 1, d);
    for (int i = 0; i < rows - 1; ++i)
        Base.row(i) = Pt.row(i + 1) - Pt.row(0);

    // Q, R = np.linalg.qr(Base.T)    // Base.T: d x (d-1), default NumPy QR is reduced
    Eigen::MatrixXd BaseT = Base.transpose();              // d x (d-1)
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(BaseT);
    Eigen::MatrixXd Q = qr.householderQ() *
                        Eigen::MatrixXd::Identity(d, rows - 1); // Q: d x (d-1)  (thin)

    // v = Q @ Q.T @ (x - Pt[0])
    Eigen::VectorXd w = x - Pt.row(0).transpose();
    Eigen::VectorXd v = Q * (Q.transpose() * w);

    // return v + Pt[0], ||v||
    c = v + Pt.row(0).transpose();
    r = v.norm();
}


// ---------- is_inside_simplex (barycentric) ----------
// Returns: 1 inside, 0 outside, -1 singular
static inline int is_inside_simplex_bary(const VectorXd& point, const MatrixXd& simplex){
    const int d = (int)point.size();
    if (simplex.rows() != d+1 || simplex.cols()!=d) return -1;

    // Build augmented system:
    // A ( (d+1)x(d+1) ) = [ simplex^T ; 1 ... 1 ]
    // rhs = [ point ; 1 ]
    MatrixXd A(d+1, d+1);
    A.topRows(d) = simplex.transpose();
    A.row(d).setOnes();
    VectorXd rhs(d+1); rhs.head(d) = point; rhs(d) = 1.0;

    Eigen::FullPivLU<MatrixXd> lu(A);
    if (!lu.isInvertible()) return -1;
    VectorXd bary = lu.solve(rhs);

    for(int j=0;j<d+1;++j){
        const double y = bary(j);
        if (y < 0.0 || y > 1.0) return 0;
    }
    return 1;
}

// ---------- FindNewVertex ----------


// KD-aware core: reuse KD built elsewhere (used inside VA_RunLoop)
static inline int FindNewVertex_core(
    const double* Points, unsigned long long N, int d,
    const int* Aa, int Aa_len, const std::string& order,
    KDAdaptor& kdt,
    std::vector<int>& out_Aa_new,
    Eigen::VectorXd& NV, double& NR, int& k_out);

// KD-less shim: builds a temporary KD, forwards to KD-aware core (used by C API)
static inline int FindNewVertex_core_shim(
    const double* Points, unsigned long long N, int d,
    const int* Aa, int Aa_len, const std::string& order,
    std::vector<int>& out_Aa_new,
    Eigen::VectorXd& NV, double& NR, int& k_out)
{
    RowMajorPoints cloud{Points, static_cast<size_t>(N), d};
    KDAdaptor kdt(d, cloud, {10});
    kdt.buildIndex();
    return FindNewVertex_core(Points, N, d, Aa, Aa_len, order,
                              kdt, out_Aa_new, NV, NR, k_out);
}

static inline int FindNewVertex_core(
    const double* Points, unsigned long long N, int d,
    const int* Aa, int Aa_len, const std::string& order,
    KDAdaptor& kdt,                                  // <-- KD passed in
    std::vector<int>& out_Aa_new,
    Eigen::VectorXd& NV, double& NR, int& k_out)
{
    // Assemble PtSel = Points[Aa]
    Eigen::MatrixXd PtSel(Aa_len, d);
    for (int i = 0; i < Aa_len; ++i)
        for (int j = 0; j < d; ++j)
            PtSel(i,j) = Points[(size_t)Aa[i]*(size_t)d + (size_t)j];

    // V, R for current Aa
    Eigen::VectorXd V; double R = 0.0;
    CircumAll(PtSel, V, R);
    const int DIM = Aa_len - 1;

    k_out = 0;
    std::vector<int> I;
    std::vector<double> Rat;

    // Build Rat, I
    for (int i = 0; i < Aa_len; ++i) {
        std::vector<int> Am; Am.reserve(Aa_len - 1);
        for (int t = 0; t < Aa_len; ++t) if (t != i) Am.push_back(Aa[t]);

        Eigen::MatrixXd PtAm(DIM, d);
        for (int r = 0; r < DIM; ++r)
            for (int j = 0; j < d; ++j)
                PtAm(r,j) = Points[(size_t)Am[r]*(size_t)d + (size_t)j];

        Eigen::VectorXd c; double r = 0.0;
        CircumPart(PtAm, c, r);

        const double cV = (V - c).norm();
        if (cV == 0.0) continue;
        const double rat = r / cV;

        Eigen::VectorXd VmC = V - c;
        Eigen::VectorXd Pi  = Eigen::Map<const Eigen::VectorXd>(
                                  Points + (size_t)Aa[i]*(size_t)d, d) - V;
        Eigen::VectorXd Pj0 = Eigen::Map<const Eigen::VectorXd>(
                                  Points + (size_t)Am[0]*(size_t)d, d) - V;

        if (VmC.dot(Pi) < VmC.dot(Pj0)) {
            Rat.push_back(rat);
            I.push_back(i);
            ++k_out;
        }
    }
    if (I.empty()) return 2;

    // Stable argsort (Python behavior)
    auto ord = stable_argsort_ascending(Rat);
    int chosen_pos = (order == "min") ? ord.back() : ord.front();
    const int i_max = I[chosen_pos];

    // Face for i_max
    std::vector<int> Amx; Amx.reserve(Aa_len - 1);
    for (int t = 0; t < Aa_len; ++t) if (t != i_max) Amx.push_back(Aa[t]);

    Eigen::MatrixXd PtAmx(DIM, d);
    for (int r = 0; r < DIM; ++r)
        for (int j = 0; j < d; ++j)
            PtAmx(r,j) = Points[(size_t)Amx[r]*(size_t)d + (size_t)j];

    Eigen::VectorXd cmax; double rmax = 0.0;
    CircumPart(PtAmx, cmax, rmax);

    const double Vc = (V - cmax).norm();
    if (Vc == 0.0) return 3;
    Eigen::VectorXd vdir = (V - cmax) / Vc;

    // Radius expansion
    double rr = (DIM != 2) ? std::pow(4.0 * DIM, 1.0 / double(DIM)) : 10.0;
    double a_min = 1000.0;
    int    pt_min = -1;

    while (a_min == 1000.0 && rr < 101.0) {
        const double rad2 = (R * rr) * (R * rr);

        using RI = typename NlfiResultItem<KDAdaptor>::type;
        std::vector<RI> matches;

        #if defined(NANOFLANN_HAS_SEARCH_PARAMS)
            nanoflann::SearchParameters params; params.sorted = false;
            kdt.radiusSearch(V.data(),
                             (KDAdaptor::DistanceType)rad2,
                             matches, params);
        #else
            kdt.radiusSearch(V.data(),
                             (KDAdaptor::DistanceType)rad2,
                             matches);
        #endif

        // sort by index to mimic Python's pt=0..N-1 order
        std::sort(matches.begin(), matches.end(),
                  [](const RI& a, const RI& b){ return nlfi_idx(a) < nlfi_idx(b); });

        // evaluate candidates
        for (const auto& m : matches) {
            const int pi = (int)nlfi_idx(m);

            bool inAa = false;
            for (int t = 0; t < Aa_len; ++t) if (Aa[t] == pi) { inAa = true; break; }
            if (inAa) continue;

            Eigen::VectorXd Pv = Eigen::Map<const Eigen::VectorXd>(
                                     Points + (size_t)pi*(size_t)d, d) - V;

            // inclusive radius, like query_ball_point
            if (Pv.squaredNorm() > rad2) continue;

            const double Pj = vdir.dot(Pv);
            if (Pj == 0.0) continue;

            const double rd = (Pv - Pj*vdir).norm();

            const double num = rd*rd + Pj*Pj - rmax*rmax - Vc*Vc;
            const double den = 2.0*(Pj + Vc);
            if (den == 0.0) continue;

            const double a = num / den;
            if (a > 0.0 && a < a_min) {
                a_min = a;
                pt_min = pi;
            }
        }
        rr *= 1.1;
    }

    if (rr > 100.0) k_out = 9999;
    if (pt_min < 0) return 4;

    // New Aa (replace i_max)
    out_Aa_new.assign(Aa, Aa + Aa_len);
    out_Aa_new[i_max] = pt_min;

    // NV, NR on new set
    Eigen::MatrixXd PtNew(Aa_len, d);
    for (int i = 0; i < Aa_len; ++i)
        for (int j = 0; j < d; ++j)
            PtNew(i,j) = Points[(size_t)out_Aa_new[i]*(size_t)d + (size_t)j];

    CircumAll(PtNew, NV, NR);
    return 0;
}



// ---- C API (only when building the shared library) ----
#if defined(BUILD_LIB) && !defined(BUILD_EXECUTABLE)

// ================== COMPONENT TEST EXPORTS (append to VAedge.cpp) ==================
extern "C" {

int VA_FindNewVertex(
    const double* Points, unsigned long long N, int d,
    const int* Aa, int Aa_len,
    const char* order,
    int* out_Aa_new,
    double* out_NV,
    double* out_NR,
    int* out_k)
{
    if (!Points || !Aa || !out_Aa_new || !out_NV || !out_NR || !out_k) return 5;

    std::vector<int> newAa;
    Eigen::VectorXd NV(d); double NR = 0.0; int k = 0;

    // Call the KD-less shim (builds a temporary KD and forwards to the KD-aware core)
    int rc = FindNewVertex_core_shim(
        Points, N, d, Aa, Aa_len,
        order ? std::string(order) : std::string("max"),
        newAa, NV, NR, k);
    if (rc != 0) return rc;

    for (int i = 0; i < Aa_len; ++i) out_Aa_new[i] = newAa[i];
    for (int j = 0; j < d;   ++j) out_NV[j] = NV(j);
    *out_NR = NR; *out_k = k;
    return 0;
}

// EqOfPlanes: X (m x d), X0 (d) -> A (m x d), b (m)
int VA_EqOfPlanes(const double* X, int m, int d,
                  const double* X0,
                  double* A_out, double* b_out)
{
    if (!X || !X0 || !A_out || !b_out || m<=0 || d<=0) return 1;
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> Xmap(X, m, d);
    Eigen::Map<const Eigen::RowVectorXd> X0map(X0, d);

    Eigen::MatrixXd A; Eigen::VectorXd b;
    EqOfPlanes(Xmap, X0map, A, b);

    // row-major write
    for (int i=0;i<m;i++) for (int j=0;j<d;j++) A_out[(size_t)i*(size_t)d + (size_t)j] = A(i,j);
    for (int i=0;i<m;i++) b_out[i] = b(i);
    return 0;
}

// Find_Direction: Pt ((d+1) x d) -> x (d)
int VA_Find_Direction(const double* Pt, int d, double* x_out)
{
    if (!Pt || !x_out || d <= 0) return 1;
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> P(Pt, d+1, d);
    Eigen::VectorXd x = Find_Direction(P);
    for (int j = 0; j < d; ++j) x_out[j] = x(j);
    return 0;
}

// CircumAll: Pt ((d+1) x d) -> center (d), R
int VA_CircumAll(const double* Pt, int d, double* center_out, double* R_out)
{
    if (!Pt || !center_out || !R_out || d<=0) return 1;
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> P(Pt, d+1, d);
    Eigen::VectorXd c; double R=0.0;
    CircumAll(P, c, R);
    for (int j=0;j<d; ++j) center_out[j] = c(j);
    *R_out = R;
    return 0;
}

// CircumPart: Pt ((d+1) x d) -> c (d), r
int VA_CircumPart(const double* Pt, int d, double* c_out, double* r_out)
{
    if (!Pt || !c_out || !r_out || d<=0) return 1;
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> P(Pt, d, d);
    Eigen::VectorXd c; double r=0.0;
    CircumPart(P, c, r);
    for (int j=0;j<d; ++j) c_out[j] = c(j);
    *r_out = r;
    return 0;
}

// is_inside_simplex_bary: point (d), simplex ((d+1) x d) -> 1/0/-1
int VA_IsInsideSimplex_Bary(const double* point, const double* simplex, int d)
{
    if (!point || !simplex || d<=0) return -1;
    Eigen::Map<const Eigen::VectorXd> p(point, d);
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> S(simplex, d+1, d);
    return is_inside_simplex_bary(p, S);
}


// VA_FindNewVertex is already exported in your build; keep it as is.

} // extern "C"
#endif

// ---------- C API: VA_RunLoop ----------
extern "C" int VA_RunLoop(
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
    int*  out_hit_inside
){
    if (!Points || !Aa_init || !out_Aa_final || !out_p_last || !out_r_last ||
        !out_iters || !out_path_length || !out_K_sum || !out_k_last || !out_hit_inside) return 5;

    const int Aa_len = d + 1;
    std::vector<int> Aa(Aa_len);
    for (int i=0;i<Aa_len;++i) Aa[i] = Aa_init[i];

    // initial circumcenter
    Eigen::MatrixXd PtSel(Aa_len, d);
    for (int i=0;i<Aa_len;++i)
        for (int j=0;j<d;++j)
            PtSel(i,j)=Points[(size_t)Aa[i]*(size_t)d+(size_t)j];
    Eigen::VectorXd p_curr(d); double r_curr=0.0;
    CircumAll(PtSel, p_curr, r_curr);

    // Build KD-tree once
    RowMajorPoints cloud{Points, (size_t)N, d};
    KDAdaptor kdt(d, cloud, {10});
    kdt.buildIndex();

    int i = 0, Ksum = 0, k_last = 0;
    double L = 0.0;
    const int DIM = Aa_len - 1;

    std::ofstream flog;
    if (log_path && std::strlen(log_path)>0) flog.open(log_path);

    Eigen::MatrixXd simplex(Aa_len, d);
    *out_hit_inside = 0;

    while (true) {
        // build simplex and check inside
        for (int t=0;t<Aa_len;++t)
            for (int j=0;j<d;++j)
                simplex(t,j)=Points[(size_t)Aa[t]*(size_t)d+(size_t)j];

        int inside = is_inside_simplex_bary(p_curr, simplex);
        if (inside == 1) { *out_hit_inside = 1; break; }
        if (i >= max_iter) break;
        if (k_last >= 9999) break;

        // step
        std::vector<int> Aa_new(Aa_len);
        Eigen::VectorXd NV(d); double NR=0.0; int k=0;
        int rc = FindNewVertex_core(Points, N, d,
                                    Aa.data(), (int)Aa.size(),
                                    order ? order : "max",
                                    kdt,                    // <-- pass KD
                                    Aa_new, NV, NR, k);
        if (rc != 0) { k_last = k; break; }

        if (k>1 && flog.good()){
            double dist = (NV - p_curr).norm();
            flog << i << ' ' << dist << ' ' << NR << ' ' << phiJ(DIM, NR) << ' '
                 << L << ' ' << k << ' ' << p_curr.norm() << '\n';
        }

        L += (NV - p_curr).norm();
        ++i; Ksum += k; k_last = k;
        Aa.swap(Aa_new);
        p_curr = NV;
        r_curr = NR;
    }

    for (int t=0;t<Aa_len;++t) out_Aa_final[t] = Aa[t];
    for (int j=0;j<d;++j) out_p_last[j] = p_curr(j);
    *out_r_last = r_curr;
    *out_iters  = i;
    *out_path_length = L;
    *out_K_sum  = Ksum;
    *out_k_last = k_last;
    return 0;
}

#if defined(BUILD_EXECUTABLE)
// CLI:
//   ./VAedge <points.txt> <d> <order:min|max> [max_iter] [log_path]
// It will read indices from "<points.txt>.idx" (created by initial_proj executable).
int main(int argc, char** argv){
    if (argc < 4){
        std::cerr << "Usage: " << argv[0] << " <points.txt> <d> <order:min|max> [max_iter] [log_path]\n";
        return 1;
    }
    const char* path = argv[1];
    const int d = std::stoi(argv[2]);
    const std::string order = argv[3];
    if (order!="min" && order!="max"){
        std::cerr << "order must be 'min' or 'max'\n"; return 1;
    }

    // read points
    std::ifstream fin(path);
    if (!fin){ std::cerr << "Cannot open " << path << "\n"; return 1; }
    std::vector<double> P; P.reserve(1<<20);
    std::string line;
    int cols=-1, rows=0;
    while (std::getline(fin, line)){
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::vector<double> row; double x;
        while (iss >> x) row.push_back(x);
        if (cols<0) cols = (int)row.size();
        if ((int)row.size()!=cols){ std::cerr << "Inconsistent row width\n"; return 1; }
        P.insert(P.end(), row.begin(), row.end());
        ++rows;
    }
    if (cols!=d){ std::cerr << "File has d="<<cols<<" but provided d="<<d<<"\n"; return 1; }

    // read indices sidecar
    std::vector<int> Aa(d+1);
    const std::string idx_path = std::string(path) + ".idx";
    std::vector<int> idx_loaded;
    if (!read_indices_file(idx_path, d+1, idx_loaded)) {
        std::cerr << "Failed to read indices file: " << idx_path
                  << " (run initial_proj executable first).\n";
        return 1;
    }
    if ((int)idx_loaded.size() != d+1) {
        std::cerr << "Indices file does not contain d+1 indices.\n"; return 1;
    }
    for (int i=0;i<d+1;++i) Aa[i] = idx_loaded[i];

    int max_iter = 100000;
    const char* log_path = nullptr;
    if (argc >= 5) max_iter = std::stoi(argv[4]);
    if (argc >= 6) log_path = argv[5];

    std::vector<int> Aa_final(d+1);
    std::vector<double> p_last(d, 0.0);
    double r_last=0.0, L=0.0;
    int iters=0, Ksum=0, k_last=0, hit=0;

    int rc = VA_RunLoop(
        P.data(), (unsigned long long)rows, d,
        Aa.data(), order.c_str(), max_iter, log_path,
        Aa_final.data(), p_last.data(), &r_last, &iters, &L, &Ksum, &k_last, &hit
    );
    if (rc!=0){ std::cerr << "VA_RunLoop failed rc="<<rc<<"\n"; return 1; }

    std::cout << "Aa_init:";
    for (int v : Aa) std::cout << ' ' << v;
    std::cout << "\nAa_final:";
    for (int v : Aa_final) std::cout << ' ' << v;
    std::cout << "\np_last:";
    for (double v : p_last) std::cout << ' ' << v;
    std::cout << "\nr_last: " << r_last
              << "\niters: " << iters
              << "\npath_length: " << L
              << "\nK_sum: " << Ksum
              << "\nk_last: " << k_last
              << "\nhit_inside: " << hit << "\n";
    return 0;
}
#endif
