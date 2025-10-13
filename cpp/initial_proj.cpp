// initial_proj.cpp — Eigen-based implementation of InitialProjection
// Build:
//   Shared lib: g++ -O3 -std=c++17 -shared -fPIC initial_proj.cpp -o libinitial_proj.so -DBUILD_LIB -I/path/to/eigen3
//   Executable: g++ -O3 -std=c++17 initial_proj.cpp -o initial_proj -DBUILD_EXECUTABLE -I/path/to/eigen3

#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <vector>
#include <cmath>
#include <cstdint>
#include <limits>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

// write one index per line, with a header
static bool write_indices_file(const std::string& path, const std::vector<int>& idx) {
    std::ofstream fout(path);
    if (!fout) return false;
    fout << "# d_plus_1 " << idx.size() << "\n";
    for (int v : idx) fout << v << "\n";
    return true;
}

// -------- helpers --------
static inline double vec_norm(const double* x, int d) {
    double s = 0.0; for (int i = 0; i < d; ++i) s += x[i]*x[i]; return std::sqrt(s);
}

// Build A0 (m x d) and b0 (m) from Points (N x d), with m = N-1:
//   A0[i,:] = normalize(Points[i+1,:])
//   b0[i]   = 0.5 * ||Points[i+1,:]||
static void build_A0_b0_from_points(const double* Points, unsigned long long N, int d,
                                    MatrixXd& A0, VectorXd& b0)
{
    const unsigned long long m = N - 1ULL;
    A0.resize((int)m, d);
    b0.resize((int)m);
    for (unsigned long long i = 0; i < m; ++i) {
        const double* p = Points + ((size_t)i + 1U) * (size_t)d;
        double nrm = 0.0; for (int j = 0; j < d; ++j) nrm += p[j]*p[j];
        nrm = std::sqrt(nrm);
        b0((int)i) = 0.5 * nrm;
        if (nrm > 0.0) {
            for (int j = 0; j < d; ++j) A0((int)i, j) = p[j] / nrm;
        } else {
            for (int j = 0; j < d; ++j) A0((int)i, j) = 0.0;
        }
    }
}

// One InitialProjection step (Eigen version).
// Inputs:
//   A0: (m x d), b0: (m), Aa: selected face indices (size k), x0: (d), dir: (d)
// Outputs:
//   x_out: stepped point (d), v_out: new direction (d), picked_idx: chosen face
// Returns 0 on success, nonzero if no valid face found.
static int initial_projection_step_eigen(const MatrixXd& A0, const VectorXd& b0,
                                         const std::vector<int>& Aa,
                                         const VectorXd& x0, const VectorXd& dir,
                                         VectorXd& x_out, VectorXd& v_out, int& picked_idx)
{
    const int m = (int)A0.rows();
    const int d = (int)A0.cols();
    const double INF = std::numeric_limits<double>::infinity();
    const double TOL = 1e-12;

    // Safe distances: D = (b - A0*x0) / (A0*dir), with mask before divide
    VectorXd den = A0 * dir;      // m
    VectorXd num = b0 - A0 * x0;  // m
    std::vector<double> D(m, INF);

    std::vector<char> used(m, 0);
    for (int idx : Aa) if (idx >= 0 && idx < m) used[idx] = 1;

    for (int i = 0; i < m; ++i) {
        if (used[i]) continue;
        double di = den(i);
        if (!(di > 0.0)) continue; // mask den <= 0
        double val = num(i) / di;
        if (val >= 0.0 && std::isfinite(val)) D[i] = val;
    }

    // argmin with deterministic tie-break (smallest index within tol)
    picked_idx = -1;
    double best = INF;
    for (int i = 0; i < m; ++i) {
        double val = D[i];
        if (!std::isfinite(val)) continue;
        if (picked_idx < 0 || val < best - TOL || (std::abs(val - best) <= TOL && i < picked_idx)) {
            best = val; picked_idx = i;
        }
    }
    if (picked_idx < 0) return 1;

    // Step: x = x0 + D[picked] * dir
    x_out = x0 + best * dir;

    // Build M = A0[Aa ∪ {picked}].T  (d x k)
    std::vector<int> idxs = Aa;
    idxs.push_back(picked_idx);
    const int k = (int)idxs.size();

    MatrixXd M(d, k);
    for (int c = 0; c < k; ++c) {
        M.col(c) = A0.row(idxs[c]).transpose();
    }

    // Householder QR with column pivoting (rank-aware), like NumPy’s reduced QR
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(M);
    const int rnk = qr.rank();
    Eigen::MatrixXd Qthin = qr.householderQ() * Eigen::MatrixXd::Identity(d, rnk);

    // v = dir - Q(Q^T dir)
    if (rnk > 0) {
        v_out = dir - Qthin * (Qthin.transpose() * dir);
    } else {
        v_out = dir;
    }

    // sign flip by sign(v·x0)
    double Dv = v_out.dot(x0);
    if (Dv < 0.0) v_out = -v_out;

    return 0;
}

// Select indices {0} ∪ {1 + Aa[i]} of size d+1 from Points (N x d)
static int select_indices_from_points_eigen(const double* Points, unsigned long long N, int d,
                                            std::vector<int>& out_idx)
{
    if (N < 2ULL || d <= 0) return 1;
    const unsigned long long m = N - 1ULL;

    // Build A0, b0
    MatrixXd A0; VectorXd b0;
    build_A0_b0_from_points(Points, N, d, A0, b0);

    // x0, dir
    VectorXd x0(d);
    for (int j = 0; j < d; ++j) x0(j) = -Points[j];
    double nx = x0.norm();
    VectorXd dir = (nx > 0.0 ? x0 / nx : x0);

    // iterate
    std::vector<int> Aa; Aa.reserve(d);
    VectorXd x(d), v(d);
    for (int s = 0; s < d; ++s) {
        int pick = -1;
        if (initial_projection_step_eigen(A0, b0, Aa, x0, dir, x, v, pick) != 0) return 2;
        Aa.push_back(pick);
        x0 = x; dir = v;
    }

    // map to original indices
    out_idx.resize(d + 1);
    out_idx[0] = 0;
    for (int i = 0; i < d; ++i) out_idx[i + 1] = 1 + Aa[i];
    return 0;
}

// --------- LIB API ---------
#if defined(BUILD_LIB)
extern "C" {

// Fill out_indices[0..d] with selected indices.
int proj_select_from_points(const double* Points, unsigned long long N, int d,
                            int* out_indices /*size d+1*/)
{
    if (!Points || !out_indices) return 3;
    std::vector<int> idx;
    int rc = select_indices_from_points_eigen(Points, N, d, idx);
    if (rc != 0) return rc;
    for (int i = 0; i < d + 1; ++i) out_indices[i] = idx[i];
    return 0;
}


} // extern "C"
#endif // BUILD_LIB

// --------- EXE main ---------
#if defined(BUILD_EXECUTABLE)
int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <points.txt> <d>\n";
        return 1;
    }
    const char* path = argv[1];
    const int d = std::stoi(argv[2]);

    std::ifstream fin(path);
    if (!fin) { std::cerr << "Cannot open: " << path << "\n"; return 1; }

    std::vector<double> P;
    P.reserve(1<<20);
    std::string line;
    int cols = -1;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::vector<double> row;
        double x;
        while (iss >> x) row.push_back(x);
        if (cols < 0) cols = (int)row.size();
        if ((int)row.size() != cols) { std::cerr << "Inconsistent row width\n"; return 1; }
        for (double v : row) P.push_back(v);
    }
    if (cols != d) { std::cerr << "File has d=" << cols << " but provided d=" << d << "\n"; return 1; }
    if (P.empty()) { std::cerr << "Empty input\n"; return 1; }

    const unsigned long long N = (unsigned long long)(P.size() / (size_t)d);
    if (N < (unsigned long long)(d + 1)) {
        std::cerr << "Need at least d+1 points\n"; return 1;
    }

    std::vector<int> idx;
    int rc = select_indices_from_points_eigen(P.data(), N, d, idx);
    if (rc != 0) { std::cerr << "Projection failed (rc=" << rc << ")\n"; return 1; }

    // print selected points (unchanged)
    for (int i = 0; i < d + 1; ++i) {
        const double* rowp = P.data() + (size_t)idx[i] * (size_t)d;
        for (int j = 0; j < d; ++j) {
            if (j) std::cout << ' ';
            std::cout.setf(std::ios::fmtflags(0), std::ios::floatfield);
            std::cout.precision(12);
            std::cout << rowp[j];
        }
        std::cout << "\n";
    }
    std::cout << "# indices:";
    for (int i = 0; i < d + 1; ++i) std::cout << ' ' << idx[i];
    std::cout << "\n";

    // --- write indices sidecar file ---
    const std::string idx_path = std::string(path) + ".idx";
    if (!write_indices_file(idx_path, idx)) {
        std::cerr << "Failed to write indices file: " << idx_path << "\n";
        return 1;
    }

    return 0;
}
#endif // BUILD_EXECUTABLE
