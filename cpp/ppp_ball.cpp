// ppp_ball.cpp  — build exe with -DBUILD_EXECUTABLE, library with -DBUILD_LIB
//
//g++ -O3 -std=c++17 ppp_ball.cpp -o ppp_ball -DBUILD_EXECUTABLE
//g++ -O3 -std=c++17 -shared -fPIC ppp_ball.cpp -o libppp_ball.so -DBUILD_LIB

#include <random>
#include <vector>
#include <cmath>
#include <cstdint>
#include <stdexcept>

static inline double unit_ball_volume(int d) {
    return std::pow(M_PI, 0.5 * d) / std::tgamma(0.5 * d + 1.0);
}

static inline double compute_radius_from_lambda_M(int d, double lambda, std::size_t M) {
    if (d <= 0 || lambda <= 0.0) throw std::invalid_argument("bad d or lambda");
    const double Cd = unit_ball_volume(d);
    const double base = static_cast<double>(M) / (lambda * Cd);
    return std::pow(base, 1.0 / static_cast<double>(d));
}

// Flat array output: size = M * d (row-major)
static inline void generate_points_flat(std::size_t M, int d, double R, uint64_t seed,
                                        std::vector<double>& out_flat)
{
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uni01(0.0, 1.0);
    std::normal_distribution<double> normal(0.0, 1.0);

    // 1) Exponential(1) spacings for Uniform order stats (Rényi representation)
    const std::size_t K = M + 1;
    std::vector<double> E(K);
    double S = 0.0;
    for (std::size_t j = 0; j < K; ++j) {
        double u;
        do { u = uni01(rng); } while (u <= 0.0);  // avoid log(0)
        E[j] = -std::log(u);                      // Exp(1)
        S += E[j];
    }

    // 2) Build sorted uniforms U_(i) via cumulative sums, convert to radii r_i
    out_flat.clear();
    out_flat.reserve(M * static_cast<std::size_t>(d));

    const double inv_d = 1.0 / static_cast<double>(d);
    double cum = 0.0;

    // Reuse direction buffer
    std::vector<double> dir(d);

    for (std::size_t i = 0; i < M; ++i) {
        cum += E[i];
        const double Ui = cum / S;                     // i-th order statistic
        const double ri = R * std::pow(Ui, inv_d);     // radius in ascending order

        // 3) Independent random direction (normalized Gaussian)
        double norm2 = 0.0;
        do {
            norm2 = 0.0;
            for (int j = 0; j < d; ++j) {
                double g = normal(rng);
                dir[j] = g;
                norm2 += g * g;
            }
        } while (norm2 == 0.0);
        const double inv = 1.0 / std::sqrt(norm2);

        for (int j = 0; j < d; ++j)
            out_flat.push_back(ri * (dir[j] * inv));
    }
}

#if defined(BUILD_LIB)
// --------------------- Shared library API ---------------------
#include <cstring>

extern "C" {

int ppp_generate(
    double lambda,                 // density
    unsigned long long M,          // number of points
    int d,                         // dimension
    unsigned long long seed,       // RNG seed
    double* out_points,            // size M*d (row-major)
    double* out_R                  // ball radius
);

// Compute radius from (lambda, M, d)
double ppp_compute_radius(int d, double lambda, unsigned long long M) {
    return compute_radius_from_lambda_M(d, lambda, static_cast<std::size_t>(M));
}

// Generate M points. `out` must have size M*d. Writes R to *out_R.
// Returns 0 on success, nonzero on error.
int ppp_generate(double lambda, unsigned long long M, int d, unsigned long long seed,
                 double* out, double* out_R)
{
    if (!out || !out_R || d <= 0 || M == 0 || lambda <= 0.0) return 1;
    const double R = compute_radius_from_lambda_M(d, lambda, static_cast<std::size_t>(M));
    std::vector<double> flat;
    generate_points_flat(static_cast<std::size_t>(M), d, R, static_cast<uint64_t>(seed), flat);
    std::memcpy(out, flat.data(), flat.size() * sizeof(double));
    *out_R = R;
    return 0;
}

} // extern "C"
#endif // BUILD_LIB

#if defined(BUILD_EXECUTABLE)
// --------------------- Executable main ---------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <lambda> <M> <d> <seed>\n";
        return 1;
    }
    const double lambda = std::stod(argv[1]);
    const std::size_t M = static_cast<std::size_t>(std::stoull(argv[2]));
    const int d = std::stoi(argv[3]);
    const uint64_t seed = std::stoull(argv[4]);

    if (lambda <= 0.0 || d <= 0 || M == 0) {
        std::cerr << "lambda>0, d>=1, M>=1 required\n";
        return 1;
    }

    const double R = compute_radius_from_lambda_M(d, lambda, M);

    std::vector<double> flat;
    generate_points_flat(M, d, R, seed, flat);

    std::ostringstream rstr;
    rstr << std::setprecision(8) << std::scientific << R;

    std::ostringstream fname;
    fname << "points_d" << d
          << "_lambda" << std::setprecision(6) << std::fixed << lambda
          << "_M" << M
          << "_R" << rstr.str()
          << "_seed" << seed << ".dat";
    const std::string filename = fname.str();

    std::ofstream fout(filename);
    fout << std::setprecision(12);
    for (std::size_t i = 0; i < M; ++i) {
        for (int j = 0; j < d; ++j) {
            fout << flat[i * static_cast<std::size_t>(d) + j];
            if (j < d - 1) fout << ' ';
        }
        fout << '\n';
    }
    fout.close();

    std::cout << filename << std::endl;
    return 0;
}
#endif // BUILD_EXECUTABLE
