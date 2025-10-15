// ppp_ball.cpp : 
// point Poisson process ball sampling of M points in d dimensions with density lambda ordered according with the distance from the origin
// CLI an API
// Build executable:   g++ -O3 -std=c++17 ppp_ball.cpp -o ppp_ball -DBUILD_EXECUTABLE
// Build shared lib:   g++ -O3 -std=c++17 -shared -fPIC ppp_ball.cpp -o libppp_ball.so -DBUILD_LIB

#include <random>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <cstring>
#include <iomanip>
#include <limits>


// --------- math helpers ---------
static inline double unit_ball_volume(int d) {
    if (d <= 0) throw std::invalid_argument("dimension must be >= 1");
    const double pi = std::acos(-1.0);
    return std::pow(pi, 0.5 * d) / std::tgamma(0.5 * d + 1.0);
}

static inline double compute_radius_from_lambda_M(int d, double lambda, int M) {
    if (d <= 0 || lambda <= 0.0 || M <= 0) throw std::invalid_argument("bad d, lambda, or M");
    const double Cd = unit_ball_volume(d);
    const double base = static_cast<double>(M) / (lambda * Cd);
    return std::pow(base, 1.0 / static_cast<double>(d));
}

// Generate M points uniformly in d-ball of radius R; flat row-major M*d
static inline void generate_points_flat(int d, int M, double R, unsigned int seed,
                                        std::vector<double>& out_flat)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uni01(0.0, 1.0);
    std::normal_distribution<double> normal(0.0, 1.0);

    // Rényi: exponential spacings to generate sorted uniforms
    const int K = M + 1;
    std::vector<double> E(K);
    double S = 0.0;
    for (int j = 0; j < K; ++j) {
        double u;
        do { u = uni01(rng); } while (u <= 0.0); // guard log(0)
        E[j] = -std::log(u);
        S += E[j];
    }

    out_flat.clear();
    out_flat.reserve(static_cast<size_t>(M) * static_cast<size_t>(d));

    const double inv_d = 1.0 / static_cast<double>(d);
    double cum = 0.0;
    std::vector<double> dir(d);

    for (int i = 0; i < M; ++i) {
        cum += E[i];
        const double Ui = cum / S;                 // i-th order statistic of Uniform(0,1)
        const double ri = R * std::pow(Ui, inv_d); // correct radius distribution

        // Random direction via normalized Gaussian
        double norm2;
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
extern "C" {

// Compute radius from (d, lambda, M)
double ppp_compute_radius(int d, double lambda, int M) {
    return compute_radius_from_lambda_M(d, lambda, M);
}

// Generate M points. `out_points` must have size M*d (row-major). Returns 0 on success.
int ppp_generate(int d, double lambda, int M, unsigned int seed, double* out_points) {
    if (!out_points || d <= 0 || M <= 0 || lambda <= 0.0) return 1;
    const double R = compute_radius_from_lambda_M(d, lambda, M);
    std::vector<double> flat;
    generate_points_flat(d, M, R, seed, flat);
    std::memcpy(out_points, flat.data(), flat.size() * sizeof(double));
    return 0;
}

} // extern "C"
#endif // BUILD_LIB

#if defined(BUILD_EXECUTABLE)
// --------------------- Executable main ---------------------
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <d> <lambda> <M> <seed>\n";
        return 1;
    }

    const int           d      = std::stoi(argv[1]);
    const double        lambda = std::stod(argv[2]);
    const int           M      = std::stoi(argv[3]);
    const unsigned int  seed   = static_cast<unsigned int>(std::stoul(argv[4]));

    if (lambda <= 0.0 || d <= 0 || M <= 0) {
        std::cerr << "lambda>0, d>=1, M>=1 required\n";
        return 1;
    }

    const double R = compute_radius_from_lambda_M(d, lambda, M);

    std::vector<double> samples;
    generate_points_flat(d, M, R, seed, samples);

    std::cout << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);

    // M lines, d columns.
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < d; ++j) {
            std::cout << samples[static_cast<size_t>(i) * static_cast<size_t>(d) + j];
            if (j < d - 1) std::cout << ' ';
        }
        std::cout << '\n';
    }

    return 0;
}
#endif // BUILD_EXECUTABLE
