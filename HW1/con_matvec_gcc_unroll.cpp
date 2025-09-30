
#include <iostream>
#include <iomanip>
#include <random>
#include <sys/time.h> // for gettimeofday
#include <stdlib.h>

// timer 
double wall_time() {
    timeval tv; gettimeofday(&tv, nullptr);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

// matrix calculation
void matrix_vect_1d(double* A, double* x, double* y, int Nrow, int Mcol) {
    for (int i = 0; i < Nrow; i++) {
        const double* a = A + i*Mcol;
        double s = 0.0;

       // #pragma GCC unroll 4
        for (int j=0; j<Mcol; ++j) {
            s += a[j] * x[j];
        }
        y[i] = s;
    }
}

void run_case(int Nrow, int Mcol) {
    // allocate contiguous memory
    double* A = new double[(size_t)Nrow * Mcol];
    double* x = new double[Mcol];
    double* y = new double[Nrow];

    // initialize 
    for (int i = 0; i < Nrow; i++) {
        double* row = A + (size_t)i * Mcol;
        #pragma GCC unroll 4
        for (int j = 0; j < Mcol; j++) {
            row[j] = drand48();
        }
    }
    for (int j = 0; j < Mcol; j++) x[j] = drand48();

    // timing
    matrix_vect_1d(A, x, y, Nrow, Mcol); // warm up

    double t0 = wall_time();
    for (int iter = 0; iter < 3; ++iter) matrix_vect_1d(A, x, y, Nrow, Mcol);
    double t1 = wall_time();
    double elapsed = (t1 - t0) / 3.0; // average time

    // add check to avoid having all zeros
    double checksum = 0.0;
    for (int i = 0; i < Nrow; ++i) checksum += y[i];
    std::cout << "checksum=" << checksum << "\n";

    // FLOPs
    long double flops  = 2.0L * Nrow * (long double)Mcol;
    long double tflops = (flops / elapsed) / 1.0e12L;

    // Calculate Bytes Moved
    double bytes = sizeof(double) * ( (long double)Nrow * Mcol + (long double)Mcol + (long double)Nrow );
    double bandwidth = (bytes / elapsed) / 1.0e9L; //
    
    // calculate arithmetic intensity
    double ai = flops / bytes;

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "N=" << Nrow << " M=" << Mcol
                << " | time=" << elapsed << " s"
                << " | FLOPs=" << (long long)flops
                << " | Perf=" << tflops << " TFLOP/s"
                << " | Bandwidth=" << bandwidth << " GB/s"
                << " | AI=" << ai << " FLOP/Byte"
                << "\n";

    delete[] A;
    delete[] x;
    delete[] y;
}


int main() {
    std::vector<int> sizes = {10, 100, 1000, 10000};
    std::cout << "Continguous Memory Allocation\n";

    for (int n : sizes) {
        run_case(n, n); 
    }

    return 0;
}