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

// calculation
void matrix_vect_sep(double** A, double* x, double* y, int Nrow, int Mcol) {
    for (int i = 0; i < Nrow; ++i) {
        double s = 0.0;

        #pragma GCC unroll 4
        for (int j = 0; j < Mcol; ++j) {
            s += A[i][j] * x[j]; // 1 mul + 1 add = 2 FLOPs
        }
        y[i] = s;
    }
}

// init, timing, measure, and delete
void run_case_separate(int Nrow, int Mcol) {
    // allocate separate rows
    double** A = new double*[Nrow];
    for (int i = 0; i < Nrow; i++) {
        A[i] = new double[Mcol];
    }
    double* x = new double[Mcol];
    double* y = new double[Nrow];

    // initialize 
    for (int i = 0; i < Nrow; i++) {
        for (int j = 0; j < Mcol; j++) {
            A[i][j] = drand48();
        }
    }
    for (int j = 0; j < Mcol; j++) x[j] = drand48();

    // timing
    matrix_vect_sep(A, x, y, Nrow, Mcol); // warm up

    double t0 = wall_time();
    for (int iter = 0; iter < 3; ++iter) matrix_vect_sep(A, x, y, Nrow, Mcol);
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

    // free memory
    for (int i = 0; i < Nrow; i++) delete[] A[i];
    delete[] A;
    delete[] x;
    delete[] y;
}

int main() {
    std::vector<int> sizes = {10,100,1000,10000};

    std::cout << "Separate Memory Allocation\n";
    for (int n : sizes) {
        run_case_separate(n, n); 
    }

    return 0;
}