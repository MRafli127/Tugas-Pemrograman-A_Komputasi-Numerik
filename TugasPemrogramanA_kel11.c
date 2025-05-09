#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_DATA 100

typedef struct {
    int year;
    double internet_pct;
    double population;
} DataPoint;

// Data dari Data_Tugas_Pemrograman_A.csv
DataPoint dataset[] = {
    {1990, 0.0, 178200000},
    {1991, 0.0, 181400000},
    {1992, 0.0, 184700000},
    {1993, 0.0, 188000000},
    {1994, 0.0, 191400000},
    {1995, 0.01, 194800000},
    {1996, 0.02, 198300000},
    {1997, 0.04, 201900000},
    {1998, 0.11, 205500000},
    {1999, 0.29, 209200000},
    {2000, 0.92, 212900000},
    {2001, 1.84, 216700000},
    {2002, 3.02, 220600000},
    {2003, 4.42, 224500000},
    {2004, 6.28, 228400000},
    {2007, 8.49, 236400000},
    {2008, 10.92, 240500000},
    {2009, 12.52, 244600000},
    {2010, 13.94, 248800000},
    {2011, 15.36, 253000000},
    {2012, 16.72, 257300000},
    {2013, 18.15, 261600000},
    {2014, 19.56, 265900000},
    {2017, 25.38, 278700000},
    {2018, 30.05, 283200000},
    {2019, 34.29, 287700000},
    {2020, 39.03, 292200000}
};
int data_count = sizeof(dataset) / sizeof(dataset[0]);

// Fungsi untuk interpolasi kubik
double cubic_interpolation(double x0, double y0, double x1, double y1,
                           double x2, double y2, double x3, double y3, double x) {
    double t = (x - x1)/(x2 - x1);
    double a = -0.5*y0 + 1.5*y1 - 1.5*y2 + 0.5*y3;
    double b = y0 - 2.5*y1 + 2.0*y2 - 0.5*y3;
    double c = -0.5*y0 + 0.5*y2;
    double d = y1;
    
    return a*t*t*t + b*t*t + c*t + d;
}

// Fungsi untuk interpolasi spline kubik
double spline_interpolation(int x[], double y[], int n, double xi) {
    // Cek jika tahun yang dicari di luar range data
    if (xi <= x[0]) return y[0];
    if (xi >= x[n-1]) return y[n-1];
    
    // Cari indeks tahun yang tepat
    int i;
    for (i = 0; i < n-1; i++) {
        if (xi >= x[i] && xi <= x[i+1]) {
            break;
        }
    }
    
    // Pastikan kita memiliki cukup titik untuk interpolasi kubik
    if (i == 0) {
        // Gunakan 4 titik pertama
        return cubic_interpolation(x[0], y[0], x[1], y[1], 
                                  x[2], y[2], x[3], y[3], xi);
    } else if (i >= n-3) {
        // Gunakan 4 titik terakhir
        return cubic_interpolation(x[n-4], y[n-4], x[n-3], y[n-3],
                                  x[n-2], y[n-2], x[n-1], y[n-1], xi);
    } else {
        // Gunakan 2 titik di kiri dan 2 titik di kanan
        return cubic_interpolation(x[i-1], y[i-1], x[i], y[i],
                                  x[i+1], y[i+1], x[i+2], y[i+2], xi);
    }
}

// Fungsi untuk fitting polinomial menggunakan metode least squares
void polynomial_fit(int x[], double y[], int n, int degree, double coeffs[]) {
    double X[2*degree+1];  // Array untuk menyimpan sigma(x^j)
    double B[degree+1];    // Array untuk sigma(y*x^j)
    
    // Inisialisasi arrays
    for (int i = 0; i <= 2*degree; i++) {
        X[i] = 0;
        if (i <= degree) {
            B[i] = 0;
        }
    }
    
    // Calculating the X and B arrays
    for (int i = 0; i < n; i++) {
        double temp = 1;
        for (int j = 0; j <= degree; j++) {
            B[j] += y[i] * temp;
            temp *= x[i];
        }
        
        temp = 1;
        for (int j = 0; j <= 2*degree; j++) {
            X[j] += temp;
            temp *= x[i];
        }
    }
    
    // Mengatur nilai normal equations dalam bentuk matriks
    double A[degree+1][degree+1];
    for (int i = 0; i <= degree; i++) {
        for (int j = 0; j <= degree; j++) {
            A[i][j] = X[i+j];
        }
    }
    
    // Menyelesaikan persamaan dengan metode Gaussian elimination
    for (int i = 0; i <= degree; i++) {
        for (int j = i+1; j <= degree; j++) {
            double ratio = A[j][i] / A[i][i];
            
            for (int k = i; k <= degree; k++) {
                A[j][k] -= ratio * A[i][k];
            }
            
            B[j] -= ratio * B[i];
        }
    }
    
    // Back substitution
    for (int i = degree; i >= 0; i--) {
        coeffs[i] = B[i];
        for (int j = i+1; j <= degree; j++) {
            coeffs[i] -= A[i][j] * coeffs[j];
        }
        coeffs[i] /= A[i][i];
    }
}

// Fungsi untuk mengevaluasi polinomial pada nilai x tertentu
double evaluate_poly(double coeffs[], int degree, double x) {
    double result = 0;
    double temp = 1;
    
    for (int i = 0; i <= degree; i++) {
        result += coeffs[i] * temp;
        temp *= x;
    }
    
    return result;
}

// Fungsi untuk mengembalikan persamaan polinomial dalam format string
void get_polynomial_equation(double coeffs[], int degree, char* equation) {
    char temp[100];
    sprintf(equation, "y = ");
    
    for (int i = degree; i >= 0; i--) {
        if (fabs(coeffs[i]) < 1e-10) continue; // Skip koefisien yang sangat kecil
        
        if (i < degree && coeffs[i] > 0) {
            strcat(equation, " + ");
        } else if (i < degree && coeffs[i] < 0) {
            strcat(equation, " - ");
            coeffs[i] = -coeffs[i]; // Ubah tanda untuk output
        }
        
        if (i > 0) {
            if (fabs(coeffs[i] - 1.0) < 1e-10) {
                // Koefisien 1, cukup tulis x^i
                if (i == 1) {
                    sprintf(temp, "x");
                } else {
                    sprintf(temp, "x^%d", i);
                }
            } else {
                // Koefisien bukan 1
                if (i == 1) {
                    sprintf(temp, "%.6gx", coeffs[i]);
                } else {
                    sprintf(temp, "%.6gx^%d", coeffs[i], i);
                }
            }
        } else {
            // Konstanta
            sprintf(temp, "%.6g", coeffs[i]);
        }
        
        strcat(equation, temp);
    }
}

// Fungsi untuk menulis hasil ke file CSV
void write_csv(const char* filename, double pop_2005, double pop_2006, double pop_2015, double pop_2016,
               double internet_2005, double internet_2006, double internet_2015, double internet_2016,
               double pop_2030, double internet_2035) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        perror("Gagal membuat file output");
        exit(1);
    }
    
    // Header - Menggunakan format yang lebih deskriptif
    fprintf(file, "Year,Percentage_Internet_User,Population\n");
    
    // Data interpolasi untuk tahun 2005, 2006, 2015, 2016
    fprintf(file, "2005,%.4f,%.0f\n", internet_2005, pop_2005);
    fprintf(file, "2006,%.4f,%.0f\n", internet_2006, pop_2006);
    fprintf(file, "2015,%.4f,%.0f\n", internet_2015, pop_2015);
    fprintf(file, "2016,%.4f,%.0f\n", internet_2016, pop_2016);
    
    // Spasi pemisah antara data interpolasi dan prediksi
    fprintf(file, "\n");
    
    // Buat bagian khusus untuk prediksi
    fprintf(file, "# Prediksi Masa Depan\n");
    
    // Prediksi populasi tahun 2030
    fprintf(file, "Population_2030,,%.0f\n", pop_2030);
    
    // Prediksi persentase pengguna internet tahun 2035
    // Pastikan tidak melebihi 100%
    if (internet_2035 > 100.0) internet_2035 = 100.0;
    fprintf(file, "Internet_Percentage_2035,%.4f,\n", internet_2035);
    
    fclose(file);
    printf("Hasil berhasil disimpan di '%s'\n", filename);
}

int main() {
    int years[MAX_DATA], internet_years[MAX_DATA];
    double populations[MAX_DATA], internet_pcts[MAX_DATA];
    int internet_count = 0;
    
    // Mengekstrak data tahun, populasi, dan persentase internet
    for (int i = 0; i < data_count; i++) {
        years[i] = dataset[i].year;
        populations[i] = dataset[i].population;
        
        // Ekstra hati-hati dengan data internet, hanya ambil yang ada nilainya
        if (dataset[i].internet_pct > 0) {
            internet_years[internet_count] = dataset[i].year;
            internet_pcts[internet_count] = dataset[i].internet_pct;
            internet_count++;
        }
    }
    
    // 1. Interpolasi untuk nilai yang hilang
    double pop_2005 = spline_interpolation(years, populations, data_count, 2005);
    double pop_2006 = spline_interpolation(years, populations, data_count, 2006);
    double pop_2015 = spline_interpolation(years, populations, data_count, 2015);
    double pop_2016 = spline_interpolation(years, populations, data_count, 2016);
    
    double internet_2005 = spline_interpolation(internet_years, internet_pcts, internet_count, 2005);
    double internet_2006 = spline_interpolation(internet_years, internet_pcts, internet_count, 2006);
    double internet_2015 = spline_interpolation(internet_years, internet_pcts, internet_count, 2015);
    double internet_2016 = spline_interpolation(internet_years, internet_pcts, internet_count, 2016);
    
    // 2. Fitting polinomial
    // Untuk populasi tetap gunakan derajat 3
    double pop_coeffs[4] = {0};
    polynomial_fit(years, populations, data_count, 3, pop_coeffs);
    
    // Untuk internet gunakan derajat lebih rendah (2) untuk menghindari overfitting
    double internet_coeffs[3] = {0};
    polynomial_fit(internet_years, internet_pcts, internet_count, 2, internet_coeffs);
    
    // 3. Estimasi untuk tahun 2030 dan 2035
    double pop_2030 = evaluate_poly(pop_coeffs, 3, 2030);
    double internet_2035 = evaluate_poly(internet_coeffs, 2, 2035);
    
    // Batasi persentase internet maksimum 100%
    if (internet_2035 > 100.0) internet_2035 = 100.0;
    
    // Menampilkan hasil
    printf("1. Nilai yang hilang:\n");
    printf("   a. Populasi Indonesia 2005: %.0f orang\n", pop_2005);
    printf("   b. Populasi Indonesia 2006: %.0f orang\n", pop_2006);
    printf("   c. Populasi Indonesia 2015: %.0f orang\n", pop_2015);
    printf("   d. Populasi Indonesia 2016: %.0f orang\n", pop_2016);
    printf("   e. Persentase Internet Indonesia 2005: %.4f%%\n", internet_2005);
    printf("   f. Persentase Internet Indonesia 2006: %.4f%%\n", internet_2006);
    printf("   g. Persentase Internet Indonesia 2015: %.4f%%\n", internet_2015);
    printf("   h. Persentase Internet Indonesia 2016: %.4f%%\n", internet_2016);
    
    // Dapatkan persamaan polinomial
    char internet_equation[200], population_equation[200];
    get_polynomial_equation(internet_coeffs, 2, internet_equation);
    get_polynomial_equation(pop_coeffs, 3, population_equation);
    
    printf("\n2. Persamaan polinomial:\n");
    printf("   a. Persentase Internet: %s\n", internet_equation);
    printf("   b. Populasi Indonesia: %s\n", population_equation);
    
    printf("\n3. Estimasi:\n");
    printf("   a. Populasi Indonesia 2030: %.0f orang\n", pop_2030);
    printf("   b. Persentase Internet Indonesia 2035: %.4f%%\n", internet_2035);
    
    // Tulis hasil ke file CSV
    write_csv("hasil_estimasi.csv", pop_2005, pop_2006, pop_2015, pop_2016,
              internet_2005, internet_2006, internet_2015, internet_2016,
              pop_2030, internet_2035);
    
    return 0;
}