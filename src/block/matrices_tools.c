#include "matrices_tools.h"

void mat_copy(const double* src, double* dest, int n){
    for (int i = 0; i < n * n; i++) dest[i] = src[i];
}

void mat_add(const double* A, const double* B, double* C, int n){
    for (int i = 0; i < n * n; i++) C[i] = A[i] + B[i];
}

void mat_sub(const double* A, const double* B, double* C, int n){
    for (int i = 0; i < n * n; i++) C[i] = A[i] - B[i];
}

void mat_mul(const double* A, const double* B, double* C, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            double ele = 0.0;
            for (int k = 0; k < n; k++) ele += A[i*n + k] * B[k*n + j];
            C[i*n+j] = ele;
        }
    }
}

// Produit y = A * x
void mat_vec(const double* A, const double* x, double* y, int n){
    for (int i = 0; i < n ; i++){
        double ele = 0.0;
        for (int j = 0; j < n; j++) ele += A[i*n + j] * x[j];
        y[i] = ele;
    }
}

void mat_iden(double* I, int n){
    for (int i = 0; i < n * n; i++) I[i] = 0.0;
    for (int i = 0; i < n; i++) I[i * n + i] = 1.0;
}

void mat_inv(const double* A, double* Ainv, int n){
    double* aug = malloc(n * 2 * n * sizeof(double));

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            aug[i*(2 * n) + j] = A[i*n + j];
            aug[i*(2 * n) + n +j] = (i==j) ? 1.0 : 0.0;
        }
    }

    for (int k = 0; k < n; k++){
        double pivot = aug[k*(2*n) + k];

        for (int j = 0; j < 2 * n; j++) aug[k * (2*n) + j] /= pivot;

        for (int i = 0; i < n; i++){
            if (i == k) continue;
            double f = aug[i*(2*n) + k];
            for (int j = 0; j < 2*n; j++) aug[i*(2*n)+j] -= f*aug[k*(2*n)+j];
        }
    }

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
        Ainv[i*n +j] = aug[i*(2*n) + n+ j];
    }

    free(aug);
}

void lire_matrice(FILE* f, double* M, int n){
    for (int i = 0; i < n*n; i++) fscanf(f, "%lf", &M[i]);
}

void afficher_matrice(FILE* f, const double* M, int n){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) fprintf(f, "%.15f ", M[i*n +j]);
        fprintf(f, "\n");
    }
}

void afficher_vecteur(FILE* f, const double* v, int n){
    for (int i = 0; i < n; i++) fprintf(f, "%.15f ", v[i]);
    fprintf(f, "\n");
}