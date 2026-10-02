#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrices_tools.h"

int main(int argc, char** argv){
    const char* fin = argv[1];
    const char* fout = argv[2];

    FILE* in = fopen(fin, "r");

    int n, N;
    fscanf(in, "%d %d", &N, &n);

    double** B = malloc(N* sizeof(double*));
    double** A = malloc(N* sizeof(double*));
    double** C = malloc(N* sizeof(double*));

    for (int i = 0; i < N; i++){
        B[i] = malloc(n*n*sizeof(double));
        A[i] = malloc(n*n*sizeof(double));
        C[i] = malloc(n*n*sizeof(double));

    }

    for (int i = 0; i < N;i++) lire_matrice(in, B[i], n);
    for (int i = 1; i < N;i++) lire_matrice(in, A[i], n);
    for (int i = 0; i < N-1; i++) lire_matrice(in, C[i], n);

    double K1A;
    fscanf(in, "%lf", &K1A);
    fclose(in);

    double* ones = malloc(N*n*sizeof(double));
    double* b = malloc(N*n*sizeof(double));

    for (int i = 0; i < N*n; i++) ones[i] = 1.0;

    for (int i = 0; i < N; i++){
        double* bi = &b[i*n];
        for (int k = 0; k < n; k++) bi[k] = 0.0;

        if (i > 0){
            double tmpA[n];
            mat_vec(A[i], &ones[(i-1)*n], tmpA, n);
            for (int k = 0; k < n; k++) bi[k] += tmpA[k];
        }

        double tmpB[n];
        mat_vec(B[i], &ones[i*n], tmpB, n);
        for (int k = 0; k < n; k++) bi[k] += tmpB[k];

        if (i < N-1){
            double tmpC[n];
            mat_vec(C[i], &ones[(i+1)*n], tmpC, n);
            for (int k = 0; k < n; k++) bi[k] += tmpC[k];
        }
    }

    //Factorisation LU
    double** L = malloc(N * sizeof(double*));
    double** U = malloc(N * sizeof(double*));

    for (int i = 0; i < N; i++){
        L[i] = malloc(n*n*sizeof(double));
        U[i] = malloc(n*n*sizeof(double));
    }

    mat_copy(B[0], U[0], n);

    for(int i = 1; i < N; i++){
        double Uinv[n*n];
        mat_inv(U[i-1], Uinv, n);

        // Li = Ai * Uinv_i-1
        mat_mul(A[i], Uinv, L[i],n);

        //temp =Li * Ci-1
        double temp[n*n];
        mat_mul(L[i], C[i-1], temp, n);

        // Ui = Bi - temp
        mat_sub(B[i], temp, U[i], n);

    }

    double* y = malloc(N*n* sizeof(double));

    // y0 = b0
    for (int k = 0; k < n; k++) y[k] = b[k];
    
    for (int i = 1; i < N; i++){
        double tmp[n];
        mat_vec(L[i], &y[(i-1)*n], tmp, n);
        for (int k = 0; k < n; k++) y[i*n + k] = b[i*n + k] - tmp[k];
    }

    double* x = malloc(N*n * sizeof(double));

    double Uinv[n*n];
    mat_inv(U[N-1], Uinv, n);
    mat_vec(Uinv, &y[(N-1)*n], &x[(N-1)*n], n);

    for (int i = N-2; i >= 0; i--){
        double tmp[n];
        mat_vec(C[i], &x[(i+1)*n], tmp, n);
        
        double r[n];
        for (int k = 0; k < n; k++) r[k] = y[i*n + k] - tmp[k];

        double Uinv[n*n];
        mat_inv(U[i], Uinv, n);
        mat_vec(Uinv, r, &x[i*n], n);

    }

    double err = 0.0;
    for (int i = 0; i < N*n; i++){
        double d = x[i] - 1.0;
        err += d*d;
    }
    err = sqrt(err);

    FILE* out = fopen(fout, "w");
    fprintf(out, "N = %d, n = %d\n", N, n);
    fprintf(out, "K1A = %.15f\n\n", K1A);
    fprintf(out, "Erreur ||X-1||2 = %.14f\n\n", err);

    // Afficher Matrice A
    fprintf(out, "Matrice A par blocs: \n");

    for (int i = 0; i < N; i++){
        fprintf(out,"Ligne de blocs %d : \n", i+1);
        
        if (i > 0){
            fprintf(out, " A_%d :\n", i+1);
            afficher_matrice(out, A[i], n);
        }
        fprintf(out, " B_%d :\n",i+1);
        afficher_matrice(out, B[i], n);

        if (i < N-1){
            fprintf(out, " C_%d :\n", i+1);
            afficher_matrice(out, C[i], n);
        }
        fprintf(out,"\n");
    }

    // Afficher Matrice L
    fprintf(out, "Matrice L par blocs\n");
    for (int i = 0; i < N; i++){

        fprintf(out, "Bloc diagonal L_%d = Id : \n", i+1);
        double I[n*n];  
        mat_iden(I,n);
        afficher_matrice(out, I, n);
        fprintf(out, "\n");

        if (i > 0){
            fprintf(out, "Bloc sous-diagonale L_%d :\n", i+1);
            afficher_matrice(out, L[i], n);
            fprintf(out, "\n");
        }

    }
    // Afficher Matrice U
    fprintf(out, "Matrice U par blocs\n");

    for (int i = 0; i < N; i++){
        fprintf(out,"Bloc diagonal U_%d : \n", i+1);
        afficher_matrice(out, U[i], n);
        fprintf(out, "\n");

        if (i < N-1){
            fprintf(out, "Bloc sur-diagonal C_%d :\n", i+1);
            afficher_matrice(out, C[i],n);
            fprintf(out,"\n");
        }
    }

    //Afficher X
    fprintf(out, "Solution numérique X\n");

    for (int i = 0; i < N; i++){
        fprintf(out, "Bloc %d :\n", i+1);
        afficher_vecteur(out, &x[i*n], n);
        fprintf(out, "\n");

    }

    fclose(out);

    for (int i = 0; i < N; i++){
        free(A[i]);
        free(B[i]);
        free(C[i]);
        free(U[i]);
        free(L[i]);
    }
    free(A);
    free(B);
    free(C);
    free(L);
    free(U);
    free(ones);
    free(b);
    free(x);
    free(y);

    return 0;

}