#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv){

    const char* fin = argv[1];
    const char* fout = argv[2];
    FILE* in = fopen(fin, "r");
    int N ;
    double K1A;

    fscanf(in, "%d ", &N);

    double* diag = (double*)malloc(N * sizeof(double));
    double* diag_inf = (double*)malloc((N-1) * sizeof(double));
    double* diag_sup = (double*)malloc((N-1) * sizeof(double));
    double* vecteur_b = (double*)malloc((N)* sizeof(double));

    for (int i = 0; i < N; i++) fscanf(in, "%lf", &diag[i]);
    for (int i = 0; i < N - 1; i++) fscanf(in, "%lf", &diag_inf[i]);
    for (int i = 0; i < N - 1; i++) fscanf(in, "%lf", &diag_sup[i]);
    for (int i = 0; i < N; i++) fscanf(in, "%lf", &vecteur_b[i]);


    fscanf(in, "%lf", &K1A);
    fclose(in);

    // Factorisation LU
    double* l = (double*)malloc((N-1) * sizeof(double));
    double* d = (double*)malloc(N * sizeof(double));
    
    d[0] = diag[0];  
    for (int i = 1; i < N; i++){
        l[i-1] = diag_inf[i-1] / d[i-1];
        d[i] = diag[i] - l[i-1] * diag_sup[i-1];

    }

    // Descente : Ly = b
    double* y = (double*)malloc(N* sizeof(double));
    y[0] = vecteur_b[0];
    for (int i = 1; i < N; i++) y[i] = vecteur_b[i] - l[i-1] * y[i-1];

    // Remonte : Ux = y
    double* x = (double*)malloc(N* sizeof(double));
    x[N-1] = y[N-1] / d[N-1];
    for (int i = N-2; i >= 0; i--) x[i] = (y[i] - diag_sup[i] * x[i+1]) / d[i];

    FILE* out = fopen(fout, "w");

    fprintf(out, "N = %d\n", N);
    fprintf(out, "K1(A) = %.14f\n\n", K1A);

    //Afficher A

    fprintf(out, "Matrice A :\n");
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            double val = 0.0;
            if (j == i) val = diag[i];
            else if (j == i-1) val = diag_inf[i-1];
            else if (j == i+1) val = diag_sup[i];
            fprintf(out, "%.15f ", val);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\n");

    // Afficher L

    fprintf(out, "Matrice L: \n");
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            double val = 0.0;
            if (j == i) val = 1.0;
            else if (j == i-1) val = l[i-1];
            fprintf(out, "%.15f ", val);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\n");

    // AFficher U

    fprintf(out, "Matrice U: \n");
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            double val = 0.0;
            if (j == i) val = d[i];
            else if (j == i+1) val = diag_sup[i];
            fprintf(out, "%.15f ", val);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "\n");

    // Afficher X 
    fprintf(out, "Solution numerique de X: \n");
    for (int i = 0; i < N; i++) fprintf(out, "%.15f ", x[i]);
    fprintf(out, "\n");

    fclose(out);
    free(diag); 
    free(diag_inf);
    free(diag_sup);
    free(vecteur_b);
    free(x);
    free(l);
    free(d);
    free(y);
    return 0;

}