#ifndef MATRICE_TOOLS_H
#define MATRICE_TOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void mat_copy(const double* src, double* dest, int n);
void mat_add(const double* A, const double* B, double* C, int n);
void mat_sub(const double* A, const double* B, double* C, int n);
void mat_mul(const double* A, const double* B, double* C, int n);
void mat_vec(const double* A, const double* x, double* y, int n);
void mat_iden(double* I, int n);
void mat_inv(const double *A, double* Ainv, int n);
void lire_matrice(FILE* f, double* M, int n);
void afficher_matrice(FILE* f, const double* M, int n);
void afficher_vecteur(FILE* f, const double* v, int n);
#endif