#pragma once
#include "Matrice.hpp"

int sign(float x);
void givens(float x, float z, float &c, float &s);
Vecteur householder(const Vecteur &x, float &beta);
Matrice reductridiag(Matrice &D);
void qrsym(const Matrice &A, Matrice &Q);
Matrice qrpivot(const Matrice &A, Matrice &Q, Matrice *AOut = nullptr);
void svd(const Matrice &A, Matrice *svdMatrices);