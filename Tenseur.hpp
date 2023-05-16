#pragma once
#include "Matrice.hpp"


class Tenseur
{
public:
    Tenseur();
    Tenseur(const int *dims, int d);
    Tenseur(const int *dims, int d, const Vecteur &vecTen);
    Tenseur(const int *dims, int d, int k, const Matrice &matTen);
    Tenseur(const Tenseur &ten);
    ~Tenseur();

    Tenseur operator=(const Tenseur &ten);
    float operator[](int i) const;
    float &operator[](int i);
    float operator[](int *indexes) const;
    float &operator[](int *indexes);

    Tenseur operator+(const Tenseur &ten) const;
    Tenseur operator-(const Tenseur &ten) const;

    Matrice mode(int k) const;
    int getD() const;

    void initDims(const int *dims, int d);
    void affiche() const;
    void deinit();

    friend Tenseur pmod(const Tenseur &ten, const Matrice &mat, int k);
    friend void swapTens(Tenseur &ten1, Tenseur &ten2);

protected:
    int m_d;
    int *m_dims;
    int m_nbelts;
    Vecteur m_ten;
};


int tenToVecIndex(const int *dims, int d, const int *tenIndexes);
void vecToTenIndex(const int *dims, int d, int vecIndex, int *tenIndexes);
int tenToMatIndex(const int *dims, int d, const int *tenIndexes, int k);
int indexDivision(int a, int b);
