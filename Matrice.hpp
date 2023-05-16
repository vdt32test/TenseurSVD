#pragma once
#include "Vecteur.hpp"


class Matrice
{
public:
    Matrice();
    Matrice(int lig, int col);
    Matrice(const Vecteur &vec);
    Matrice(const Vecteur *vecs, int col);
    Matrice(const Matrice &mat);
    ~Matrice();

    Matrice operator=(const Matrice &mat);
    Matrice operator+(const Matrice &vec) const;
    Matrice operator-(const Matrice &vec) const;
    Matrice operator*(const Matrice &vec) const;
    Vecteur &operator[](int index);
    const Vecteur &operator[](int index) const;

    Vecteur mvprod(const Vecteur &vec) const;
    Matrice transpose() const;
    Matrice submat(int i1, int i2, int j1, int j2) const;
    int getLigDim() const;
    int getColDim() const;
    bool isDiagonale() const;
    bool isTridiagonale() const;
    void affiche() const;
    void setValues(const Matrice &mat, int i, int j);
    void swap(int index_lig_begin_1, int index_lig_end_1, int index_col_begin_1, 
        int index_col_end_1, int index_lig_begin_2, int index_col_begin_2);
    void zeroLowValues(float epsilon = EPSILON);

    friend float norm(const Matrice &mat);
    friend Matrice operator*(float numb, const Matrice &mat);
    friend Matrice outer(const Vecteur &vec1, const Vecteur &vec2);
    //friend void swap(Matrice &mat1, Matrice &mat2, int i1, int i2, int j1, int j1);
    friend void swapMats(Matrice &mat1, Matrice &mat2);


private:
    Vecteur *m_mat;
    int m_dims[2];
};

Matrice vecToMat(const Vecteur &vec);