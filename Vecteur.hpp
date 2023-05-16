#pragma once
#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdio>
#include "utils.hpp"


class Vecteur
{
public:
    Vecteur();
    Vecteur(int dim, float default_value = 0);
    Vecteur(const float *vec, int dim);
    Vecteur(const Vecteur &vec);
    ~Vecteur();

    Vecteur operator=(const Vecteur &vec);
    Vecteur operator+(const Vecteur &vec);
    Vecteur operator-(const Vecteur &vec);
    float &operator[](int index);
    const float &operator[](int index) const;

    Vecteur subvec(int i, int j) const;
    float getMin() const;
    float getMax() const;
    //float getMin(int begin_index, int end_index, int &min_index) const;
    //float getMax(int begin_index, int end_index, int &max_index) const;
    int getDim() const;
    void affiche() const;
    void setValues(const Vecteur &vec, int i);

    friend float dot(const Vecteur &vec1, const Vecteur &vec2);
    friend float norm(const Vecteur &vec);
    friend Vecteur operator*(float numb, const Vecteur &vec);
    friend void swap(Vecteur &vec1, Vecteur &vec2);

private:
    float *m_tab;
    int m_dim;
};
