#include "Vecteur.hpp"


Vecteur::Vecteur():
    m_tab(nullptr),
    m_dim(0)
{}


Vecteur::Vecteur(int dim, float default_value):
    m_tab(nullptr),
    m_dim(dim)
{
    assert(dim != 0);

    m_tab = new float[m_dim]();
    if(default_value != 0)
    {
        for(int i = 0; i < m_dim; i++)
        {
            m_tab[i] = default_value;
        }
    }
}


Vecteur::Vecteur(const float *vec, int dim):
    m_tab(nullptr),
    m_dim(dim)
{
    m_tab = new float[m_dim]();
    for(int i = 0; i < m_dim; i++)
    {
        m_tab[i] = vec[i];
    }
}


Vecteur::Vecteur(const Vecteur &vec):
    m_tab(nullptr),
    m_dim(vec.m_dim)
{
    m_tab = new float[m_dim]();
    for(int i = 0; i < m_dim; i++)
    {
        m_tab[i] = vec.m_tab[i];
    }
}


Vecteur::~Vecteur()
{
    if(m_dim != 0)
    {
        assert(m_tab != nullptr);
        delete[] m_tab;
    }
    else
    {
        assert(m_tab == nullptr);
    }
}


Vecteur Vecteur::operator=(const Vecteur &vec)
{
    Vecteur tmp(vec);
    swap(*this, tmp);
    
    return *this;
}


Vecteur Vecteur::operator+(const Vecteur &vec)
{
    assert(m_dim == vec.m_dim);

    Vecteur new_vec(m_dim);
    for(int i = 0; i < m_dim; i++)
    {
        new_vec.m_tab[i] = m_tab[i]+vec.m_tab[i];
    }

    return new_vec;
}


Vecteur Vecteur::operator-(const Vecteur &vec)
{
    assert(m_dim == vec.m_dim);

    Vecteur new_vec(m_dim);
    for(int i = 0; i < m_dim; i++)
    {
        new_vec.m_tab[i] = m_tab[i]-vec.m_tab[i];
    }

    return new_vec;
}


float &Vecteur::operator[](int index)
{
    TEST_ASSERT(index, >=, 0);
    TEST_ASSERT(index, <, m_dim);

    return m_tab[index];
}


const float &Vecteur::operator[](int index) const
{
    TEST_ASSERT(index, >=, 0);
    TEST_ASSERT(index, <, m_dim);

    return m_tab[index];
}


Vecteur Vecteur::subvec(int i, int j) const
{
    assert(i >= 0 && i <= j && j < m_dim);

    Vecteur vec(j-i+1);

    for(int k = i; k <= j; k++)
    {
        vec.m_tab[k-i] = m_tab[k];
    }

    return vec;
}


float Vecteur::getMin() const
{
    float min = m_tab[0];
    for(int i = 1; i < m_dim; i++)
    {
        if(m_tab[i] < min)
        {
            min = m_tab[i];
        }
    }

    return min;
}


float Vecteur::getMax() const
{
    float max = m_tab[0];
    for(int i = 1; i < m_dim; i++)
    {
        if(m_tab[i] > max)
        {
            max = m_tab[i];
        }
    }

    return max;
}


void Vecteur::affiche() const
{
    std::cout << "dim = " << m_dim << std::endl;
    for(int i = 0; i < m_dim; i++)
    {
        std::cout << "tab[" << i << "] = " << m_tab[i] << std::endl;
    }
}


void Vecteur::setValues(const Vecteur &vec, int i)
{
    assert(m_dim >= i + vec.getDim());

    for(int j = 0; j < vec.getDim(); j++)
    {
        m_tab[j+i] = vec[j];
    }
}


int Vecteur::getDim() const
{
    return m_dim;
}


float dot(const Vecteur &vec1, const Vecteur &vec2)
{
    assert(vec1.m_dim == vec2.m_dim);
    float dot_pr = 0;
    
    for(int i = 0; i < vec1.m_dim; i++)
    {
        dot_pr += vec1.m_tab[i]*vec2.m_tab[i];
    }

    return dot_pr;
}


float norm(const Vecteur &vec)
{
    float vec_norm = 0;

    for(int i = 0; i < vec.m_dim; i++)
    {
        vec_norm += vec.m_tab[i]*vec.m_tab[i];
    }

    vec_norm = sqrt(vec_norm);

    return vec_norm;
}


Vecteur operator*(float numb, const Vecteur &vec)
{
    Vecteur new_vec(vec);
    for(int i = 0; i < vec.m_dim; i++)
    {
        new_vec.m_tab[i] *= numb;
    }

    return new_vec;
}


void swap(Vecteur &vec1, Vecteur &vec2)
{
    std::swap(vec1.m_tab, vec2.m_tab);
    std::swap(vec1.m_dim, vec2.m_dim);
}