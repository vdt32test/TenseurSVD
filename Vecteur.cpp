#include "Vecteur.hpp"


Vecteur::Vecteur():
    m_tab(nullptr),
    m_dim(0)
{
    //std::cout << "Vecteur::Vecteur()" << std::endl;
    //std::cout << this << std::endl;
    //std::cout << "Vecteur::Vecteur() end" << std::endl;
}


Vecteur::Vecteur(int dim, float default_value):
    m_tab(nullptr),
    m_dim(dim)
{
    //std::cout << "Vecteur::Vecteur(int dim, float default_value)" << std::endl;
    //std::cout << this << std::endl;

    assert(dim != 0);

    m_tab = new float[m_dim]();
    if(default_value != 0)
    {
        for(int i = 0; i < m_dim; i++)
        {
            m_tab[i] = default_value;
        }
    }
    
    //std::cout << "Vecteur::Vecteur(int dim, float default_value) end" << std::endl;
}


Vecteur::Vecteur(const float *vec, int dim):
    m_tab(nullptr),
    m_dim(dim)
{
    //std::cout << "Vecteur::Vecteur(const float *vec, int dim)" << std::endl;
    //std::cout << this << std::endl;

    m_tab = new float[m_dim]();
    for(int i = 0; i < m_dim; i++)
    {
        m_tab[i] = vec[i];
    }
    
    //std::cout << "Vecteur::Vecteur(const float *vec, int dim) end" << std::endl;
}


Vecteur::Vecteur(const Vecteur &vec):
    m_tab(nullptr),
    m_dim(vec.m_dim)
{
    //std::cout << "Vecteur::Vecteur(const Vecteur &vec)" << std::endl;
    //std::cout << this << std::endl;

    m_tab = new float[m_dim]();
    for(int i = 0; i < m_dim; i++)
    {
        m_tab[i] = vec.m_tab[i];
    }
    
    //std::cout << "Vecteur::Vecteur(const Vecteur &vec) end" << std::endl;
}


Vecteur::~Vecteur()
{
    //std::cout << "Vecteur::~Vecteur() begin" << m_tab << std::endl;
    //std::cout << this << std::endl;
    //std::cout << "m_tab = " << m_tab << std::endl;
    
    if(m_dim != 0)
    {
        assert(m_tab != nullptr);
        delete[] m_tab;
    }
    else
    {
        assert(m_tab == nullptr);
    }
    //std::cout << "Vecteur::~Vecteur() end" << std::endl;
}


Vecteur Vecteur::operator=(const Vecteur &vec)
{
    //std::cout << "Vecteur::operator=(const Vecteur &vec)" << std::endl;

    Vecteur tmp(vec);
    swap(*this, tmp);
    
    //std::cout << "Vecteur::operator=(const Vecteur &vec) end" << std::endl;
    return *this;
}


Vecteur Vecteur::operator+(const Vecteur &vec)
{
    //std::cout << "Vecteur::operator+(const Vecteur &vec)" << std::endl;

    assert(m_dim == vec.m_dim);

    Vecteur new_vec(m_dim);
    for(int i = 0; i < m_dim; i++)
    {
        new_vec.m_tab[i] = m_tab[i]+vec.m_tab[i];
    }

    //std::cout << "Vecteur::operator+(const Vecteur &vec) end" << std::endl;
    return new_vec;
}


Vecteur Vecteur::operator-(const Vecteur &vec)
{
    //std::cout << "Vecteur::operator-(const Vecteur &vec)" << std::endl;

    assert(m_dim == vec.m_dim);

    Vecteur new_vec(m_dim);
    for(int i = 0; i < m_dim; i++)
    {
        new_vec.m_tab[i] = m_tab[i]-vec.m_tab[i];
    }

    //std::cout << "Vecteur::operator-(const Vecteur &vec) end" << std::endl;
    return new_vec;
}


float &Vecteur::operator[](int index)
{
    //std::cout << "float &Vecteur::operator[](int index) begin" << std::endl;

    //std::cout << "index = " << index << std::endl;
    TEST_ASSERT(index, >=, 0);
    TEST_ASSERT(index, <, m_dim);
    
    //std::cout << "float &Vecteur::operator[](int index) end" << std::endl;

    return m_tab[index];
}


const float &Vecteur::operator[](int index) const
{
    //std::cout << "const float &Vecteur::operator[](int index) const begin" << std::endl;

    TEST_ASSERT(index, >=, 0);
    TEST_ASSERT(index, <, m_dim);
    
    //std::cout << "const float &Vecteur::operator[](int index) const end" << std::endl;

    return m_tab[index];
}


Vecteur Vecteur::subvec(int i, int j) const
{
    //std::cout << "Vecteur Vecteur::subvec(int i, int j) const begin" << std::endl;

    assert(i >= 0 && i <= j && j < m_dim);

    Vecteur vec(j-i+1);

    for(int k = i; k <= j; k++)
    {
        vec.m_tab[k-i] = m_tab[k];
    }
    
    //std::cout << "Vecteur Vecteur::subvec(int i, int j) const end" << std::endl;

    return vec;
}


float Vecteur::getMin() const
{
    //std::cout << "float Vecteur::getMin() const begin" << std::endl;

    float min = m_tab[0];
    for(int i = 1; i < m_dim; i++)
    {
        if(m_tab[i] < min)
        {
            min = m_tab[i];
        }
    }
    
    //std::cout << "float Vecteur::getMin() const end" << std::endl;

    return min;
}


float Vecteur::getMax() const
{
    //std::cout << "float Vecteur::getMax() const begin" << std::endl;

    float max = m_tab[0];
    for(int i = 1; i < m_dim; i++)
    {
        if(m_tab[i] > max)
        {
            max = m_tab[i];
        }
    }
    
    //std::cout << "float Vecteur::getMax() const end" << std::endl;

    return max;
}


/*
float Vecteur::getMin(int begin_index, int end_index, int &min_index) const
{
    assert(begin_index >= 0 && begin_index <= end_index && end_index < m_dim);

    min_index = begin_index;
    float min = m_tab[begin_index];
    for(int i = begin_index+1; i <= end_index; i++)
    {
        if(m_tab[i] < min)
        {
            min = m_tab[i];
            min_index = i;
        }
    }

    return min;
}


float Vecteur::getMax(int begin_index, int end_index, int &max_index) const
{
    assert(begin_index >= 0 && begin_index <= end_index && end_index < m_dim);

    max_index = begin_index;
    float max = m_tab[begin_index];
    for(int i = begin_index+1; i <= end_index; i++)
    {
        if(m_tab[i] > max)
        {
            max = m_tab[i];
            max_index = i;
        }
    }

    return max;
}
*/


void Vecteur::affiche() const
{
    //std::cout << "void Vecteur::affiche() const begin" << std::endl;

    std::cout << "dim = " << m_dim << std::endl;
    for(int i = 0; i < m_dim; i++)
    {
        std::cout << "tab[" << i << "] = " << m_tab[i] << std::endl;
    }
    
    //std::cout << "void Vecteur::affiche() const end" << std::endl;
}


void Vecteur::setValues(const Vecteur &vec, int i)
{
    //std::cout << "void Vecteur::setValues(const Vecteur &vec, int i) begin" << std::endl;

    assert(m_dim >= i + vec.getDim());

    for(int j = 0; j < vec.getDim(); j++)
    {
        m_tab[j+i] = vec[j];
    }
    
    //std::cout << "void Vecteur::setValues(const Vecteur &vec, int i) end" << std::endl;
}


int Vecteur::getDim() const
{
    //std::cout << "int Vecteur::getDim() const begin" << std::endl;
    //std::cout << "int Vecteur::getDim() const end" << std::endl;

    return m_dim;
}


float dot(const Vecteur &vec1, const Vecteur &vec2)
{
    //std::cout << "float dot(const Vecteur &vec1, const Vecteur &vec2) begin" << std::endl;

    assert(vec1.m_dim == vec2.m_dim);
    float dot_pr = 0;
    
    for(int i = 0; i < vec1.m_dim; i++)
    {
        dot_pr += vec1.m_tab[i]*vec2.m_tab[i];
    }
    
    //std::cout << "float dot(const Vecteur &vec1, const Vecteur &vec2) end" << std::endl;

    return dot_pr;
}


float norm(const Vecteur &vec)
{
    //std::cout << "float norm(const Vecteur &vec) begin" << std::endl;

    float vec_norm = 0;

    for(int i = 0; i < vec.m_dim; i++)
    {
        vec_norm += vec.m_tab[i]*vec.m_tab[i];
    }

    vec_norm = sqrt(vec_norm);

    //std::cout << "float norm(const Vecteur &vec) end" << std::endl;
    return vec_norm;
}


Vecteur operator*(float numb, const Vecteur &vec)
{
    //std::cout << "Vecteur operator*(float numb, const Vecteur &vec) begin" << std::endl;

    Vecteur new_vec(vec);
    for(int i = 0; i < vec.m_dim; i++)
    {
        new_vec.m_tab[i] *= numb;
    }

    //std::cout << "Vecteur operator*(float numb, const Vecteur &vec) end" << std::endl;
    return new_vec;
}


void swap(Vecteur &vec1, Vecteur &vec2)
{
    std::swap(vec1.m_tab, vec2.m_tab);
    std::swap(vec1.m_dim, vec2.m_dim);
}