#include "Matrice.hpp"


Matrice::Matrice():
    m_mat(nullptr),
    m_dims()
{
    //std::cout << __PRETTY_FUNCTION__ << std::endl;
}


Matrice::Matrice(int lig, int col):
    m_mat(nullptr),
    m_dims()
{
    //std::cout << __PRETTY_FUNCTION__ << std::endl;

    m_dims[0] = lig;
    m_dims[1] = col;
    //std::cout << "dims[0] = " << m_dims[0] << "; dims[1] = " << m_dims[1] << std::endl;

    if(lig == 0)
    {
        assert(col == 0);
        return;
    }

    Vecteur tmp(lig);
    //tmp.affiche();
    m_mat = new Vecteur[col];
    for(int i = 0; i < col; i++)
    {   
        m_mat[i] = tmp;
    }
}


Matrice::Matrice(const Vecteur &vec):
    m_mat(nullptr),
    m_dims()
{
    //std::cout << __PRETTY_FUNCTION__ << std::endl;
    
    int size = vec.getDim();
    m_dims[0] = m_dims[1] = size;

    Vecteur tmp(size);
    m_mat = new Vecteur[size];
    for(int i = 0; i < size; i++)
    {   
        m_mat[i] = tmp;
        m_mat[i][i] = vec[i];
    }
}


Matrice::Matrice(const Vecteur *vecs, int col):
    m_mat(nullptr),
    m_dims()
{
    //std::cout << __PRETTY_FUNCTION__ << std::endl;

    assert(vecs != nullptr && col != 0);
    m_dims[0] = vecs[0].getDim();
    m_dims[1] = col;

    m_mat = new Vecteur[col];
    for(int i = 0; i < col; i++)
    {   
        assert(vecs[i].getDim() == m_dims[0]);
        m_mat[i] = vecs[i];
    }   
}


Matrice::Matrice(const Matrice &mat):
    m_mat(nullptr),
    m_dims()
{
    //std::cout << __PRETTY_FUNCTION__ << std::endl;

    m_dims[0] = mat.m_dims[0];
    m_dims[1] = mat.m_dims[1];
    if(m_dims[0] == 0)
    {
        assert(m_dims[1] == 0);
        return;
    }

    m_mat = new Vecteur[m_dims[1]];
    for(int i = 0; i < m_dims[1]; i++)
    {   
        m_mat[i] = mat.m_mat[i];
    }
}


Matrice::~Matrice()
{
    //std::cout << "Matrice::~Matrice" << std::endl;
    if(m_dims[0] != 0)
    {
        assert(m_dims[1] != 0);
        assert(m_mat != nullptr);

        delete[] m_mat;
    }
    else
    {
        assert(m_dims[1] == 0);
        assert(m_mat == nullptr);
    }
    //std::cout << "Matrice::~Matrice end" << std::endl;
}


Matrice Matrice::operator=(const Matrice &mat)
{
    Matrice tmp(mat);
    swapMats(*this, tmp);

    return *this;
}


Matrice Matrice::operator+(const Matrice &mat) const
{
    assert(m_dims[0] == mat.m_dims[0] && m_dims[1] == mat.m_dims[1]);

    Matrice new_mat(m_dims[0], m_dims[1]);
    for(int i = 0; i < m_dims[1]; i++)
    {
        new_mat[i] = m_mat[i]+mat.m_mat[i];
    }

    return new_mat;
}


Matrice Matrice::operator-(const Matrice &mat) const
{
    assert(m_dims[0] == mat.m_dims[0] && m_dims[1] == mat.m_dims[1]);

    Matrice new_mat(m_dims[0], m_dims[1]);
    for(int i = 0; i < m_dims[1]; i++)
    {
        new_mat[i] = m_mat[i]-mat.m_mat[i];
    }

    return new_mat;
}


Matrice Matrice::operator*(const Matrice &mat) const
{
    //puts("operator*:1");   
    assert(m_dims[1] == mat.m_dims[0]);

    Matrice new_mat(m_dims[0], mat.m_dims[1]);
    for(int i = 0; i < m_dims[0]; i++)
    {
        for(int j = 0; j < mat.m_dims[1]; j++)
        {
            float sum = 0;
            for(int k = 0; k < m_dims[1]; k++)
            {
                sum += m_mat[k][i]*mat.m_mat[j][k];
            }
            new_mat[j][i] = sum;
        }
    }

    //puts("operator*:2"); 
    return new_mat;
}

Vecteur &Matrice::operator[](int index)
{
    assert(index >= 0 && index < m_dims[1]);
    return m_mat[index];
}


const Vecteur &Matrice::operator[](int index) const
{
    assert(index >= 0 && index < m_dims[1]);
    return m_mat[index];
}


Vecteur Matrice::mvprod(const Vecteur &vec) const
{
    assert(vec.getDim() == m_dims[1]);
    Vecteur new_vec(m_dims[0]);

    for(int i = 0; i < m_dims[0]; i++)
    {
        float sum = 0;
        for(int j = 0; j < m_dims[1]; j++)
        {
            sum += m_mat[i][j]*vec[j];
        }
        new_vec[i] = sum;
    }

    return new_vec;
}


Matrice Matrice::transpose() const
{
    //puts("transpose:1"); 

    Matrice new_mat(m_dims[1], m_dims[0]);
    //affiche();
    //new_mat.affiche();

    for(int i = 0; i < m_dims[1]; i++)
    {
        for(int j = 0; j < m_dims[0]; j++)
        {
            new_mat.m_mat[j][i] = m_mat[i][j];
        }
    }

    //puts("transpose:2"); 
    return new_mat;
}


Matrice Matrice::submat(int i1, int i2, int j1, int j2) const
{
    //puts("submat:1");
    TEST_ASSERT(i1, >=, 0);
    TEST_ASSERT(i1, <=, i2);
    TEST_ASSERT(i2, <, m_dims[0]);
    assert(i1 >= 0 && i1 <= i2 && i2 < m_dims[0]);
    assert(j1 >= 0 && j1 <= j2 && j2 < m_dims[1]);

    //affiche();
    Matrice mat(i2-i1+1, j2-j1+1);

    for(int i = i1; i <= i2; i++)
    {
        for(int j = j1; j <= j2; j++)
        {
            //std::cout << "i = " << i << ", j = " << j << std::endl;
            mat.m_mat[j-j1][i-i1] = m_mat[j][i];
        }
    }

    //puts("submat:10");
    return mat;
}


int Matrice::getLigDim() const
{
    return m_dims[0];
}


int Matrice::getColDim() const
{
    return m_dims[1];
}


bool Matrice::isDiagonale() const
{
    assert(m_dims[0] == m_dims[1]);

    for(int i = 0; i < m_dims[0]; i++)
    {
        for(int j = 0; j < m_dims[1]; j++)
        {
            if(i != j && m_mat[j][i] != 0)
            {
                return false;
            }
        }
    }

    return true;
}


bool Matrice::isTridiagonale() const
{
    assert(m_dims[0] == m_dims[1]);

    for(int i = 0; i < m_dims[0]; i++)
    {
        for(int j = 0; j < m_dims[1]; j++)
        {
            if(i != j && i != j+1 && i != j-1 && m_mat[j][i] != 0)
            {
                return false;
            }
        }
    }

    return true;    
}


void Matrice::affiche() const
{
    std::cout << "dims[0] = " << m_dims[0] << "; dims[1] = " << m_dims[1] << std::endl;
    for(int i = 0; i < m_dims[0]; i++)
    {
        for(int j = 0; j < m_dims[1]; j++)
        {
            if(i == 0)
            {
                TEST_ASSERT(m_mat[j].getDim(), ==, m_dims[0]);
            }

            std::cout << m_mat[j][i] << " ";
        }
        std::cout << std::endl;
    }
}


void Matrice::setValues(const Matrice &mat, int i, int j)
{
    assert(i >= 0 && i + mat.getLigDim() <= m_dims[0]);
    assert(j >= 0 && j + mat.getColDim() <= m_dims[1]);

    //puts("1");
    //std::cout << "i = " << i << std::endl;
    //std::cout << "j = " << j << std::endl;    
    for(int k = 0; k < mat.getLigDim(); k++)
    {
        for(int l = 0; l < mat.getColDim(); l++)
        {
            //std::cout << "k = " << k << std::endl;
            //std::cout << "l = " << l << std::endl;  
            m_mat[l+j][k+i] = mat[l][k];
        }
    }
    //puts("5");
}


void Matrice::swap(int index_lig_begin_1, int index_lig_end_1, int index_col_begin_1, int index_col_end_1, int index_lig_begin_2, int index_col_begin_2)
{
    int delta_lig = index_lig_end_1 - index_lig_begin_1;
    int delta_col = index_col_end_1 - index_col_begin_1;
    int index_lig_end_2 = index_lig_begin_2 + delta_lig;
    int index_col_end_2 = index_col_begin_2 + delta_col;
    assert(index_lig_begin_1 >= 0 && index_lig_begin_1 <= index_lig_end_1 && index_lig_end_1 < m_dims[0]);
    assert(index_col_begin_1 >= 0 && index_col_begin_1 <= index_col_end_1 && index_col_end_1 < m_dims[1]);
    assert(index_lig_begin_2 >= 0 && index_lig_begin_2 <= index_lig_end_2 && index_lig_end_2 < m_dims[0]);
    assert(index_col_begin_2 >= 0 && index_col_begin_2 <= index_col_end_2 && index_col_end_2 < m_dims[1]);

    // submatrices do not intersect (by lines or by columns)
    //assert(index_lig_end_1 < index_lig_begin_2 || index_lig_end_2 < index_lig_begin_1 || 
    //    index_col_end_1 < index_col_begin_2 || index_col_end_2 < index_col_begin_1);

    for(int i = 0; i <= delta_lig; i++)
    {
        for(int j = 0; j <= delta_col; j++)
        {
            float tmp = m_mat[j+index_col_begin_1][i+index_lig_begin_1];
            m_mat[j+index_col_begin_1][i+index_lig_begin_1] = m_mat[j+index_col_begin_2][i+index_lig_begin_2];
            m_mat[j+index_col_begin_2][i+index_lig_begin_2] = tmp;
        }
    }
}


void Matrice::zeroLowValues(float epsilon)
{
    for(int i = 0; i < m_dims[0]; i++)
    {
        for(int j = 0; j < m_dims[1]; j++)
        {
            if(abs(m_mat[j][i]) < epsilon)
            {
                m_mat[j][i] = 0;
            }
        }
    }
}


float norm(const Matrice &mat)
{
    float sum = 0;
    for(int i = 0; i < mat.m_dims[0]; i++)
    {
        for(int j = 0; j < mat.m_dims[1]; j++)
        {
            sum += mat.m_mat[j][i]*mat.m_mat[j][i];
        }
    }
    sum = sqrt(sum);
    return sum;
}


Matrice operator*(float numb, const Matrice &mat)
{
    Matrice new_mat(mat);
    for(int i = 0; i < mat.m_dims[0]; i++)
    {
        for(int j = 0; j < mat.m_dims[1]; j++)
        {
            new_mat.m_mat[j][i] *= numb;
        }
        
    }

    return new_mat;    
}


Matrice outer(const Vecteur &vec1, const Vecteur &vec2)
{
    assert(vec1.getDim() == vec2.getDim());
    int size = vec1.getDim();
    Matrice mat(size, size);
    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
        {
            mat.m_mat[i][j] = vec1[i]*vec2[j];
        }
    }

    return mat;
}


void swapMats(Matrice &mat1, Matrice &mat2)
{
    std::swap(mat1.m_mat, mat2.m_mat);
    std::swap(mat1.m_dims[0], mat2.m_dims[0]);
    std::swap(mat1.m_dims[1], mat2.m_dims[1]);
}


Matrice vecToMat(const Vecteur &vec)
{
    Matrice tmp = Matrice(&vec, 1);
    return tmp;
}
