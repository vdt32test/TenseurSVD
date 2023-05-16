#include "Tenseur.hpp"


Tenseur::Tenseur():
    m_d(0),
    m_dims(nullptr),
    m_nbelts(0),
    m_ten()
{}


Tenseur::Tenseur(const int *dims, int d):
    m_d(d),
    m_dims(nullptr),
    m_nbelts(0),
    m_ten()
{
    initDims(dims, d);
    m_ten = Vecteur(m_nbelts);
}


Tenseur::Tenseur(const int *dims, int d, const Vecteur &vecTen):
    m_d(d),
    m_dims(nullptr),
    m_nbelts(0),
    m_ten(vecTen)
{
    initDims(dims, d);
    assert(m_nbelts == m_ten.getDim());
}


Tenseur::Tenseur(const int *dims, int d, int k, const Matrice &matTen):
    m_d(d),
    m_dims(nullptr),
    m_nbelts(0),
    m_ten()
{
    assert(k >= 1 && k <= d);
    initDims(dims, d);
    assert(m_nbelts == matTen.getLigDim()*matTen.getColDim());
    m_ten = Vecteur(m_nbelts);
    
    assert(matTen.getLigDim() == m_dims[k-1]);
    assert(matTen.getColDim() == m_nbelts/m_dims[k-1] && m_nbelts%m_dims[k-1] == 0);
    
    for(int j = 0; j < matTen.getColDim(); j++)
    {
        for(int i = 0; i < matTen.getLigDim(); i++)
        {
            m_ten[i + matTen.getLigDim()*j] = matTen[j][i];
        }
    }
}


Tenseur::Tenseur(const Tenseur &ten):
    m_d(ten.m_d),
    m_dims(nullptr),
    m_nbelts(0),
    m_ten(ten.m_ten)
{
    initDims(ten.m_dims, m_d);
    assert(m_nbelts == m_ten.getDim());
}


Tenseur::~Tenseur()
{
    if(m_d != 0)
    {
        assert(m_dims != nullptr && m_nbelts != 0 && m_ten.getDim() != 0);
        delete[] m_dims;
    }
    else
    {
        assert(m_dims == nullptr && m_nbelts == 0 && m_ten.getDim() == 0);
    }
}


Tenseur Tenseur::operator=(const Tenseur &ten)
{
    Tenseur tmp(ten);
    swapTens(*this, tmp);

    assert(m_nbelts == m_ten.getDim());

    return *this;
}


float Tenseur::operator[](int i) const
{
    assert(i >= 0 && i < m_nbelts && m_nbelts == m_ten.getDim());
    return m_ten[i];
}


float &Tenseur::operator[](int i)
{
    assert(i >= 0 && i < m_nbelts && m_nbelts == m_ten.getDim());
    return m_ten[i];
}


float Tenseur::operator[](int *indexes) const
{
    int i = tenToVecIndex(m_dims, m_d, indexes);
    assert(i >= 0 && i < m_nbelts && m_nbelts == m_ten.getDim());
    return m_ten[i];
}


float &Tenseur::operator[](int *indexes)
{
    int i = tenToVecIndex(m_dims, m_d, indexes);
    assert(i >= 0 && i < m_nbelts && m_nbelts == m_ten.getDim());
    return m_ten[i];
}


Tenseur Tenseur::operator+(const Tenseur &ten) const
{
    assert(m_d == ten.m_d && m_nbelts == ten.m_nbelts);
    for(int i = 0; i < m_d; i++)
    {
        assert(m_dims[i] == ten.m_dims[i]);
    }

    Tenseur new_ten(m_dims, m_d);
    for(int i = 0; i < m_nbelts; i++)
    {
        new_ten.m_ten[i] = m_ten[i]+ten.m_ten[i];
    }

    return new_ten;
}


Tenseur Tenseur::operator-(const Tenseur &ten) const
{
    assert(m_d == ten.m_d && m_nbelts == ten.m_nbelts);
    for(int i = 0; i < m_d; i++)
    {
        assert(m_dims[i] == ten.m_dims[i]);
    }

    Tenseur new_ten(m_dims, m_d);
    for(int i = 0; i < m_nbelts; i++)
    {
        new_ten.m_ten[i] = m_ten[i]-ten.m_ten[i];
    }

    return new_ten;
}


Matrice Tenseur::mode(int k) const
{
    assert(k >= 1 && k <= m_d);
    int m = m_dims[k-1]; // lig dim
    int n = m_nbelts/m; // col dim
    int *tenIndexes = new int[m_d];
    
    assert(n*m == m_nbelts);

    Matrice mat(m, n);
    for(int i = 0; i < m_nbelts; i++)
    {
        vecToTenIndex(m_dims, m_d, i, tenIndexes);
        int matIndex1 = tenIndexes[k-1]-1;
        int matIndex2 = tenToMatIndex(m_dims, m_d, tenIndexes, k);

        mat[matIndex2][matIndex1] = m_ten[i];
    }

    delete[] tenIndexes;

    return mat;
}


int Tenseur::getD() const
{
    return m_d;
}


void Tenseur::initDims(const int *dims, int d)
{
    m_d = d;
    m_dims = new int[d];
    m_nbelts = dims[0];

    for(int i = 0; i < d; i++)
    {
        m_dims[i] = dims[i];
        if(i != 0)
        {
            m_nbelts *= m_dims[i];
        }
    }
    assert(m_nbelts != 0);
}


void Tenseur::affiche() const
{
    std::cout << "m_d = " << m_d << std::endl;
    std::cout << "m_dims = ";
    for(int i = 0; i < m_d-1; i++)
    {
        std::cout << m_dims[i] << "x";
    }
    std::cout << m_dims[m_d-1];

    std::cout << std::endl;
    std::cout << "m_nbelts = " << m_nbelts << std::endl;
    std::cout << "m_ten: " << std::endl;
    m_ten.affiche();
}


void Tenseur::deinit()
{
    if(m_d != 0)
    {
        assert(m_dims != nullptr && m_nbelts != 0 && m_ten.getDim() != 0);
        delete[] m_dims;
        m_nbelts = 0;
        m_d = 0;
        m_ten = Vecteur();
    }
    else
    {
        assert(m_dims == nullptr && m_nbelts == 0 && m_ten.getDim() == 0);
    }
}


Tenseur pmod(const Tenseur &S, const Matrice &M, int k)
{
    assert(k >= 1 && k <= S.m_d);
    assert(S.m_dims[k-1] == M.getColDim());
    int nk = S.m_dims[k-1];

    int *new_dims = new int[S.m_d];

    for(int i = 0; i < S.m_d; i++)
    {
        new_dims[i] = S.m_dims[i];

    }
    new_dims[k-1] = M.getLigDim();

    int new_nbelts = S.m_nbelts/nk*M.getLigDim();

    int new_nbelts_test = new_dims[0];
    for(int i = 1; i < S.m_d; i++)
    {
        new_nbelts_test *= new_dims[i];
    }

    assert(new_nbelts == new_nbelts_test);

    Tenseur new_ten(new_dims, S.m_d);
    TEST_ASSERT(new_ten.m_ten.getDim(), ==, new_nbelts);
    delete[] new_dims;

    int *indexes = new int[S.m_d];

    for(int i = 0; i < new_nbelts; i++)
    {
        vecToTenIndex(new_ten.m_dims, S.m_d, i, indexes);
        int kIndex = indexes[k-1]-1;
        for(int j = 0; j < nk; j++)
        {
            indexes[k-1] = j+1;
            int s_index = tenToVecIndex(S.m_dims, S.m_d, indexes);
            new_ten[i] += M[j][kIndex]*S[s_index];
        }
    }

    delete[] indexes;
    
    return new_ten;
}


void swapTens(Tenseur &ten1, Tenseur &ten2)
{
    std::swap(ten1.m_d, ten2.m_d);
    std::swap(ten1.m_dims, ten2.m_dims);
    std::swap(ten1.m_nbelts, ten2.m_nbelts);
    swap(ten1.m_ten, ten2.m_ten);
}


int tenToVecIndex(const int *dims, int d, const int *tenIndexes)
{
    int vecIndex = tenIndexes[0]-1;
    for(int i  = 1; i < d; i++)
    {
        assert(tenIndexes[i] > 0 && tenIndexes[i] <= dims[i]);
        vecIndex = vecIndex * dims[i] + tenIndexes[i] - 1;
    }

    return vecIndex;
}


void vecToTenIndex(const int *dims, int d, int vecIndex, int *tenIndexes)
{
    int ft = vecIndex + 1;
    tenIndexes[d-1] = indexDivision(ft, dims[d-1]);
    for(int i = d-2; i >= 0; i--)
    {
        TEST_ASSERT((ft - tenIndexes[i+1]) % dims[i+1], ==, 0);
        ft = (ft - tenIndexes[i+1]) / dims[i+1] + 1;
        tenIndexes[i] = indexDivision(ft, dims[i]);
    }

    for(int i = 0; i < d; i++)
    {
        assert(tenIndexes[i] > 0);
        assert(tenIndexes[i] <= dims[i]);
    }
}


int tenToMatIndex(const int *dims, int d, const int *tenIndexes, int k)
{
    int *matIndexes = new int[d-1];
    int *matDims = new int[d-1];
    for(int i = 0; i < k-1; i++)
    {
        matDims[i] = dims[i];
        matIndexes[i] = tenIndexes[i];
    }

    for(int i = k; i < d; i++)
    {
        matDims[i-1] = dims[i];
        matIndexes[i-1] = tenIndexes[i];
    }

    int matIndex = tenToVecIndex(matDims, d-1, matIndexes);

    delete[] matIndexes;
    delete[] matDims;

    return matIndex;
}


int indexDivision(int a, int b)
{
    assert(b != 0);

    int ret = a % b;
    if(ret == 0)
    {
        ret = b;
    }

    return ret;
}
