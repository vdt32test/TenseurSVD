#include "TenseurSVD.hpp"


TenseurSVD::TenseurSVD():
    Tenseur(),
    m_facteurs(nullptr)
{}


TenseurSVD::TenseurSVD(const int *dims, int d, const Matrice *facteurs):
    Tenseur(dims, d),
    m_facteurs(nullptr)
{
    m_facteurs = new Matrice[d];
    for(int i = 0; i < d; i++)
    {
        m_facteurs[i] = facteurs[i];
    }
}


TenseurSVD::TenseurSVD(const int *dims, int d, const Vecteur &vecTen, const Matrice *facteurs):
    Tenseur(dims, d, vecTen),
    m_facteurs(nullptr)
{
    m_facteurs = new Matrice[d];
    for(int i = 0; i < d; i++)
    {
        m_facteurs[i] = facteurs[i];
    }
}


TenseurSVD::TenseurSVD(const int *dims, int d, int k, const Matrice &matTen, const Matrice *facteurs):
    Tenseur(dims, d, k, matTen),
    m_facteurs(nullptr)
{
    m_facteurs = new Matrice[d];
    for(int i = 0; i < d; i++)
    {
        m_facteurs[i] = facteurs[i];
    }
}


TenseurSVD::TenseurSVD(const Tenseur &ten, const Matrice *facteurs):
    Tenseur(ten),
    m_facteurs(nullptr)
{
    m_facteurs = new Matrice[m_d];
    for(int i = 0; i < m_d; i++)
    {
        m_facteurs[i] = facteurs[i];
    }
}


TenseurSVD::TenseurSVD(const TenseurSVD &ten):
    Tenseur(ten),
    m_facteurs(nullptr)
{
    m_facteurs = new Matrice[m_d];
    for(int i = 0; i < m_d; i++)
    {
        m_facteurs[i] = ten.m_facteurs[i];
    }
}


TenseurSVD::~TenseurSVD()
{
    if(m_d != 0)
    {
        assert(m_facteurs != nullptr);
        delete[] m_facteurs;
    }
    else
    {
        assert(m_facteurs == nullptr);
    }
}


TenseurSVD TenseurSVD::operator=(const TenseurSVD &ten)
{
    TenseurSVD tmp(ten);
    swapTensSVD(*this, tmp);

    return *this;
}


Tenseur TenseurSVD::tenseurtotal() const
{
    return *this;
}


Tenseur TenseurSVD::getS() const
{
    Tenseur S(*this);

    for(int i = 0; i < m_d; i++)
    {
        S = pmod(S, m_facteurs[i].transpose(), i+1);
    }

    Tenseur Tau(S);
    for(int i = 0; i < m_d; i++)
    {
        Tau = pmod(Tau, m_facteurs[i], i+1);
    }

    std::cout << "SxU1xU2x...xUn:" << std::endl;
    Tau.affiche();

    return S;
}


void TenseurSVD::affiche() const
{
    Tenseur::affiche();
    for(int i = 0; i < m_d; i++)
    {
        std::cout << "U[" << i << "]:" << std::endl;
        m_facteurs[i].affiche();
    }
}


void TenseurSVD::showErrorMatrices() const
{
    int d = m_d;
    Matrice *Us = new Matrice[d];
    Matrice *Sigmas = new Matrice[d];
    Matrice *Vs = new Matrice[d];
    Matrice *Ms = new Matrice[d];
    Matrice mats[3];

    for(int i = 0; i < d; i++)
    {
        Matrice Tk = mode(i+1);
        svd(Tk, mats);
        Us[i] = mats[0];
        Sigmas[i] = mats[1];
        Vs[i] = mats[2];
        Matrice TkHat = Us[i] * Sigmas[i] * (Vs[i].transpose());
        Matrice M = (1/norm(Tk))*(Tk - TkHat);

        std::cout << "M (k = " << i+1 << "):" << std::endl;
        M.affiche();          
    }

    delete[] Us;
    delete[] Sigmas;
    delete[] Vs;
    delete[] Ms;
}


void swapTensSVD(TenseurSVD &ten1, TenseurSVD &ten2)
{
    swapTens(ten1, ten2);
    std::swap(ten1.m_facteurs, ten2.m_facteurs);
}


TenseurSVD hosvd(const Tenseur &ten)
{
    int d = ten.getD();
    Matrice *Us = new Matrice[d];
    Matrice mats[3];
    for(int i = 0; i < d; i++)
    {
        Matrice Tk = ten.mode(i+1);
        svd(Tk, mats);
        Us[i] = mats[0];
        std::cout << "Us[i]:" << std::endl;
        Us[i].affiche(); 
        std::cout << "Sigma:" << std::endl;
        mats[1].affiche();        
        std::cout << "V:" << std::endl;
        mats[2].affiche();
        std::cout << "U*Sigma*V^T:" << std::endl;
        (mats[0]*mats[1]*(mats[2].transpose())).affiche();          
    }

    TenseurSVD tenSVD(ten, Us);
    delete[] Us;
    return tenSVD;
}
