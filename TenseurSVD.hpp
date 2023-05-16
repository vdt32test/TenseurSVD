#pragma once
#include "Tenseur.hpp"
#include "SVD.hpp"


class TenseurSVD : private Tenseur
{
public:
    TenseurSVD();
    TenseurSVD(const int *dims, int d, const Matrice *facteurs);
    TenseurSVD(const int *dims, int d, const Vecteur &vecTen, const Matrice *facteurs);
    TenseurSVD(const int *dims, int d, int k, const Matrice &matTen, const Matrice *facteurs);
    TenseurSVD(const Tenseur &ten, const Matrice *facteurs);
    TenseurSVD(const TenseurSVD &ten);
    ~TenseurSVD();

    TenseurSVD operator=(const TenseurSVD &ten);

    Tenseur tenseurtotal() const;
    Tenseur getS() const;

    void affiche() const;
    void showErrorMatrices() const;
    friend void swapTensSVD(TenseurSVD &ten1, TenseurSVD &ten2);

private:
    Matrice *m_facteurs;
};


TenseurSVD hosvd(const Tenseur &ten);