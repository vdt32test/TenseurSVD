#include "TenseurSVD.hpp"


void vec_test();
void mat_test();
void svd_test();
void ten_test();
void tensvd_test();


constexpr float sqrtOf2 = 1.41421356237;
constexpr float sqrtOf3 = 1.73205080757;


int main(int argc, char *argv[])
{
    assert(argc == 2);

    int cmd = atoi(argv[1]);

    switch(cmd)
    {
    case(1):
        vec_test();
        break;

    case(2):
        mat_test();
        break;

    case(3):
        svd_test();
        break;

    case(4):    
        ten_test();
        break;

    case(5):
        tensvd_test();
        break;

    default:
        assert(false);
    }

    return 0;
}


void vec_test()
{
    float tab1[3] = {1.0, 1.0, 1.0};
    float tab2[4] = {3.0, 4.0, 0.0, 0.0};
    Vecteur U(tab1, 3), V(tab2, 4);

    std::cout << "1.";
    U.affiche();
    V.affiche();

    Vecteur t(U);
    
    float tab3[3] = {1.0, 1.0, 0.0};
    U = Vecteur(tab3, 3);

    std::cout << "\n3.";
    U.affiche();
    t.affiche();

    float v_dot = dot(V, V);
    float v_norm = norm(V);
    std::cout << "\n4.V^T*V = " << v_dot << "\n|V|=" << v_norm << std::endl;

    Vecteur tmp;
    tmp = (1/v_norm)*V;
    std::cout << "\n5.";
    tmp.affiche();

    Vecteur W;
    W = V.subvec(1,3);
    std::cout << "\n6.";
    V.affiche();
    W.affiche();

    std::cout << "\n7.";
    tmp = U+W;
    tmp.affiche();
    tmp = U-W;
    tmp.affiche();
}


void mat_test()
{
    float Aarr1[3] = {1.0, 1.0, 0.0};
    float Aarr2[3] = {-0.5, 2.0, -1.0};
    float Aarr3[3] = {0.0, -1.0, 1.0};
    Vecteur Avec[3] = {Vecteur(Aarr1, 3), Vecteur(Aarr2, 3), Vecteur(Aarr3, 3)};
    Matrice A(Avec, 3);

    float Barr1[2] = {-2.0, 0.0};
    float Barr2[2] = {3.0, 1.0};
    Vecteur Bvec[2] = {Vecteur(Barr1, 2), Vecteur(Barr2,2 )};
    Matrice B(Bvec, 2);

    std::cout << "1.";
    A.affiche();
    B.affiche();

    Matrice C(B);
    B[1][0] = 0;

    std::cout << "\n2.";
    B.affiche();
    C.affiche();

    Matrice D(3,2);
    D = A.submat(0,2,0,1);

    std::cout << "\n3.";
    D.affiche();

    float varr[3] = {3,2,1};
    Matrice E(Vecteur(varr, 3));

    std::cout << "\n4.";
    E.affiche();

    std::cout << "\n5.";
    (B+C).affiche();
    (B-C).affiche();
    (D*C).affiche();

    std::cout << "\n6." << norm(C);

    std::cout << "\n7.";
    B[1][0] = 3;
    (0.5*(B+B.transpose())).affiche();
}


void svd_test()
{
    float c, s;
    givens(1, 2, c, s);
    std::cout << "1.c = " << c << ", s =  " << s << std::endl;

    
    float xFl[2] = { -1, 0 };
    float yFl[2] = { 1/sqrtOf2, 1/sqrtOf2 };
    float zFl[1] = { -4 };
    Vecteur x(xFl, 2), y(yFl, 2), z(zFl, 1);
    
    float beta;
    Vecteur v = householder(x, beta);
    std::cout << "\n2.1)beta = " << beta << std::endl;
    v.affiche();

    v = householder(y, beta);
    std::cout << "\n2)beta = " << beta << std::endl;
    v.affiche();

    std::cout << "\n3)";
    v = householder(z, beta);
    std::cout << "beta = " << beta << std::endl;
    v.affiche();    
    
    std::cout << "\n3.";
    float MFl[2][2] = { {10, -6}, 
                        {-6, 10} };
    Vecteur MVec[2] = { Vecteur(MFl[0], 2), Vecteur(MFl[1], 2) };
    Matrice M(MVec, 2);
    Matrice Q = reductridiag(M);

    Q.affiche();

    std::cout << "\n4.";
    float M2Fl[3][3] = { {     1, 1.0/2, 1.0/3 }, 
                         { 1.0/2, 1.0/3, 1.0/4 },
                         { 1.0/3, 1.0/4, 1.0/5 } };

    Vecteur M2Vec[3] = { Vecteur(M2Fl[0], 3), Vecteur(M2Fl[1], 3), Vecteur(M2Fl[2], 3) };
    M = Matrice(M2Vec, 3);
    
    Matrice R;
    Matrice Pi = qrpivot(M, Q, &R);
    std::cout << "R:" << std::endl;
    R.affiche();
    std::cout << "Q:" << std::endl;
    Q.affiche();
    std::cout << "Pi:" << std::endl;
    Pi.affiche();

    std::cout << "\n5.";
    float AFl[2][2] = { {1, 0 },
                        {0, -1}};
    Vecteur AVec[2] = { Vecteur(AFl[0], 2), Vecteur(AFl[1], 2) };
    Matrice A(AVec, 2);

    float BFl[2][2] = { {  2*sqrtOf2, -sqrtOf2 },
                        { -2*sqrtOf2, -sqrtOf2 } };
    Vecteur BVec[2] = { Vecteur(BFl[0], 2), Vecteur(BFl[1], 2) };
    Matrice B(BVec, 2);
    
    float CFl[3][2] = { {       1.0/2, sqrtOf3/2 },
                        { 3*sqrtOf3/2,    -3.0/2 },
                        {           0,         0 } };
    Vecteur CVec[3] = { Vecteur(CFl[0], 2), Vecteur(CFl[1], 2), Vecteur(CFl[2], 2) };
    Matrice C(CVec, 3);

    Matrice svdMatrices[3];
    Matrice &U = svdMatrices[0];
    Matrice &Sigma = svdMatrices[1];
    Matrice &V = svdMatrices[2];
    
    svd(A, svdMatrices);
    std::cout << "1)U:" << std::endl;
    U.affiche();
    std::cout << "Sigma:" << std::endl;
    Sigma.affiche();
    std::cout << "V:" << std::endl;
    V.affiche();

    std::cout << "A:" << std::endl;
    A.affiche();
    std::cout << "U*Sigma*V^T:" << std::endl;
    (U*Sigma*(V.transpose())).affiche();

    svd(B, svdMatrices);
    std::cout << "2)U:" << std::endl;
    U.affiche();
    std::cout << "Sigma:" << std::endl;
    Sigma.affiche();
    std::cout << "V:" << std::endl;
    V.affiche();

    std::cout << "B:" << std::endl;
    B.affiche();
    std::cout << "U*Sigma*V^T:" << std::endl;
    (U*Sigma*(V.transpose())).affiche();

    svd(C, svdMatrices);
    std::cout << "3)U:" << std::endl;
    U.affiche();
    std::cout << "Sigma:" << std::endl;
    Sigma.affiche();
    std::cout << "V:" << std::endl;
    V.affiche();

    std::cout << "C:" << std::endl;
    C.affiche();
    std::cout << "U*Sigma*V^T:" << std::endl;
    (U*Sigma*(V.transpose())).affiche();

    svd(M, svdMatrices);
    std::cout << "4)U:" << std::endl;
    U.affiche();
    std::cout << "Sigma:" << std::endl;
    Sigma.affiche();
    std::cout << "V:" << std::endl;
    V.affiche();

    std::cout << "M:" << std::endl;
    M.affiche();
    std::cout << "U*Sigma*V^T:" << std::endl;
    (U*Sigma*(V.transpose())).affiche();
}


void ten_test()
{
    int dims1[3] = { 2, 2, 2 };
    Tenseur Tau(dims1, 3);
    std::cout << "1.";
    Tau.affiche();

    Vecteur UVec(8, 1);
    Tenseur U(dims1, 3, UVec);
    std::cout << "\n2.";
    U.affiche();

    Tenseur V(U+Tau), W(U-Tau);
    std::cout << "\n3.";
    V.affiche();
    W.affiche();

    std::cout << "\n4.U[2][2][2] = " << std::endl;
    std::cout << U[dims1] << std::endl;
    U[dims1] = -1;
    U.affiche();
    V.affiche();

    std::cout << "\n5.";
    Tau.mode(1).affiche();

    std::cout << "\n6.";
    Tau[0] = 1;
    Tau[1] = 4;
    Tau[2] = 3;
    Tau[3] = 1.0/3;
    Tau[4] = 0;
    Tau[5] = 1.5;
    Tau[6] = -1;
    Tau[7] = 2;

    Tau.mode(2).affiche();

    std::cout << "\n7.";
    float AFl[2][3] = { {  3, 0,  0 }, 
                        { -1, 6, -3 } };
    Vecteur AVec[2] = { Vecteur(AFl[0], 3), Vecteur(AFl[1], 3) };
    Matrice A(AVec, 2);
    //A.affiche();

    Tenseur S = pmod(Tau, A, 3);
    S.affiche();

    std::cout << "\n8.";
    Tenseur R = S + S;
    R.affiche();
}


void tensvd_test()
{                         
    float TauFl[27] = { 0.9073,  0.7158, -0.3698, 1.7842,  1.6970, 0.0151,  2.1236, -0.0740,  1.4429, 
                        0.8924, -0.4898,  2.4288, 1.7753, -1.5077, 4.0337, -0.6631,  1.9103, -1.7495,
                        2.1488,  0.3054,  2.3753, 4.2495,  0.3207, 4.7146,  1.8260,  2.1335, -0.2716 };

    int dims[3] = { 3, 3, 3 };
    
    Tenseur Tau(dims, 3, Vecteur(TauFl, 27));
    Tau.affiche();

    Tau.mode(1).affiche();
    Tau.mode(2).affiche();
    Tau.mode(3).affiche();

    TenseurSVD TauSVD = hosvd(Tau);
    TauSVD.affiche();

    TauSVD.getS().affiche();
    TauSVD.showErrorMatrices();
}
