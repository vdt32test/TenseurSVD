#include "SVD.hpp"


int sign(float x)
{
    if(x < 0) 
    {
        return -1;
    }

    return 1;
}


void givens(float x, float z, float &c, float &s)
{
    if(z == 0)
    {
        c = 1;
        s = 0;
        return;
    }

    float tau;
    if(abs(z) > abs(x))
    {
        tau = -x/z;
        s = 1/sqrt(1+tau*tau);
        c = s*tau;
    }
    else
    {
        tau = -z/x;
        c = 1/sqrt(1+tau*tau);
        s = c*tau;
    }
}


Vecteur householder(const Vecteur &x, float &beta)
{
    assert(x.getDim() > 0);

    float sigma;
    if(x.getDim() > 1)
    {
        Vecteur tmp = x.subvec(1, x.getDim()-1);
        sigma = dot(tmp, tmp);
    }
    else
    {
        sigma = 0;
    }
    
    Vecteur V(x);
    V[0] = 1;

    if(sigma == 0)
    {
        if(x[0] >= 0)
        {
            beta = 0;
        }
        else
        {
            beta = 2;
        }
    }
    else
    {
        float mu = sqrt(x[0]*x[0]+sigma);
        if(x[0] <= 0)
        {
            V[0] = x[0] - mu;
        }
        else
        {
            V[0] = - sigma/(x[0] + mu);
        }
        beta = 2*V[0]*V[0]/(sigma+V[0]*V[0]);
        V = (1.0/V[0])*V;
    }

    return V;
}


Matrice reductridiag(Matrice &D)
{
    assert(D.getColDim() == D.getLigDim());
    int DSize = D.getColDim();

    float d = (D[DSize-2][DSize-2] - D[DSize-1][DSize-1])/2;
    float mu = D[DSize-1][DSize-1] - D[DSize-2][DSize-1]*D[DSize-2][DSize-1]/
        (d+sign(d)*sqrt(d*d+D[DSize-2][DSize-1]*D[DSize-2][DSize-1]));

    float x = D[0][0] - mu;
    float z = D[0][1];

    Matrice Z(Vecteur(DSize, 1));
    for(int k = 0; k < DSize-1; k++)
    {
        float c, s;
        givens(x, z, c, s);

        for(int j = 0; j < DSize; j++)
        {
            float t1 = D[k][j];
            float t2 = D[k+1][j];
            D[k][j]=c*t1 - s*t2;
            D[k+1][j]=s*t1 + c*t2;
            t1=Z[k][j];
            t2=Z[k+1][j];
            Z[k][j]=c*t1 - s*t2;
            Z[k+1][j]=s*t1 + c*t2;
        }

        for(int j = 0; j < DSize; j++)
        {
            float t1 = D[j][k];
            float t2 = D[j][k+1];
            D[j][k] = c*t1 - s*t2;
            D[j][k+1] = s*t1 + c*t2;
        }

        if(k < DSize-2)
        {
            x = D[k][k+1];
            z = D[k][k+2];
        }
    }

    D.zeroLowValues();
    Z.zeroLowValues();

    return Z;
}


void qrsym(const Matrice &Ain, Matrice &Q)
{
    Matrice A(Ain);

    assert(A.getColDim() == A.getLigDim());
    int n = A.getLigDim();
    Q = Matrice(Vecteur(n, 1));

    for(int k = 0; k < n-2; k++)
    {
        float beta;
        Matrice tmpMat = A.submat(k+1, n-1, k, k);
        Vecteur tmpVec = tmpMat[0];

        Vecteur v = householder(tmpVec, beta);
        Matrice tmp1 = A.submat(k+1, n-1, k+1, n-1);
        Matrice tmp2 = tmp1*Matrice(&v,1);
        assert(tmp2.getColDim() == 1);
        
        Vecteur p = beta*tmp2[0];

        Vecteur w = p - (beta/2)*dot(p, v)*v;
        A[k][k+1] = norm(A.submat(k+1,n-1,k,k));
        A[k+1][k] = A[k][k+1];

        Matrice vMat(&v, 1), wMat(&w, 1);
        Matrice tmp(A.submat(k+1, n-1, k+1, n-1));
        tmp = tmp - vMat*(wMat.transpose()) - wMat*(vMat.transpose());

        for(int i = k+1; i < n; i++)
        {
            for(int j = k+1; j < n; j++)
            {
                A[j][i] = tmp[j-k-1][i-k-1];
            }
        }

        tmp = Q.submat(k+1, n-1, k+1, n-1);
        tmp = tmp - beta*vMat*(vMat.transpose())*tmp;

        for(int i = k+1; i < n; i++)
        {
            for(int j = k+1; j < n; j++)
            {
                Q[j][i] = tmp[j-k-1][i-k-1];
            }
        }    
    }

    Matrice T(n);
    for(int j = n-1; j >= 0; j--)
    {
        T[j][j] = A[j][j];
        if(j > 0)
        {
            T[j][j-1] = A[j-1][j];
            T[j-1][j] = T[j][j-1];
        }
    }

    while(!T.isDiagonale())
    {
        for(int i = 0; i < n-1; i++)
        {
            if(abs(T[i+1][i])+abs(T[i][i+1]) <= 0.000000001*(abs(T[i][i])+abs(T[i+1][i+1])))
            {
                T[i+1][i] = 0;
                T[i][i+1] = 0;
            }
        }

        if(!T.isTridiagonale())
        {
            T.zeroLowValues();
        }

        if(!T.isTridiagonale())
        {
            T.affiche();
        }
        
        assert(T.isTridiagonale());
        
        int p = 0;
        for(int i = 0; i < n-2; i++)
        {
            Matrice tmpMat(T.submat(0,i+1, 0,i+1));
            if(tmpMat.isDiagonale())
            {
                p = i+1;
            }
        }

        
        int q = 0;
        for(int i = n-1; i >= 2; i--)
        {
            Matrice tmpMat(T.submat(i-1, n-1, i-1, n-1));
            if(tmpMat.isDiagonale())
            {
                q = n-i;
            }
        }

        assert(p + q <= n);
        Matrice T1;
        if(p != 0)
        {
            T1 = T.submat(0, p-1, 0, p-1);
        }

        Matrice T3;
        if(q != 0)
        {
            T3 = T.submat(n-q, n-1, n-q, n-1);
        }

        assert(T1.isDiagonale());
        assert(T3.isDiagonale());
        assert(T1.getColDim() == p);
        assert(T3.getColDim() == q);

        Matrice T2;
        T2 = T.submat(p, n-q-1, p, n-q-1);

        assert(T2.isTridiagonale()); 
        assert(T1.getColDim()+T2.getColDim()+T3.getColDim() == T.getColDim());

        if(p + q < n)
        {
            assert(T2.getColDim() > 0);
            Matrice Z = reductridiag(T2);

            assert(T2.getColDim() == n - p - q);
            Matrice THat(T);
            for(int i = 0; i < T2.getColDim(); i++)
            {
                for(int j = 0; j < T2.getLigDim(); j++)
                {
                    THat[j+p][i+p] = T2[j][i];
                }
            }

            T = (1.0/2)*(THat + THat.transpose());

            Matrice tmp(Vecteur(n,1));

            for(int i = 0; i < Z.getColDim(); i++)
            {
                for(int j = 0; j < Z.getLigDim(); j++)
                {
                    tmp[j+p][i+p] = Z[j][i];
                }
            }

            Q = Q*tmp;    
        }

    }
}


Matrice qrpivot(const Matrice &Ain, Matrice &Q, Matrice *AOut)
{
    // m, n - sizes
    // i, j, k, r - indexes 
    Matrice A(Ain); 

    int m = A.getLigDim();
    int n = A.getColDim();
    assert(m >= n);

    Matrice Pi(Vecteur(n,1));
    Vecteur c(n);

    for(int j = 0; j < n; j++)
    {
        Matrice tmp(A.submat(0, m-1, j, j));
        tmp = tmp.transpose()*tmp;
        assert(tmp.getColDim() == 1 && tmp.getLigDim() == 1);
        c[j] = tmp[0][0];
    }

    int r = -1;
    float tau = c.getMax();
    while(tau > 0 && r < n-1)
    {
        r++;
        int k;
        bool found = false;

        for(int i = r; i < n; i++)
        {
            if(c[i] == tau)
            {
                k = i;
                found = true;
                break;
            }
        }
        assert(found);

        A.swap(0, m-1, r, r, 0, k);

        float tmp = c[r];
        c[r] = c[k];
        c[k] = tmp;

        Pi.swap(0, n-1, r, r, 0, k);

        float beta;
        assert(A.submat(r, m-1, r, r).getColDim() == 1);
        Vecteur v = householder(A.submat(r, m-1, r, r)[0], beta);

        Matrice tmpMat = A.submat(r, m-1, r, n-1) - beta*(Matrice(&v,1)*(Matrice(&v,1).transpose()))*A.submat(r, m-1, r, n-1);
        A.setValues(tmpMat, r, r);
        if(v.getDim() > 1)
        {
            A.setValues(vecToMat(v.subvec(1, m-r-1)), r+1, r);
        }

        for(int i = r+1; i < n; i++)
        {
            c[i] -= A[i][r]*A[i][r];
        }

        if(r < n-1)
        {
            tau = c.subvec(r+1, n-1).getMax();
        }
        else
        {
            tau = 0;
        }
    }

    Q = Matrice(Vecteur(m, 1));
    Vecteur v(m);

    for(int j = n-1; j >= 0; j--)
    {
        v[j] = 1;
        float normVal;

        if(j == n-1)
        {
            normVal = 0;
        }
        else
        {
            Matrice tmpMat = A.submat(j+1, m-1, j, j);
            v.setValues(tmpMat[0], j+1);
            normVal = norm(tmpMat);
        }
        
        float beta = 2/(1 + normVal*normVal);
        
        Vecteur tmpVec(v.subvec(j, m-1));
        Matrice tmpVecMat(&tmpVec, 1);
        Matrice tmpMat = tmpVecMat*(tmpVecMat.transpose());
        tmpMat = tmpMat*Q.submat(j, m-1, j, m-1);
        tmpMat = beta*tmpMat;
        tmpMat = Q.submat(j, m-1, j, m-1) - tmpMat;
        Q.setValues(tmpMat, j, j);
    }

    if(AOut != nullptr)
    {
        *AOut = A;
    }
    return Pi;
}


void svd(const Matrice &A, Matrice *svdMatrices)
{
    int m = A.getLigDim();
    int n = A.getColDim();
    Matrice &U = svdMatrices[0];
    Matrice &Sigma = svdMatrices[1];
    Matrice &V = svdMatrices[2];
    Matrice Q1, Q2, Pi, R;
    
    if(m >= n)
    {
        qrsym(A.transpose()*A, Q1);
        Pi = qrpivot(A*Q1, Q2);

        R = Q2.transpose()*(A*Q1)*Pi;
        for(int j = 0; j < n; j++)
        {
            if(R[j][j] < 0)
            {
                Q2.setValues((-1)*Q2.submat(0, m-1, j, j), 0, j);
            }
        }
        R = Q2.transpose()*(A*Q1)*Pi;
        U = Q2;
        Sigma = R;
        V = Q1*Pi;
    }
    else
    {
        qrsym(A*(A.transpose()), Q1);
        Pi = qrpivot(A.transpose()*Q1, Q2);
        R = Q2.transpose()*(A.transpose()*Q1)*Pi;
        for(int i = 0; i < m; i++)
        {
            if(R[i][i] < 0)
            {
                Q2.setValues((-1)*Q2.submat(0, n-1, i, i), 0, i);
            }
        }
        R = Q2.transpose()*(A.transpose()*Q1)*Pi;
        U = Q1*Pi;
        Sigma = R.transpose();
        V = Q2;
    }

    U.zeroLowValues();
    Sigma.zeroLowValues();
    V.zeroLowValues();
}