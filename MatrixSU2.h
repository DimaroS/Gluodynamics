#ifndef MATRIXSU2_H
#define MATRIXSU2_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>


#include "main.h"

using namespace std;




template <typename Float, class PrngClass>
class MatrixSU2;

template <typename Float, class PrngClass>
MatrixSU2<Float, PrngClass> operator* (MatrixSU2<Float, PrngClass> A,
                                            MatrixSU2<Float, PrngClass> B);

template <typename Float, class PrngClass>
MatrixSU2<Float, PrngClass> operator/ (MatrixSU2<Float, PrngClass> A,
                                        MatrixSU2<Float, PrngClass> B);

template <typename Float, class PrngClass>
ofstream &operator<< (ofstream &stream, MatrixSU2<Float, PrngClass> A);

template <typename Float, class PrngClass>
ifstream &operator>> (ifstream &stream, MatrixSU2<Float, PrngClass> &A);



template <typename Float, class PrngClass>
class MatrixSU2
{
    public:
        MatrixSU2(int mode, unsigned int seed = 0);
        MatrixSU2() {}
        MatrixSU2 UniformlyRandomMatrix(unsigned int seed);
        MatrixSU2 AroundIdentityMatrix(unsigned int seed, Float epsilon);
        MatrixSU2 HeatBathMatrix(unsigned int seed, Float beta, MatrixSU2<Float, PrngClass> A[6]);


        Float Trace();
        MatrixSU2 Inversed();

        Float& operator[] (int k);
        Float& Get (int k);
        void Set_inacc_step(int k);
        unsigned int Get_inacc_step();
        int GetRepresentationDimension();

        friend MatrixSU2<Float, PrngClass> operator* <Float, PrngClass>
                (MatrixSU2<Float, PrngClass> A, MatrixSU2<Float, PrngClass> B);
        friend MatrixSU2<Float, PrngClass> operator/ <Float, PrngClass>
                (MatrixSU2<Float, PrngClass> A, MatrixSU2<Float, PrngClass> B);
        MatrixSU2 operator*= (MatrixSU2 B);
        MatrixSU2 operator/= (MatrixSU2 B);

        friend ofstream &operator<< <Float, PrngClass> (ofstream &stream,
                                                MatrixSU2<Float, PrngClass> A);
        friend ifstream &operator>> <Float, PrngClass> (ifstream &stream,
                                                MatrixSU2<Float, PrngClass> &A);

    private:
        int inacc_step = 1;
        Float a[4];
};



class MatrixSU2_InitError_WrongMode
{
    public:
        MatrixSU2_InitError_WrongMode(int mode);
        int ret_mode;

};

class MatrixSU2_WrongIndex
{
    public:
        MatrixSU2_WrongIndex(int k);
        int ret_index;
};



template <typename Float, class PrngClass>
ofstream &operator<< (ofstream &stream, MatrixSU2<Float, PrngClass> A) {
    stream.write((char *) &A, sizeof(MatrixSU2<Float, PrngClass>));

    return stream;
}

template <typename Float, class PrngClass>
ifstream &operator>> (ifstream &stream, MatrixSU2<Float, PrngClass> &A) {
    MatrixSU2<Float, PrngClass> tmp;
    stream.read((char *) &tmp, sizeof(MatrixSU2<Float, PrngClass>));
    A = tmp;

    return stream;
}









































template <typename Float, class PrngClass>
int MatrixSU2<Float, PrngClass>::GetRepresentationDimension() {
    return 2;
}


template <typename Float, class PrngClass>
//0 stands for identity matrix
//1 stands for uniformly random matrix
MatrixSU2<Float, PrngClass>::MatrixSU2(int mode, unsigned int seed)
{
    inacc_step = 1;
    if (mode == 0) {
        a[0] = 1.0;
        a[1] = a[2] = a[3] = 0.0;
        return;
    }
    if (mode == 1) {
        *this = MatrixSU2().UniformlyRandomMatrix(seed);
        return;
    }

    throw(MatrixSU2_InitError_WrongMode(mode));
}

MatrixSU2_InitError_WrongMode::MatrixSU2_InitError_WrongMode(int Mode) {
    ret_mode = Mode;
}



template <typename Float, class PrngClass>
Float& MatrixSU2<Float, PrngClass>::operator[] (int k) {
    if (k < 0 || k > 3) {
        throw(MatrixSU2_WrongIndex(k));
    }

    return a[k];
}

template <typename Float, class PrngClass>
Float& MatrixSU2<Float, PrngClass>::Get (int k) {
    if (k < 0 || k > 3) {
        throw(MatrixSU2_WrongIndex(k));
    }

    return a[k];
}

MatrixSU2_WrongIndex::MatrixSU2_WrongIndex(int k) {
    ret_index = k;
}

template <typename Float, class PrngClass>
void MatrixSU2<Float, PrngClass>::Set_inacc_step(int k){
    inacc_step = k;
}

template <typename Float, class PrngClass>
unsigned int MatrixSU2<Float, PrngClass>::Get_inacc_step() {
    return inacc_step;
}


template <typename Float, class PrngClass>
Float MatrixSU2<Float, PrngClass>::Trace() {
    return 2*a[0];
}



template <typename Float, class PrngClass>
MatrixSU2<Float, PrngClass> MatrixSU2<Float, PrngClass>::Inversed() {
    MatrixSU2<Float, PrngClass> temp;
    temp[0] = a[0];
    temp[1] = -a[1];
    temp[2] = -a[2];
    temp[3] = -a[3];

    return temp;
}





template <typename Float, class PrngClass>
MatrixSU2<Float, PrngClass> operator*(MatrixSU2<Float, PrngClass> A,
                                            MatrixSU2<Float, PrngClass> B) {
    MatrixSU2<Float, PrngClass> C;

    C[0] = A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
    C[1] = A[0]*B[1] + B[0]*A[1] - A[2]*B[3] + A[3]*B[2];
    C[2] = A[0]*B[2] + B[0]*A[2] - A[3]*B[1] + A[1]*B[3];
    C[3] = A[0]*B[3] + B[0]*A[3] - A[1]*B[2] + A[2]*B[1];

    C.Set_inacc_step(A.Get_inacc_step() +  B.Get_inacc_step());
    if (C.Get_inacc_step() > 10000) {
        Float r = sqrt(C[0]*C[0] + C[1]*C[1] + C[2]*C[2] + C[3]*C[3]);
        for (int i = 0; i < 4; i++) {
            C[i] /= r;
        }
        C.Set_inacc_step(1);
    }

    return C;
}

template <typename Float, class PrngClass>
MatrixSU2<Float, PrngClass> operator/ (MatrixSU2<Float, PrngClass> A,
                                            MatrixSU2<Float, PrngClass> B) {
    return A*B.Inversed();
}

template <typename Float, class PrngClass>
MatrixSU2<Float, PrngClass> MatrixSU2<Float, PrngClass>::operator*=
                                            (MatrixSU2<Float, PrngClass> B) {
    MatrixSU2<Float, PrngClass> temp = (*this)*B;
    *this = temp;

    return temp;
}

template <typename Float, class PrngClass>
MatrixSU2<Float, PrngClass> MatrixSU2<Float, PrngClass>::operator/=
                                            (MatrixSU2<Float, PrngClass> B) {
    MatrixSU2<Float, PrngClass> temp = (*this)/B;
    *this = temp;

    return *this;
}






template <typename Float, class PrngClass>
MatrixSU2<Float, PrngClass> MatrixSU2<Float, PrngClass>::UniformlyRandomMatrix
                                                        (unsigned int _seed) {
    Float f[4];
    Float rr;
    PrngClass rand_gen;
    rand_gen.seed(_seed);
    while (true) {
        rr = 0;
        for (int i = 0; i < 4; i++) {
            f[i] = ((Float)rand_gen())/(rand_gen.max() - rand_gen.min());
            rr += f[i]*f[i];
        }
        if (rr == 0 || rr > 1.0) {
            continue;
        }

        Float r = sqrt(rr);
        for (int i = 0; i < 4; i++) {
            a[i] = f[i]/r;
            if (rand_gen()%2) {
                a[i] = -a[i];
            }
        }
        inacc_step = 1;
        break;
    }

    return *this;
}



class MatrixSU2_AroundIdentityMatrix_WRONG_EPSILON
{
    public:
        float epsilon;
        MatrixSU2_AroundIdentityMatrix_WRONG_EPSILON() {}
};




template <typename Float, class PrngClass>
MatrixSU2<Float, PrngClass> MatrixSU2<Float, PrngClass>::AroundIdentityMatrix
                                        (unsigned int _seed, Float epsilon) {
    if (epsilon < 0.0 || epsilon > 1.0) {
        throw(MatrixSU2_AroundIdentityMatrix_WRONG_EPSILON());
    }

//    epsilon /= 1.0;

    Float f[4];
    PrngClass rand_gen;
    rand_gen.seed(_seed);
    while (true) {
        for (int i = 1; i < 4; i++) {
            f[i] = ((Float)rand_gen())/(rand_gen.max() - rand_gen.min()) - 0.5;
        }

        Float rr = 0.0;
        for (int i = 1; i < 4; i++) {
            rr += f[i]*f[i];
        }
        if (rr == 0.0) {
            continue;
        }
        Float r = sqrt(rr);
        for (int i = 1; i < 4; i++) {
            a[i] = epsilon * f[i] / r;
        }
        a[0] = sqrt(1.0 - epsilon*epsilon);
        if (rand_gen()%2) {
            a[0] = -a[0];
        }
        break;
    }



    inacc_step = 1;
    return *this;
}





template <typename Float, class PrngClass>
MatrixSU2<Float, PrngClass> MatrixSU2<Float, PrngClass>::HeatBathMatrix
            (unsigned int _seed, Float beta, MatrixSU2<Float, PrngClass> A[6]) {
    Float __a[4] = {0.0, 0.0, 0.0, 0.0};
    //cout << '!';
    for (int i = 0; i < 6; i++) {
        __a[0] += A[i][0];
        __a[1] += A[i][1];
        __a[2] += A[i][2];
        __a[3] += A[i][3];
    }
    Float _a_a = __a[0]*__a[0] + __a[1]*__a[1]
                                            + __a[2]*__a[2] + __a[3]*__a[3];
    if (_a_a == 0.0) {
        cout << '!';
        *this = UniformlyRandomMatrix(_seed);
        return *this;
    }
    Float _a = sqrt(_a_a);

    MatrixSU2<Float, PrngClass> V(0);
    for (int i = 0; i < 4; i++) {
        V[i] = - __a[i] / _a;
    }

    PrngClass rand_gen;
    rand_gen.seed(_seed);

    MatrixSU2<Float, PrngClass> X(0);
    Float lambda_sqr, r[4];
    while (true) {
        for (int i = 1; i < 4; i++) {
            r[i] = 0.0;
            while (r[i] == 0.0) {
                r[i] = ((Float)rand_gen())/(rand_gen.max() - rand_gen.min());
            }
        }
        lambda_sqr = -(log(r[1]) + cos(2*M_PI*r[2])*cos(2*M_PI*r[2])*log(r[3]))
                                                            / (2.0 * _a * beta);
        Float _rand = ((Float)rand_gen())/(rand_gen.max() - rand_gen.min());
        if (_rand*_rand < 1.0 - lambda_sqr) {
            break;
        }
    }
    X[0] = 1.0 - 2.0*lambda_sqr;

    while (true) {
        for (int i = 1; i < 4; i++) {
            r[i] = 0.0;
            while (r[i] == 0.0) {
                r[i] = 2*((Float)rand_gen())/(rand_gen.max() - rand_gen.min())
                                                                        - 1.0;
            }
        }
        if (r[1]*r[1] + r[2]*r[2] + r[3]*r[3] < 1.0) {
            break;
        }
    }

    Float rr = r[1]*r[1] + r[2]*r[2] + r[3]*r[3];
    for (int i = 1; i < 4; i++) {
        X[i] = r[i] * sqrt((1.0 - X[0]*X[0]) / rr);
    }
//    cout << X[0] << ' ' << X[1] << ' ' << X[2] << ' ' << X[3] << '\n';
//    cout << V[0]*V[0] + V[1]*V[1] + V[2]*V[2] + V[3]*V[3] <<
//            ' ' << X[0]*X[0] + X[1]*X[1] + X[2]*X[2] + X[3]*X[3]
//                                        << ' ' << V[0] << ' ' << _a;

//    cout << X[0] << '*' << beta << '+';

    *this = V.Inversed()*X;

//    (*this).UniformlyRandomMatrix(_seed);
    return (*this);
}





#endif // MATRIXSU2_H
