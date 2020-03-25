#ifndef MATRIXSU3_H
#define MATRIXSU3_H


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>


#include "main.h"

using namespace std;


template <typename Float, class PrngClass>
class MatrixSU3;

template <typename Float, class PrngClass>
MatrixSU3<Float, PrngClass> operator*
                (MatrixSU3<Float, PrngClass> A, MatrixSU3<Float, PrngClass> B);

template <typename Float, class PrngClass>
MatrixSU3<Float, PrngClass> operator/
                (MatrixSU3<Float, PrngClass> A, MatrixSU3<Float, PrngClass> B);

template <typename Float, class PrngClass>
ofstream &operator<< (ofstream &stream, MatrixSU3<Float, PrngClass> A);

template <typename Float, class PrngClass>
ifstream &operator>> (ifstream &stream, MatrixSU3<Float, PrngClass> &A);



template <typename Float, class PrngClass>
class MatrixSU3
{
    public:
        MatrixSU3(int mode, unsigned int seed = 0);
        MatrixSU3();
        MatrixSU3 UniformlyRandomMatrix(unsigned int seed);
        MatrixSU3 AroundIdentityMatrix(unsigned int seed, Float epsilon);



        Float Trace();
        MatrixSU3 Inversed();

        Float& GetRe (int k, int m);
        Float& GetIm (int k, int m);
        void Set_inacc_step(int k);
        unsigned int Get_inacc_step();
        int GetRepresentationDimension();

        friend MatrixSU3<Float, PrngClass> operator* <Float, PrngClass>
                (MatrixSU3<Float, PrngClass> A, MatrixSU3<Float, PrngClass> B);
        friend MatrixSU3<Float, PrngClass> operator/ <Float, PrngClass>
                (MatrixSU3<Float, PrngClass> A, MatrixSU3<Float, PrngClass> B);
        MatrixSU3 operator*= (MatrixSU3 B);
        MatrixSU3 operator/= (MatrixSU3 B);
        void Unitarization();

        friend ofstream &operator<< <Float, PrngClass>
                            (ofstream &stream, MatrixSU3<Float, PrngClass> A);
        friend ifstream &operator>> <Float, PrngClass>
                            (ifstream &stream, MatrixSU3<Float, PrngClass> &A);

    private:
        int inacc_step = 0;
        Float re[3][3], im[3][3];
};



class MatrixSU3_InitError_WrongMode
{
    public:
        MatrixSU3_InitError_WrongMode(int mode);
        int ret_mode;
};


class MatrixSU3_WrongIndex
{
    public:
        MatrixSU3_WrongIndex() {}
        int ret_index;
};




template <typename Float, class PrngClass>
ofstream &operator<< (ofstream &stream, MatrixSU3<Float, PrngClass> A) {
    stream.write((char *) &A, sizeof(MatrixSU3<Float, PrngClass>));

    return stream;
}

template <typename Float, class PrngClass>
ifstream &operator>> (ifstream &stream, MatrixSU3<Float, PrngClass> &A) {
    MatrixSU3<Float, PrngClass> tmp;
    stream.read((char *) &tmp, sizeof(MatrixSU3<Float, PrngClass>));
    A = tmp;

    return stream;
}




















template <typename Float, class PrngClass>
MatrixSU3<Float, PrngClass>::MatrixSU3() {
    return;
}

template <typename Float, class PrngClass>
int MatrixSU3<Float, PrngClass>::GetRepresentationDimension() {
    return 3;
}


template <typename Float, class PrngClass>
//0 stands for identity matrix
//1 stands for uniformly random matrix
MatrixSU3<Float, PrngClass>::MatrixSU3(int mode, unsigned int seed)
{
    inacc_step = 1;
    if (mode == 0) {
        im[0][0] = im[0][1] = im[0][2] = im[1][0] = im[1][1] = im[1][2]
                                    = im[2][0] = im[2][1] = im[2][2] = 0.0;
        re[0][0] = re[1][1] = re[2][2] = 1.0;
        re[0][1] = re[0][2] = re[1][0] = re[1][2] = re[2][0] = re[2][1] = 0.0;
        return;
    }
    if (mode == 1) {
        *this = MatrixSU3().UniformlyRandomMatrix(seed);
        return;
    }

    throw(MatrixSU3_InitError_WrongMode(mode));
}

MatrixSU3_InitError_WrongMode::MatrixSU3_InitError_WrongMode(int Mode) {
    ret_mode = Mode;
}


template <typename Float, class PrngClass>
Float &MatrixSU3<Float, PrngClass>::GetRe (int k, int m) {
    if (k < 0 || k > 2 || m < 0 || m > 2) {
        throw(MatrixSU3_WrongIndex());
    }

    return re[k][m];
}

template <typename Float, class PrngClass>
Float &MatrixSU3<Float, PrngClass>::GetIm (int k, int m) {
    if (k < 0 || k > 2 || m < 0 || m > 2) {
        throw(MatrixSU3_WrongIndex());
    }

    return im[k][m];
}

template <typename Float, class PrngClass>
void MatrixSU3<Float, PrngClass>::Set_inacc_step(int k) {
    inacc_step = k;
}

template <typename Float, class PrngClass>
unsigned int MatrixSU3<Float, PrngClass>::Get_inacc_step() {
    return inacc_step;
}


template <typename Float, class PrngClass>
Float MatrixSU3<Float, PrngClass>::Trace() {
    return re[0][0] + re[1][1] + re[2][2];
}



template <typename Float, class PrngClass>
MatrixSU3<Float, PrngClass> MatrixSU3<Float, PrngClass>::Inversed() {
    MatrixSU3<Float, PrngClass> temp;
    temp.GetRe(0, 0) = re[0][0];
    temp.GetRe(0, 1) = re[1][0];
    temp.GetRe(0, 2) = re[2][0];
    temp.GetRe(1, 0) = re[0][1];
    temp.GetRe(1, 1) = re[1][1];
    temp.GetRe(1, 2) = re[2][1];
    temp.GetRe(2, 0) = re[0][2];
    temp.GetRe(2, 1) = re[1][2];
    temp.GetRe(2, 2) = re[2][2];

    temp.GetIm(0, 0) = -im[0][0];
    temp.GetIm(0, 1) = -im[1][0];
    temp.GetIm(0, 2) = -im[2][0];
    temp.GetIm(1, 0) = -im[0][1];
    temp.GetIm(1, 1) = -im[1][1];
    temp.GetIm(1, 2) = -im[2][1];
    temp.GetIm(2, 0) = -im[0][2];
    temp.GetIm(2, 1) = -im[1][2];
    temp.GetIm(2, 2) = -im[2][2];

    return temp;
}



template <typename Float, class PrngClass>
void MatrixSU3<Float, PrngClass>::Unitarization() {
    MatrixSU3<Float, PrngClass> C = *this;
    Float ru = sqrt(C.GetRe(0, 0)*C.GetRe(0, 0) + C.GetRe(0, 1)*C.GetRe(0, 1)
                    + C.GetRe(0, 2)*C.GetRe(0, 2) + C.GetIm(0, 0)*C.GetIm(0, 0)
                    + C.GetIm(0, 1)*C.GetIm(0, 1) + C.GetIm(0, 2)*C.GetIm(0, 2));
    for (int k = 0; k < 3; k++) {
        C.GetRe(0, k) /= ru;
        C.GetIm(0, k) /= ru;
    }

    Float uvRe = C.GetRe(0, 0)*C.GetRe(1, 0) + C.GetRe(0, 1)*C.GetRe(1, 1)
                + C.GetRe(0, 2)*C.GetRe(1, 2) + C.GetIm(0, 0)*C.GetIm(1, 0)
                + C.GetIm(0, 1)*C.GetIm(1, 1) + C.GetIm(0, 2)*C.GetIm(1, 2);
    Float uvIm = C.GetIm(0, 0)*C.GetRe(1, 0) + C.GetIm(0, 1)*C.GetRe(1, 1)
                + C.GetIm(0, 2)*C.GetRe(1, 2) - C.GetRe(0, 0)*C.GetIm(1, 0)
                - C.GetRe(0, 1)*C.GetIm(1, 1) - C.GetRe(0, 2)*C.GetIm(1, 2);
    for (int k = 0; k < 3; k++) {
        C.GetRe(1, k) -= uvRe*C.GetRe(0, k) + uvIm*C.GetIm(0, k);
        C.GetIm(1, k) -= uvRe*C.GetIm(0, k) - uvIm*C.GetRe(0, k);
    }
    Float rv = sqrt(C.GetRe(1, 0)*C.GetRe(1, 0) + C.GetRe(1, 1)*C.GetRe(1, 1)
                    + C.GetRe(1, 2)*C.GetRe(1, 2) + C.GetIm(1, 0)*C.GetIm(1, 0)
                    + C.GetIm(1, 1)*C.GetIm(1, 1) + C.GetIm(1, 2)*C.GetIm(1, 2));
    for (int k = 0; k < 3; k++) {
        C.GetRe(1, k) /= rv;
        C.GetIm(1, k) /= rv;
    }

    Float upRe = C.GetRe(0, 0)*C.GetRe(2, 0) + C.GetRe(0, 1)*C.GetRe(2, 1)
                    + C.GetRe(0, 2)*C.GetRe(2, 2) + C.GetIm(0, 0)*C.GetIm(2, 0)
                    + C.GetIm(0, 1)*C.GetIm(2, 1) + C.GetIm(0, 2)*C.GetIm(2, 2);
    Float upIm = C.GetIm(0, 0)*C.GetRe(2, 0) + C.GetIm(0, 1)*C.GetRe(2, 1)
                    + C.GetIm(0, 2)*C.GetRe(2, 2) - C.GetRe(0, 0)*C.GetIm(2, 0)
                    - C.GetRe(0, 1)*C.GetIm(2, 1) - C.GetRe(0, 2)*C.GetIm(2, 2);
    for (int k = 0; k < 3; k++) {
        C.GetRe(2, k) -= upRe*C.GetRe(0, k) + upIm*C.GetIm(0, k);
        C.GetIm(2, k) -= upRe*C.GetIm(0, k) - upIm*C.GetRe(0, k);
    }
    Float vpRe = C.GetRe(1, 0)*C.GetRe(2, 0) + C.GetRe(1, 1)*C.GetRe(2, 1)
                    + C.GetRe(1, 2)*C.GetRe(2, 2) + C.GetIm(1, 0)*C.GetIm(2, 0)
                    + C.GetIm(1, 1)*C.GetIm(2, 1) + C.GetIm(1, 2)*C.GetIm(2, 2);
    Float vpIm = C.GetIm(1, 0)*C.GetRe(2, 0) + C.GetIm(1, 1)*C.GetRe(2, 1)
                    + C.GetIm(1, 2)*C.GetRe(2, 2) - C.GetRe(1, 0)*C.GetIm(2, 0)
                    - C.GetRe(1, 1)*C.GetIm(2, 1) - C.GetRe(1, 2)*C.GetIm(2, 2);
    for (int k = 0; k < 3; k++) {
        C.GetRe(2, k) -= vpRe*C.GetRe(1, k) + vpIm*C.GetIm(1, k);
        C.GetIm(2, k) -= vpRe*C.GetIm(1, k) - vpIm*C.GetRe(1, k);
    }
    Float rp = sqrt(C.GetRe(2, 0)*C.GetRe(2, 0) + C.GetRe(2, 1)*C.GetRe(2, 1)
                    + C.GetRe(2, 2)*C.GetRe(2, 2) + C.GetIm(2, 0)*C.GetIm(2, 0)
                    + C.GetIm(2, 1)*C.GetIm(2, 1) + C.GetIm(2, 2)*C.GetIm(2, 2));
    for (int k = 0; k < 3; k++) {
        C.GetRe(2, k) /= rp;
        C.GetIm(2, k) /= rp;
    }



    C.Set_inacc_step(1);
    *this = C;
}



template <typename Float, class PrngClass>
MatrixSU3<Float, PrngClass> operator*
                (MatrixSU3<Float, PrngClass> A, MatrixSU3<Float, PrngClass> B) {
    MatrixSU3<Float, PrngClass> C;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            C.GetRe(i, j) = C.GetIm(i, j) = 0.0;

            for (int k = 0; k < 3; k++) {
                C.GetRe(i, j) += A.GetRe(i, k)*B.GetRe(k, j)
                                                - A.GetIm(i, k)*B.GetIm(k, j);
                C.GetIm(i, j) += A.GetRe(i, k)*B.GetIm(k, j)
                                                + A.GetIm(i, k)*B.GetRe(k, j);
            }
        }
    }


    C.Set_inacc_step(A.Get_inacc_step() +  B.Get_inacc_step());
    if (C.Get_inacc_step() > 100) {
        C.Unitarization();
    }

    return C;
}

template <typename Float, class PrngClass>
MatrixSU3<Float, PrngClass> operator/
            (MatrixSU3<Float, PrngClass> A, MatrixSU3<Float, PrngClass> B) {
    return A*B.Inversed();
}

template <typename Float, class PrngClass>
MatrixSU3<Float, PrngClass> MatrixSU3<Float, PrngClass>::operator*=
                                        (MatrixSU3<Float, PrngClass> B) {
    MatrixSU3<Float, PrngClass> temp = (*this)*B;
    *this = temp;

    return temp;
}

template <typename Float, class PrngClass>
MatrixSU3<Float, PrngClass> MatrixSU3<Float, PrngClass>::operator/=
                                    (MatrixSU3<Float, PrngClass> B) {
    MatrixSU3<Float, PrngClass> temp = (*this)/B;
    *this = temp;

    return *this;
}






template <typename Float, class PrngClass>
MatrixSU3<Float, PrngClass> MatrixSU3<Float, PrngClass>::UniformlyRandomMatrix
                                                        (unsigned int _seed) {
    Float f4[4];
    Float rr;
    PrngClass rand_gen;
    rand_gen.seed(_seed);
    while (true) {
        rr = 0;
        for (int i = 0; i < 4; i++) {
            f4[i] = (Float)rand_gen()/(rand_gen.max() - rand_gen.min());
            rr += f4[i]*f4[i];
        }
        if (rr == 0 || rr > 1.0) {
            continue;
        }

        Float r = sqrt(rr);
        for (int i = 0; i < 4; i++) {
            f4[i] /= r;
            if (rand_gen()%2) {
                f4[i] = -f4[i];
            }
        }
        break;
    }


    Float f6[6];
    while (true) {
        rr = 0;
        for (int i = 0; i < 6; i++) {
            f6[i] = (Float)rand_gen()/(rand_gen.max() - rand_gen.min());
            rr += f6[i]*f6[i];
        }
        if (f6[0] == 0.0 || f6[1] == 0.0 || f6[2] == 0.0 || f6[3] == 0.0 ||
                                    f6[4] == 0.0 || f6[5] == 0.0 || rr > 1.0) {
            continue;
        }

        Float r = sqrt(rr);
        for (int i = 0; i < 6; i++) {
            f6[i] /= r;
            if (rand_gen()%2) {
                f6[i] = -f6[i];
            }
        }
        break;
    }

    MatrixSU3<Float, PrngClass> R, S;

    R.GetRe(2, 0) = R.GetRe(2, 1) = R.GetRe(0, 2) = R.GetRe(1, 2) = 0.0;
    R.GetRe(2, 2) = 1.0;
    R.GetIm(2, 0) = R.GetIm(2, 1) = R.GetIm(2, 2)
                                    = R.GetIm(0, 2) = R.GetIm(1, 2) = 0.0;

    R.GetRe(0, 0) = f4[0];
    R.GetRe(1, 1) = f4[0];
    R.GetRe(0, 1) = f4[2];
    R.GetRe(1, 0) = -f4[2];
    R.GetIm(0, 0) = f4[3];
    R.GetIm(1, 1) = -f4[3];
    R.GetIm(0, 1) = f4[1];
    R.GetIm(1, 0) = f4[1];


    Float sin_theta = sqrt(f6[0]*f6[0] + f6[1]*f6[1]);
    Float cos_theta = sqrt(1.0 - sin_theta*sin_theta);
    Float ei23re = f6[0] / sin_theta;
    Float ei23im = f6[1] / sin_theta;
    Float sin_phi = sqrt(f6[2]*f6[2] + f6[3]*f6[3])/cos_theta;
    Float cos_phi = sqrt(1.0 - sin_phi*sin_phi);
    Float ei22re = f6[2]/cos_theta/sin_phi;
    Float ei22im = f6[3]/cos_theta/sin_phi;
    Float ei21re = f6[4]/cos_theta/cos_phi;
    Float ei21im = f6[5]/cos_theta/cos_phi;

    S.GetRe(0, 0) = ei23re * cos_theta;
    S.GetRe(0, 1) = 0.0;
    S.GetRe(0, 2) = f6[0];
    S.GetRe(1, 0) = -ei22re * sin_theta * sin_phi;
    S.GetRe(1, 1) = cos_phi * (ei21re*ei23re - ei21im*ei23im);
    S.GetRe(1, 2) = f6[2];
    S.GetRe(2, 0) = -ei21re * sin_theta * cos_phi;
    S.GetRe(2, 1) = -sin_phi * (ei22re*ei23re - ei22im*ei23im);
    S.GetRe(2, 2) = f6[4];
    S.GetIm(0, 0) = ei23im * cos_theta;
    S.GetIm(0, 1) = 0.0;
    S.GetIm(0, 2) = f6[1];
    S.GetIm(1, 0) = -ei22im * sin_theta * sin_phi;
    S.GetIm(1, 1) = (-ei21re*ei23im - ei21im*ei23re) * cos_phi;
    S.GetIm(1, 2) = f6[3];
    S.GetIm(2, 0) = -ei21im * sin_theta * cos_phi;
    S.GetIm(2, 1) = -(-ei22re*ei23im - ei22im*ei23re) * sin_phi;
    S.GetIm(2, 2) = f6[5];


    *this = S*R;
    inacc_step = 1;
/*
    Float rr = 0.0;
    Float f6[6];
    while (true) {
        rr = 0.0;
        for (int i = 0; i < 6; i++) {
            f6[i] = (Float)rand()/RAND_MAX;
            rr += f6[i]*f6[i];
        }
        if (f6[0] == 0.0 || f6[1] == 0.0 || f6[2] == 0.0 || f6[3] == 0.0 || f6[4] == 0.0 || f6[5] == 0.0 || rr > 1.0) {
            continue;
        }

        Float r = sqrt(rr);
        for (int i = 0; i < 6; i++) {
            f6[i] /= r;
            if (rand()%2) {
                f6[i] = -f6[i];
            }
        }
        break;
    }


    Float g6[6];
    while (true) {
        rr = 0.0;
        for (int i = 0; i < 6; i++) {
            g6[i] = (Float)rand()/RAND_MAX;
            rr += g6[i]*g6[i];
        }
        if (g6[0] == 0.0 || g6[1] == 0.0 || g6[2] == 0.0 || g6[3] == 0.0 || g6[4] == 0.0 || g6[5] == 0.0 || rr > 1.0) {
            continue;
        }

        Float r = sqrt(rr);
        for (int i = 0; i < 6; i++) {
            g6[i] /= r;
            if (rand()%2) {
                g6[i] = -g6[i];
            }
        }
        break;
    }


    Float h6[6];
    while (true) {
        rr = 0.0;
        for (int i = 0; i < 6; i++) {
            h6[i] = (Float)rand()/RAND_MAX;
            rr += h6[i]*h6[i];
        }
        if (h6[0] == 0.0 || h6[1] == 0.0 || h6[2] == 0.0 || h6[3] == 0.0 || h6[4] == 0.0 || h6[5] == 0.0 || rr > 1.0) {
            continue;
        }

        Float r = sqrt(rr);
        for (int i = 0; i < 6; i++) {
            h6[i] /= r;
            if (rand()%2) {
                h6[i] = -h6[i];
            }
        }
        break;
    }


    for (int i = 0; i < 3; i++) {
        re[0][i] = f6[2*i];
        im[0][i] = f6[2*i + 1];
        re[1][i] = g6[2*i];
        im[1][i] = g6[2*i +1];
        re[2][i] = h6[2*i];
        im[2][i] = h6[2*i + 1];
    }

    Unitarization();

*/

    return *this;
}



class MatrixSU3_AroundIdentityMatrix_WRONGEPSILON
{
    public:
        float epsilon;
        MatrixSU3_AroundIdentityMatrix_WRONGEPSILON() {epsilon = 2.0;}
};



template <typename Float, class PrngClass>
MatrixSU3<Float, PrngClass> MatrixSU3<Float, PrngClass>::AroundIdentityMatrix
                                        (unsigned int _seed, Float epsilon) {
    if (epsilon < 0.0 || epsilon > 1.0) {
        throw(MatrixSU3_AroundIdentityMatrix_WRONGEPSILON());
    }

//    epsilon = epsilon/2.0;
    PrngClass rand_gen;
    rand_gen.seed(_seed);


    Float r4[4];
    while (true) {
        for (int i = 1; i < 4; i++) {
            r4[i] = ((Float)rand_gen())/(rand_gen.max() - rand_gen.min()) - 0.5;
        }

        Float rr = 0.0;
        for (int i = 1; i < 4; i++) {
            rr += r4[i]*r4[i];
        }
        if (rr == 0.0) {
            continue;
        }
        Float r = sqrt(rr);
        for (int i = 1; i < 4; i++) {
            r4[i] = epsilon * r4[i] / r;
        }
        r4[0] = sqrt(1 - epsilon*epsilon);
        if (rand_gen()%2) {
            r4[0] = -r4[0];
        }
        break;
    }

    Float s4[4];
    while (true) {
        for (int i = 1; i < 4; i++) {
            s4[i] = ((Float)rand_gen())/(rand_gen.max() - rand_gen.min()) - 0.5;
        }

        Float rr = 0.0;
        for (int i = 1; i < 4; i++) {
            rr += s4[i]*s4[i];
        }
        if (rr == 0.0) {
            continue;
        }
        Float r = sqrt(rr);
        for (int i = 1; i < 4; i++) {
            s4[i] = epsilon * s4[i] / r;
        }
        s4[0] = sqrt(1 - epsilon*epsilon);
        if (rand_gen()%2) {
            s4[0] = -s4[0];
        }
        break;
    }

    Float t4[4];
    while (true) {
        for (int i = 1; i < 4; i++) {
            t4[i] = ((Float)rand_gen())/(rand_gen.max() - rand_gen.min()) - 0.5;
        }

        Float rr = 0.0;
        for (int i = 1; i < 4; i++) {
            rr += t4[i]*t4[i];
        }
        if (rr == 0.0) {
            continue;
        }
        Float r = sqrt(rr);
        for (int i = 1; i < 4; i++) {
            t4[i] = epsilon * t4[i] / r;
        }
        t4[0] = sqrt(1 - epsilon*epsilon);
        if (rand_gen()%2) {
            t4[0] = -t4[0];
        }
        break;
    }




    MatrixSU3<Float, PrngClass> R, S, T;

    R.GetRe(2, 0) = R.GetRe(2, 1) = R.GetRe(0, 2) = R.GetRe(1, 2) = 0.0;
    R.GetRe(2, 2) = 1.0;
    R.GetIm(2, 0) = R.GetIm(2, 1) = R.GetIm(2, 2)
                                    = R.GetIm(0, 2) = R.GetIm(1, 2) = 0.0;

    R.GetRe(0, 0) = r4[0];
    R.GetRe(1, 1) = r4[0];
    R.GetRe(0, 1) = r4[2];
    R.GetRe(1, 0) = -r4[2];
    R.GetIm(0, 0) = r4[3];
    R.GetIm(1, 1) = -r4[3];
    R.GetIm(0, 1) = r4[1];
    R.GetIm(1, 0) = r4[1];

    S.GetRe(1, 0) = S.GetRe(1, 2) = S.GetRe(0, 1) = S.GetRe(2, 1) = 0.0;
    S.GetRe(1, 1) = 1.0;
    S.GetIm(1, 0) = S.GetIm(2, 1) = S.GetIm(1, 1)
                                    = S.GetIm(0, 1) = S.GetIm(1, 2) = 0.0;

    S.GetRe(0, 0) = s4[0];
    S.GetRe(2, 2) = s4[0];
    S.GetRe(0, 2) = s4[2];
    S.GetRe(2, 0) = -s4[2];
    S.GetIm(0, 0) = s4[3];
    S.GetIm(2, 2) = -s4[3];
    S.GetIm(0, 2) = s4[1];
    S.GetIm(2, 0) = s4[1];

    T.GetRe(0, 1) = T.GetRe(0, 2) = T.GetRe(1, 0) = T.GetRe(2, 0) = 0.0;
    T.GetRe(0, 0) = 1.0;
    T.GetIm(0, 1) = T.GetIm(0, 2) = T.GetIm(0, 0)
                                    = T.GetIm(1, 0) = T.GetIm(2, 0) = 0.0;

    T.GetRe(1, 1) = t4[0];
    T.GetRe(2, 2) = t4[0];
    T.GetRe(1, 2) = t4[2];
    T.GetRe(2, 1) = -t4[2];
    T.GetIm(1, 1) = t4[3];
    T.GetIm(2, 2) = -t4[3];
    T.GetIm(1, 2) = t4[1];
    T.GetIm(2, 1) = t4[1];

    switch ((rand_gen() - rand_gen.min())%6) {
        case 0:
            *this = R*S*T;
        break;

        case 1:
            *this = R*T*S;
        break;

        case 2:
            *this = S*R*T;
        break;

        case 3:
            *this = S*T*R;
        break;

        case 4:
            *this = T*R*S;
        break;

        case 5:
            *this = T*S*R;
        break;
    }

    inacc_step = 1;
    return *this;
}









#endif // MATRIXSU3_H
