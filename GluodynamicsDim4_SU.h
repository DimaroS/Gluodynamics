#ifndef GLUODYNAMICSDIM4_SU_H
#define GLUODYNAMICSDIM4_SU_H


#include <cmath>
#include <cstdlib>
#include <utility>
#include <exception>

using namespace std;







template <class Link>
class LinkNodeDim4
{
private:
    Link m[4];
public:
    LinkNodeDim4() = default;

    Link& up (int k) {
        if (k < 1 || k > 4) {
            throw invalid_argument("LinkNodeDim4_DimentionError");
        }

        return m[k-1];
    }

    Link copy_of_up (int k) const {
        if (k < 1 || k > 4) {
            throw invalid_argument("LinkNodeDim4_DimentionError");
        }

        return m[k-1];
    }
};























template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class DynamicArray,
                typename Float, class PrngClass>
class GluodynamicsWithoutBorderDim4_SU_Base
{
    public:
        GluodynamicsWithoutBorderDim4_SU_Base() = default;
        GluodynamicsWithoutBorderDim4_SU_Base(int N1, int N2, int N3, int N4,
                                              int mode, float g0_input, unsigned int seed);

        DynamicArray<LinkNodeDim4<MatrixSU<Float, PrngClass>>> m;

        virtual MatrixSU<Float, PrngClass> GetLink(int x, int y, int z, int t,
                                                            int direction) const = 0;
        MatrixSU<Float, PrngClass> * PrecalculatedPlaquettes
                                        (int i, int j, int k, int l, int direction) const ;


        bool MetropolisHit(int i, int j, int k, int l, int direction,
                                            MatrixSU<Float, PrngClass> *Usim);
        unsigned int SmallMonteCarloStep(unsigned int number_of_hits);
        float MonteCarloStep(float number_of_steps,
                                    unsigned int number_of_hits);
        unsigned int HeatBath_SmallMonteCarloStep
            (unsigned int number_of_hits, int i, int j, int k, int l, int direction);
        float HeatBath_MonteCarloStep(float number_of_steps,
                                    unsigned int number_of_hits);

        unsigned int Timeslice_I_SmallMonteCarloStep
                (unsigned int number_of_hits, int mean_timeslice, int width);
        unsigned int Timeslice_II_SmallMonteCarloStep
                (unsigned int number_of_hits, int mean_timeslice, int width);
        float Timeslice_I_MonteCarloStep(float number_of_steps,
                unsigned int number_of_hits, int mean_timeslice, int width);
        float Timeslice_II_MonteCarloStep(float number_of_steps,
                unsigned int number_of_hits, int mean_timeslice, int width);

        void Base_Input(ifstream &stream);
        void Base_Output(ofstream &stream) const;

        double Action();
        double AverageWilsonLoop(int I, int J) const;
        double AverageOrientedWilsonLoop(int I, int J, int i_direction,
                                                            int j_direction) const;
        double SingleWilsonLoop(int I, int J, int i_direction, int j_direction,
                int i, int j, int k, int l) const;
        double SingleKnot(int a_direction, int b_direction, int c_direction, int d_direction,
                            int i, int j, int k, int l) const;
        double AverageFace() const;

        float Get_g0() const;
        float Get_beta() const;
        float Get_t() const;
        void Set_g0(float g0_new);
        void Set_beta(float beta_new);
        void Set_t(float t_new);
        void Seed(unsigned int key);



    protected:
        PrngClass rand_gen;
        int arr_size_1;
        int arr_size_2;
        int arr_size_3;
        int arr_size_4;
        float g0, beta;
        float t = 0.0;
        double action = 0.0;
        bool action_measured = false;
};













template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
float GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                                                            ::Get_g0() const {
    return g0;
}

template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
float GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                                                            ::Get_beta() const {
    return beta;
}

template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
float GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                                                            ::Get_t() const {
    return t;
}

template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
void GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                                                        ::Set_g0(float g0_new) {
    g0 = g0_new;
    beta = 2.0*MatrixSU<Float, PrngClass>().GetRepresentationDimension()/g0/g0;
    action_measured = false;
}

template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
void GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                                                    ::Set_beta(float beta_new) {
    beta = beta_new;
    g0 = sqrt(2.0*MatrixSU<Float, PrngClass>().GetRepresentationDimension()
                                                                        /beta);
    action_measured = false;
}

template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
void GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                                                    ::Set_t(float t_new) {
    t = t_new;
}

template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
void GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                                                    ::Seed(unsigned int key) {
    rand_gen.seed(key);
}



template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class DynamicArray,
                typename Float, class PrngClass>
void GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, DynamicArray, Float, PrngClass>::
                                        Base_Output (ofstream &stream) const {
    stream.write((const char *) &arr_size_1, sizeof(int));
    stream.write((const char *) &arr_size_2, sizeof(int));
    stream.write((const char *) &arr_size_3, sizeof(int));
    stream.write((const char *) &arr_size_4, sizeof(int));

    stream.write((const char *) &g0, sizeof(float));
    stream.write((const char *) &beta, sizeof(float));
    stream.write((const char *) &t, sizeof(float));
    stream.write((const char *) &rand_gen, sizeof(PrngClass));
    stream.write((const char *) &action, sizeof(double));
    stream.write((const char *) &action_measured, sizeof(bool));

    for (int i = 0; i < arr_size_1; i++) {
        for (int j = 0; j < arr_size_2; j++) {
            for (int k = 0; k < arr_size_3; k++) {
                for (int l = 0; l < arr_size_4; l++) {
                    stream.write((const char *) &m[i][j][k][l],
                            sizeof(LinkNodeDim4<MatrixSU<Float, PrngClass>>));
                }
            }
        }
    }
}

template <template <typename _Float, class _PrngClass> class MatrixSU,
            template <typename Node> class DynamicArray,
            typename Float, class PrngClass>
void GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, DynamicArray, Float, PrngClass>::
                                            Base_Input (ifstream &stream) {
    stream.read((char *) &arr_size_1, sizeof(int));
    stream.read((char *) &arr_size_2, sizeof(int));
    stream.read((char *) &arr_size_3, sizeof(int));
    stream.read((char *) &arr_size_4, sizeof(int));

    stream.read((char *) &g0, sizeof(float));
    stream.read((char *) &beta, sizeof(float));
    stream.read((char *) &t, sizeof(float));
    stream.read((char *) &rand_gen, sizeof(PrngClass));
    stream.read((char *) &action, sizeof(double));
    stream.read((char *) &action_measured, sizeof(bool));



    for (int i = 0; i < arr_size_1; i++) {
        for (int j = 0; j < arr_size_2; j++) {
            for (int k = 0; k < arr_size_3; k++) {
                for (int l = 0; l < arr_size_4; l++) {
                    stream.read((char *) &m[i][j][k][l],
                            sizeof(LinkNodeDim4<MatrixSU<Float, PrngClass>>));
                }
            }
        }
    }
}


















template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class DynamicArray,
                typename Float, class PrngClass>
GluodynamicsWithoutBorderDim4_SU_Base<MatrixSU, DynamicArray, Float, PrngClass>
            ::GluodynamicsWithoutBorderDim4_SU_Base (int N1, int N2, int N3, int N4,
                                                     int mode, float g0_input, unsigned int _seed)
:   m(N1, N2, N3, N4),
    arr_size_1(N1),
    arr_size_2(N2),
    arr_size_3(N3),
    arr_size_4(N4)
{
    rand_gen.seed(_seed);
    g0 = g0_input;
    action_measured = false;
    if (g0 == 0) {
        throw invalid_argument("PeriodicGluodynamicsDim4_SU_Base_InitError_ZeroG0Constant");
        return;
    }
    beta = MatrixSU<Float, PrngClass>().GetRepresentationDimension()*2/g0/g0;
    t = 0;


    for (int i = 0; i < arr_size_1; i++) {
        for (int j = 0; j < arr_size_2; j++) {
            for (int k = 0; k < arr_size_3; k++) {
                for (int l = 0; l < arr_size_4; l++) {
                    for (int d = 1; d < 5; d++) {
                        m[i][j][k][l].up(d) = MatrixSU<Float, PrngClass>(mode,
                                                                    rand_gen());
                    }
                }
            }
        }
    }
}


template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
double GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                                                                ::Action() {
    if (action_measured) {
        return action;
    }

    double sum = 0.0;

    for (int i = 0; i < arr_size_1; i++) {
        for (int j = 0; j < arr_size_2; j++) {
            for (int k = 0; k < arr_size_3; k++) {
                for (int l = 0; l < arr_size_4; l++) {
                    sum += (GetLink(i, j, k, l, 1)*
                            GetLink(i+1, j, k, l, 2)*
                            GetLink(i, j+1, k, l, 1).Inversed()*
                            GetLink(i, j, k, l, 2).Inversed()).ReTrace()
                        + (GetLink(i, j, k, l, 1)*
                            GetLink(i+1, j, k, l, 3)*
                            GetLink(i, j, k+1, l, 1).Inversed()*
                            GetLink(i, j, k, l, 3).Inversed()).ReTrace()
                        + (GetLink(i, j, k, l, 1)*
                            GetLink(i+1, j, k, l, 4)*
                            GetLink(i, j, k, l+1, 1).Inversed()*
                            GetLink(i, j, k, l, 4).Inversed()).ReTrace()
                        + (GetLink(i, j, k, l, 2)*
                            GetLink(i, j+1, k, l, 3)*
                            GetLink(i, j, k+1, l, 2).Inversed()*
                            GetLink(i, j, k, l, 3).Inversed()).ReTrace()
                        + (GetLink(i, j, k, l, 2)*
                            GetLink(i, j+1, k, l, 4)*
                            GetLink(i, j, k, l+1, 2).Inversed()*
                            GetLink(i, j, k, l, 4).Inversed()).ReTrace()
                        + (GetLink(i, j, k, l, 3)*
                            GetLink(i, j, k+1, l, 4)*
                            GetLink(i, j, k, l+1, 3).Inversed()*
                            GetLink(i, j, k, l, 4).Inversed()).ReTrace();
                }
            }
        }
    }

    sum *= -beta/MatrixSU<Float, PrngClass>().GetRepresentationDimension();
    action = sum;
    action_measured = true;
    return sum;
}








template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
MatrixSU<Float, PrngClass> * GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
            ::PrecalculatedPlaquettes(int i, int j, int k, int l, int direction) const {
    MatrixSU<Float, PrngClass> *Usim;
    Usim = new MatrixSU<Float, PrngClass>[6];
    for (int i = 0; i < 6; i++) {
        Usim[i] = MatrixSU<Float, PrngClass>(0);
    }

    switch (direction) {
        case 1:
            Usim[0] *= GetLink(i, j, k, l, 1);
            Usim[0] *= GetLink(i+1, j, k, l, 2);
            Usim[0] *= GetLink(i, j+1, k, l, 1).Inversed();
            Usim[0] *= GetLink(i, j, k, l, 2).Inversed();

            Usim[1] *= GetLink(i, j, k, l, 1);
            Usim[1] *= GetLink(i+1, j-1, k, l, 2).Inversed();
            Usim[1] *= GetLink(i, j-1, k, l, 1).Inversed();
            Usim[1] *= GetLink(i, j-1, k, l, 2);

            Usim[2] *= GetLink(i, j, k, l, 1);
            Usim[2] *= GetLink(i+1, j, k, l, 3);
            Usim[2] *= GetLink(i, j, k+1, l, 1).Inversed();
            Usim[2] *= GetLink(i, j, k, l, 3).Inversed();

            Usim[3] *= GetLink(i, j, k, l, 1);
            Usim[3] *= GetLink(i+1, j, k-1, l, 3).Inversed();
            Usim[3] *= GetLink(i, j, k-1, l, 1).Inversed();
            Usim[3] *= GetLink(i, j, k-1, l, 3);

            Usim[4] *= GetLink(i, j, k, l, 1);
            Usim[4] *= GetLink(i+1, j, k, l, 4);
            Usim[4] *= GetLink(i, j, k, l+1, 1).Inversed();
            Usim[4] *= GetLink(i, j, k, l, 4).Inversed();

            Usim[5] *= GetLink(i, j, k, l, 1);
            Usim[5] *= GetLink(i+1, j, k, l-1, 4).Inversed();
            Usim[5] *= GetLink(i, j, k, l-1, 1).Inversed();
            Usim[5] *= GetLink(i, j, k, l-1, 4);

            break;
        case 2:
            Usim[0] *= GetLink(i, j, k, l, 2);
            Usim[0] *= GetLink(i, j+1, k, l, 1);
            Usim[0] *= GetLink(i+1, j, k, l, 2).Inversed();
            Usim[0] *= GetLink(i, j, k, l, 1).Inversed();

            Usim[1] *= GetLink(i, j, k, l, 2);
            Usim[1] *= GetLink(i-1, j+1, k, l, 1).Inversed();
            Usim[1] *= GetLink(i-1, j, k, l, 2).Inversed();
            Usim[1] *= GetLink(i-1, j, k, l, 1);

            Usim[2] *= GetLink(i, j, k, l, 2);
            Usim[2] *= GetLink(i, j+1, k, l, 3);
            Usim[2] *= GetLink(i, j, k+1, l, 2).Inversed();
            Usim[2] *= GetLink(i, j, k, l, 3).Inversed();

            Usim[3] *= GetLink(i, j, k, l, 2);
            Usim[3] *= GetLink(i, j+1, k-1, l, 3).Inversed();
            Usim[3] *= GetLink(i, j, k-1, l, 2).Inversed();
            Usim[3] *= GetLink(i, j, k-1, l, 3);

            Usim[4] *= GetLink(i, j, k, l, 2);
            Usim[4] *= GetLink(i, j+1, k, l, 4);
            Usim[4] *= GetLink(i, j, k, l+1, 2).Inversed();
            Usim[4] *= GetLink(i, j, k, l, 4).Inversed();

            Usim[5] *= GetLink(i, j, k, l, 2);
            Usim[5] *= GetLink(i, j+1, k, l-1, 4).Inversed();
            Usim[5] *= GetLink(i, j, k, l-1, 2).Inversed();
            Usim[5] *= GetLink(i, j, k, l-1, 4);

            break;
        case 3:
            Usim[0] *= GetLink(i, j, k, l, 3);
            Usim[0] *= GetLink(i, j, k+1, l, 1);
            Usim[0] *= GetLink(i+1, j, k, l, 3).Inversed();
            Usim[0] *= GetLink(i, j, k, l, 1).Inversed();

            Usim[1] *= GetLink(i, j, k, l, 3);
            Usim[1] *= GetLink(i-1, j, k+1, l, 1).Inversed();
            Usim[1] *= GetLink(i-1, j, k, l, 3).Inversed();
            Usim[1] *= GetLink(i-1, j, k, l, 1);

            Usim[2] *= GetLink(i, j, k, l, 3);
            Usim[2] *= GetLink(i, j, k+1, l, 2);
            Usim[2] *= GetLink(i, j+1, k, l, 3).Inversed();
            Usim[2] *= GetLink(i, j, k, l, 2).Inversed();

            Usim[3] *= GetLink(i, j, k, l, 3);
            Usim[3] *= GetLink(i, j-1, k+1, l, 2).Inversed();
            Usim[3] *= GetLink(i, j-1, k, l, 3).Inversed();
            Usim[3] *= GetLink(i, j-1, k, l, 2);

            Usim[4] *= GetLink(i, j, k, l, 3);
            Usim[4] *= GetLink(i, j, k+1, l, 4);
            Usim[4] *= GetLink(i, j, k, l+1, 3).Inversed();
            Usim[4] *= GetLink(i, j, k, l, 4).Inversed();

            Usim[5] *= GetLink(i, j, k, l, 3);
            Usim[5] *= GetLink(i, j, k+1, l-1, 4).Inversed();
            Usim[5] *= GetLink(i, j, k, l-1, 3).Inversed();
            Usim[5] *= GetLink(i, j, k, l-1, 4);

            break;
        case 4:
            Usim[0] *= GetLink(i, j, k, l, 4);
            Usim[0] *= GetLink(i, j, k, l+1, 1);
            Usim[0] *= GetLink(i+1, j, k, l, 4).Inversed();
            Usim[0] *= GetLink(i, j, k, l, 1).Inversed();

            Usim[1] *= GetLink(i, j, k, l, 4);
            Usim[1] *= GetLink(i-1, j, k, l+1, 1).Inversed();
            Usim[1] *= GetLink(i-1, j, k, l, 4).Inversed();
            Usim[1] *= GetLink(i-1, j, k, l, 1);

            Usim[2] *= GetLink(i, j, k, l, 4);
            Usim[2] *= GetLink(i, j, k, l+1, 2);
            Usim[2] *= GetLink(i, j+1, k, l, 4).Inversed();
            Usim[2] *= GetLink(i, j, k, l, 2).Inversed();

            Usim[3] *= GetLink(i, j, k, l, 4);
            Usim[3] *= GetLink(i, j-1, k, l+1, 2).Inversed();
            Usim[3] *= GetLink(i, j-1, k, l, 4).Inversed();
            Usim[3] *= GetLink(i, j-1, k, l, 2);

            Usim[4] *= GetLink(i, j, k, l, 4);
            Usim[4] *= GetLink(i, j, k, l+1, 3);
            Usim[4] *= GetLink(i, j, k+1, l, 4).Inversed();
            Usim[4] *= GetLink(i, j, k, l, 3).Inversed();

            Usim[5] *= GetLink(i, j, k, l, 4);
            Usim[5] *= GetLink(i, j, k-1, l+1, 3).Inversed();
            Usim[5] *= GetLink(i, j, k-1, l, 4).Inversed();
            Usim[5] *= GetLink(i, j, k-1, l, 3);

            break;
    }

    return Usim;
}








template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
bool GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                        ::MetropolisHit(int i, int j, int k, int l, int direction,
                                            MatrixSU<Float, PrngClass> *Usim) {
    MatrixSU<Float, PrngClass> rand_matrix;
    Float epsilon = (1.0/beta < 0.5) ? (1.0/beta) : 0.5;
//        Float epsilon = 0.01;
    rand_matrix.AroundIdentityMatrix(rand_gen(), epsilon);
    Float action_old = 0.0;
    Float action_new = 0.0;
    for (int x = 0; x < 6; x++) {
        action_old += Usim[x].ReTrace();
        action_new += (rand_matrix * Usim[x]).ReTrace();
    }
    action_old *= -beta/MatrixSU<Float, PrngClass>()
                                            .GetRepresentationDimension();
    action_new *= -beta/MatrixSU<Float, PrngClass>()
                                            .GetRepresentationDimension();

    float probability = (action_new < action_old) ? 1.0
                                            : exp(action_old - action_new);
    float rand_float = ((float)rand_gen() - rand_gen.min())
                        /(rand_gen.max() - rand_gen.min());

    if (rand_float < probability) {
        m[i][j][k][l].up(direction) = rand_matrix * m[i][j][k][l].up(direction);
        for (int x = 0; x < 6; x++) {
            Usim[x] = rand_matrix * Usim[x];
        }
        return true;
    }

    return false;
}






template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
unsigned int GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                            ::SmallMonteCarloStep(unsigned int number_of_hits) {
    action_measured = false;
    int i = rand_gen()%arr_size_1;
    int j = rand_gen()%arr_size_2;
    int k = rand_gen()%arr_size_3;
    int l = rand_gen()%arr_size_4;
    int direction = rand_gen() % 4 + 1;
    MatrixSU<Float, PrngClass> *Usim = PrecalculatedPlaquettes(i, j, k, l, direction);


    unsigned int number_of_changes = 0;
    for (unsigned int hit = 0; hit < number_of_hits; hit++) {
        if (MetropolisHit(i, j, k, l, direction, Usim)) {
            number_of_changes++;
        }
    }

    delete [] Usim;

    return number_of_changes;
}




template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
float GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
        ::/*HeatBath_*/MonteCarloStep(float number_of_steps, unsigned int number_of_hits) {
    action_measured = false;
    unsigned int all_hits_counter = 0;
    for (int successful_hits_counter = 0; successful_hits_counter <
            4*number_of_steps*arr_size_1*arr_size_2*arr_size_3*arr_size_4; ) {
        bool flag = false;
        while (!flag) {
            unsigned int steps = SmallMonteCarloStep(number_of_hits);
            flag = (steps > 0);
            all_hits_counter += number_of_hits;
            successful_hits_counter += steps;
        }
    }
    t += number_of_steps;

    return (4.0*number_of_steps*arr_size_1*arr_size_2*arr_size_3*arr_size_4)
                                                            /all_hits_counter;
}





template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
unsigned int GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                ::HeatBath_SmallMonteCarloStep(unsigned int number_of_hits,
                int i, int j, int k, int l, int direction) {
    action_measured = false;
//    int i = rand_gen()%arr_size_1;
//    int j = rand_gen()%arr_size_2;
//    int k = rand_gen()%arr_size_3;
//    int l = rand_gen()%arr_size_4;
//    int direction = rand_gen()%4 + 1;
    MatrixSU<Float, PrngClass> *Usim = PrecalculatedPlaquettes(i, j, k, l, direction);




//    MatrixSU<Float, PrngClass> ZXC(0);
//    m[i][j][k][l].up(direction) = MatrixSU<Float, PrngClass>(0, rand_gen());
    m[i][j][k][l].up(direction).HeatBathMatrix(rand_gen(), beta, Usim);

    delete [] Usim;

    cout << '$';
    return 1;
}




template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
float GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>::HeatBath_MonteCarloStep
                        (float number_of_steps, unsigned int number_of_hits) {
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            for (int k = 0; k < N3; k++) {
                for (int l = 0; l < N4; l++) {
                    for (int direction = 1; direction < 5; direction++) {
                        HeatBath_SmallMonteCarloStep
                                                (number_of_hits, i, j, k, l, direction);
                    }
                }
            }
        }
    }





    return 1.0;
}







template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
unsigned int GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                ::Timeslice_I_SmallMonteCarloStep(unsigned int number_of_hits,
                                                int mean_timeslice, int width) {
    action_measured = false;
    int i = rand_gen()%arr_size_1;
    int j = rand_gen()%arr_size_2;
    int k = rand_gen()%arr_size_3;
    int ran_num = rand_gen()%(4*width + 3);
    int l = mean_timeslice - width/2 + ran_num/4;
    l = (l >= 0 ? l%arr_size_4 : arr_size_4 + l%arr_size_4)%arr_size_4;
    int direction = ran_num % 4 + 1;
    MatrixSU<Float, PrngClass> *Usim = PrecalculatedPlaquettes(i, j, k, l, direction);




    unsigned int number_of_changes = 0;
    for (unsigned int hit = 0; hit < number_of_hits; hit++) {
        if (MetropolisHit(i, j, k, l, direction, Usim)) {
            number_of_changes++;
        }
    }

    delete [] Usim;

    return number_of_changes;
}




template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
float GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
        ::Timeslice_I_MonteCarloStep(float number_of_steps,
                unsigned int number_of_hits, int mean_timeslice, int width) {
    action_measured = false;
    unsigned int all_hits_counter = 0;
    for (int successful_hits_counter = 0; successful_hits_counter <
            (3 + 4*width)*number_of_steps*arr_size_1*arr_size_2*arr_size_3; ) {
        bool flag = false;
        while (!flag) {
            unsigned int steps = Timeslice_I_SmallMonteCarloStep
                                        (number_of_hits, mean_timeslice, width);
            flag = steps > 0;
            all_hits_counter += number_of_hits;
            successful_hits_counter += steps;
        }
    }
    t += number_of_steps;

    return (float(3 + 4*width)*number_of_steps*arr_size_1*arr_size_2*arr_size_3)
                                                            /all_hits_counter;
}







template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
unsigned int GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                ::Timeslice_II_SmallMonteCarloStep(unsigned int number_of_hits,
                                                int mean_timeslice, int width) {
    action_measured = false;
    int i = rand_gen()%arr_size_1;
    int j = rand_gen()%arr_size_2;
    int k = rand_gen()%arr_size_3;
    int ran_num = rand_gen()%(4*width - 3) + 3;
    int l = mean_timeslice - width/2 + ran_num/4;
    l = (l >= 0 ? l%arr_size_4 : arr_size_4 + l%arr_size_4)%arr_size_4;
    int direction = ran_num % 4 + 1;
    MatrixSU<Float, PrngClass> *Usim = PrecalculatedPlaquettes(i, j, k, l, direction);




    unsigned int number_of_changes = 0;
    for (unsigned int hit = 0; hit < number_of_hits; hit++) {
        if (MetropolisHit(i, j, k, l, direction, Usim)) {
            number_of_changes++;
        }
    }

    delete [] Usim;

    return number_of_changes;
}




template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
float GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
        ::Timeslice_II_MonteCarloStep(float number_of_steps,
                unsigned int number_of_hits, int mean_timeslice, int width) {
    action_measured = false;
    unsigned int all_hits_counter = 0;
    for (int successful_hits_counter = 0; successful_hits_counter <
            (4*width - 3)*number_of_steps*arr_size_1*arr_size_2*arr_size_3; ) {
        bool flag = false;
        while (!flag) {
            unsigned int steps = Timeslice_II_SmallMonteCarloStep
                                        (number_of_hits, mean_timeslice, width);
            flag = steps > 0;
            all_hits_counter += number_of_hits;
            successful_hits_counter += steps;
        }
    }
    t += number_of_steps;

    return (float(4*width - 3)*number_of_steps*arr_size_1*arr_size_2*arr_size_3)
                                                            /all_hits_counter;
}









//Wilson loop for rectangles IxJ
template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
double GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                                            ::AverageWilsonLoop(int I, int J) const {
    double sum = 0.0;

    for (int i_direction = 1; i_direction < 5; i_direction++) {
        for (int j_direction = 1; j_direction < 5; j_direction++) {
            if (i_direction == j_direction) {
                continue;
            }
            sum += AverageOrientedWilsonLoop(I, J, i_direction, j_direction);
        }
    }


    sum /= 12;
    return sum;
}








template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
double GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
            ::AverageOrientedWilsonLoop(int I, int J,
                                        int i_direction, int j_direction) const {
    if (I <= 0 || J <= 0 || i_direction < 1
                || j_direction < 1 || i_direction > 4 || j_direction > 4
                || i_direction == j_direction) {
        throw invalid_argument("PeriodicAverageOrientedWilsonLoop_WRONGPARAMETERS");
    }


    double sum = 0.0;
    for (int i = 0; i < arr_size_1; i++) {
        for (int j = 0; j < arr_size_2; j++) {
            for (int k = 0; k < arr_size_3; k++) {
                for (int l = 0; l < arr_size_4; l++) {
                    sum += SingleWilsonLoop
                                (I, J, i_direction, j_direction, i, j, k, l);
                }
            }
        }
    }

    sum /= arr_size_1*arr_size_2*arr_size_3*arr_size_4;
    return sum;
}








template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
double GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
            ::SingleWilsonLoop(int I, int J, int i_direction, int j_direction,
                    int i, int j, int k, int l) const {
    if (I <= 0 || J <= 0 || i_direction < 1
            || j_direction < 1 || i_direction > 4 || j_direction > 4
            || i_direction == j_direction) {
        throw invalid_argument("PeriodicSingleWilsonLoop_WRONGPARAMETERS");
    }

    double loop = 0.0;
    MatrixSU<Float, PrngClass> U(0);

    const int Shift = 1000;
    switch (i_direction*Shift + j_direction) {
        case (1*Shift + 2):
            for (int x = 0; x < I; x++) {
                U *= GetLink(i+x, j, k, l, 1);
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i+I, j+y, k, l, 2);
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i+I-1-x, j+J, k, l, 1).Inversed();
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j+J-1-y, k, l, 2).Inversed();
            }
            return U.ReTrace();


        case (2*Shift + 1):
            for (int y = 0; y < J; y++) {
                U *= GetLink(i+y, j, k, l, 1);
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i+J, j+x, k, l, 2);
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i+J-1-y, j+I, k, l, 1).Inversed();
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j+I-1-x, k, l, 2).Inversed();
            }
            return U.ReTrace();




        case (1*Shift + 3):
            for (int x = 0; x < I; x++) {
                U *= GetLink(i+x, j, k, l, 1);
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i+I, j, k+y, l, 3);
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i+I-1-x, j, k+J, l, 1).Inversed();
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j, k+J-1-y, l, 3).Inversed();
            }
            return U.ReTrace();


        case (3*Shift + 1):
            for (int y = 0; y < J; y++) {
                U *= GetLink(i+y, j, k, l, 1);
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i+J, j, k+x, l, 3);
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i+J-1-y, j, k+I, l, 1).Inversed();
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j, k+I-1-x, l, 3).Inversed();
            }
            return U.ReTrace();





        case (1*Shift + 4):
            for (int x = 0; x < I; x++) {
                U *= GetLink(i+x, j, k, l, 1);
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i+I, j, k, l+y, 4);
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i+I-1-x, j, k, l+J, 1).Inversed();
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j, k, l+J-1-y, 4).Inversed();
            }
            return U.ReTrace();


        case (4*Shift + 1):
            for (int y = 0; y < J; y++) {
                U *= GetLink(i+y, j, k, l, 1);
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i+J, j, k, l+x, 4);
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i+J-1-y, j, k, l+I, 1).Inversed();
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j, k, l+I-1-x, 4).Inversed();
            }
            return U.ReTrace();





        case (2*Shift + 3):
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j+x, k, l, 2);
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j+I, k+y, l, 3);
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j+I-1-x, k+J, l, 2).Inversed();
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j, k+J-1-y, l, 3).Inversed();
            }
            return U.ReTrace();


        case (3*Shift + 2):
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j+y, k, l, 2);
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j+J, k+x, l, 3);
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j+J-1-y, k+I, l, 2).Inversed();
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j, k+I-1-x, l, 3).Inversed();
            }
            return U.ReTrace();





        case (2*Shift + 4):
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j+x, k, l, 2);
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j+I, k, l+y, 4);
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j+I-1-x, k, l+J, 2).Inversed();
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j, k, l+J-1-y, 4).Inversed();
            }
            return U.ReTrace();
            break;

        case (4*Shift + 2):
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j+y, k, l, 2);
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j+J, k, l+x, 4);
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j+J-1-y, k, l+I, 2).Inversed();
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j, k, l+I-1-x, 4).Inversed();
            }
            return U.ReTrace();





        case (3*Shift + 4):
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j, k+x, l, 3);
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j, k+I, l+y, 4);
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j, k+I-1-x, l+J, 3).Inversed();
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j, k, l+J-1-y, 4).Inversed();
            }
            return U.ReTrace();


        case (4*Shift + 3):
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j, k+y, l, 3);
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j, k+J, l+x, 4);
            }
            for (int y = 0; y < J; y++) {
                U *= GetLink(i, j, k+J-1-y, l+I, 3).Inversed();
            }
            for (int x = 0; x < I; x++) {
                U *= GetLink(i, j, k, l+I-1-x, 4).Inversed();
            }
            return U.ReTrace();
    }

    return loop;
}



template <template <typename _Float, class _PrngClass> class MatrixSU,
        template <typename Node> class CycledArrayDim4,
        typename Float, class PrngClass>
double GluodynamicsWithoutBorderDim4_SU_Base
    <MatrixSU, CycledArrayDim4, Float, PrngClass>
    ::SingleKnot(int a_direction, int b_direction, int c_direction, int d_direction,
                  int i, int j, int k, int l) const {
    if (abs(a_direction) == abs(b_direction) ||
        abs(a_direction) == abs(c_direction) ||
        abs(a_direction) == abs(d_direction) ||
        abs(b_direction) == abs(c_direction) ||
        abs(b_direction) == abs(d_direction) ||
        abs(c_direction) == abs(d_direction) ||
        abs(a_direction) < 1 || abs(a_direction) > 4 ||
        abs(b_direction) < 1 || abs(b_direction) > 4 ||
        abs(c_direction) < 1 || abs(c_direction) > 4 ||
        abs(d_direction) < 1 || abs(d_direction) > 4) {
        throw invalid_argument("PeriodicSingleKnot_WRONGPARAMETERS");
    }

    int A[4], B[4], C[4], D[4];
    for (int step = 0; step < 4; step++) {
        A[step] = B[step] = C[step] = D[step] = 0;
    }
    A[abs(a_direction) - 1] += (a_direction > 0) ? +1 : -1;
    B[abs(b_direction) - 1] += (b_direction > 0) ? +1 : -1;
    C[abs(c_direction) - 1] += (c_direction > 0) ? +1 : -1;
    D[abs(d_direction) - 1] += (d_direction > 0) ? +1 : -1;

    MatrixSU<Float, PrngClass> U(0);
    U *= GetLink(i, j, k, l, a_direction) *
         GetLink(i + A[0], j + A[1], k + A[2], l + A[3], b_direction) *
         GetLink(i + B[0], j + B[1], k + B[2], l + B[3], a_direction).Inversed() *
         GetLink(i, j, k, l, b_direction).Inversed() *
         GetLink(i, j, k, l, c_direction) *
         GetLink(i + C[0], j + C[1], k + C[2], l + C[3], d_direction) *
         GetLink(i + D[0], j + D[1], k + D[2], l + D[3], c_direction).Inversed() *
         GetLink(i, j, k, l, d_direction).Inversed();

    return U.ReTrace();
}





template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
double GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, CycledArrayDim4, Float, PrngClass>
                                                            ::AverageFace() const {
    return 1 - AverageWilsonLoop(1, 1)/MatrixSU<Float, PrngClass>()
                                                .GetRepresentationDimension();
}













template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
class PeriodicGluodynamicsDim4_SU;

template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
ofstream &operator<< (ofstream &stream,
                                    const PeriodicGluodynamicsDim4_SU<MatrixSU,
                                    CycledArrayDim4, Float, PrngClass> &System);

template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
ifstream &operator>> (ifstream &stream, PeriodicGluodynamicsDim4_SU<MatrixSU,
                            CycledArrayDim4, Float, PrngClass> &System);







template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class DynamicArray,
                typename Float, class PrngClass>
class PeriodicGluodynamicsDim4_SU : public GluodynamicsWithoutBorderDim4_SU_Base
        <MatrixSU, DynamicArray, Float, PrngClass>
{
    public:
        PeriodicGluodynamicsDim4_SU() {}
        PeriodicGluodynamicsDim4_SU(int _N1, int _N2, int _N3, int _N4,
                                int mode, float g0_input, unsigned int seed);

        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, DynamicArray, Float, PrngClass>::m;

        MatrixSU<Float, PrngClass> GetLink(int x, int y, int z, int t,
                                                            int direction) const override;


        friend ofstream &operator<<
            <MatrixSU, DynamicArray, Float, PrngClass>
            (ofstream &stream, const PeriodicGluodynamicsDim4_SU<MatrixSU,
                            DynamicArray, Float, PrngClass> &System);
        friend ifstream &operator>>
            <MatrixSU, DynamicArray, Float, PrngClass>
            (ifstream &stream, PeriodicGluodynamicsDim4_SU<MatrixSU,
                            DynamicArray, Float, PrngClass> &System);



        double AverageOrientedPolyakovLoop(int direction) const;
        double SinglePolyakovLoop(int direction, int i, int j, int k) const;


        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, DynamicArray, Float, PrngClass>::rand_gen;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, DynamicArray, Float, PrngClass>::arr_size_1;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, DynamicArray, Float, PrngClass>::arr_size_2;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, DynamicArray, Float, PrngClass>::arr_size_3;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, DynamicArray, Float, PrngClass>::arr_size_4;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, DynamicArray, Float, PrngClass>::g0;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, DynamicArray, Float, PrngClass>::beta;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, DynamicArray, Float, PrngClass>::t;
};




template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
ofstream &operator<< (ofstream &stream, const PeriodicGluodynamicsDim4_SU
                    <MatrixSU, CycledArrayDim4, Float, PrngClass> &System) {
    System.Base_Output(stream);

    return stream;
}

template <template <typename _Float, class _PrngClass> class MatrixSU,
            template <typename Node> class CycledArrayDim4,
            typename Float, class PrngClass>
ifstream &operator>> (ifstream &stream, PeriodicGluodynamicsDim4_SU<MatrixSU,
                            CycledArrayDim4, Float, PrngClass> &System) {
    System.Base_Input(stream);

    return stream;
}




template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class SafeArrayDim4,
                typename Float, class PrngClass>
MatrixSU<Float, PrngClass> PeriodicGluodynamicsDim4_SU<MatrixSU,
                                    SafeArrayDim4, Float, PrngClass>
                ::GetLink   (int x, int y, int z, int t, int direction) const {
    if (direction < 0) {
        switch (direction) {
            case -1:
                x--;
                break;
            case -2:
                y--;
                break;
            case -3:
                z--;
                break;
            case -4:
                t--;
                break;
        }
    } else if (direction == 0) {
        throw runtime_error("GetLinkPeriodic");
    }

    x = (x >= 0 ? x%arr_size_1 : arr_size_1 + x%arr_size_1)%arr_size_1;
    y = (y >= 0 ? y%arr_size_2 : arr_size_2 + y%arr_size_2)%arr_size_2;
    z = (z >= 0 ? z%arr_size_3 : arr_size_3 + z%arr_size_3)%arr_size_3;
    t = (t >= 0 ? t%arr_size_4 : arr_size_4 + t%arr_size_4)%arr_size_4;

    if (direction > 0) {
        return m[x][y][z][t].copy_of_up(direction);
    } else {
        return m[x][y][z][t].copy_of_up(-direction).Inversed();
    }
}





















template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
PeriodicGluodynamicsDim4_SU<MatrixSU, CycledArrayDim4, Float, PrngClass>
            ::PeriodicGluodynamicsDim4_SU (int N1, int N2, int N3, int N4,
                            int mode, float g0_input, unsigned int _seed)
:GluodynamicsWithoutBorderDim4_SU_Base<MatrixSU, CycledArrayDim4, Float, PrngClass>
                                        (N1, N2, N3, N4, mode, g0_input, _seed)
{
}









template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
double PeriodicGluodynamicsDim4_SU<MatrixSU, CycledArrayDim4, Float, PrngClass>
                        ::AverageOrientedPolyakovLoop(int direction) const {
    if (direction < 1 || direction > 4) {
        throw invalid_argument("PeriodicAverageOrientedPolyakovLoop_WRONGPARAMETERS");
    }


    double sum = 0;
    switch (direction) {
        case 1:
            for (int j = 0; j < arr_size_2; j++) {
                for (int k = 0; k < arr_size_3; k++) {
                    for (int l = 0; l < arr_size_4; l++) {
                        MatrixSU<Float, PrngClass> U(0);

                        for (int x = 0; x < arr_size_1; x++) {
                            U *= GetLink(x, j, k, l, direction);
                        }
                        sum += U.ReTrace();
                    }
                }
            }
            sum /= arr_size_2*arr_size_3*arr_size_4;
            break;


        case 2:
            for (int i = 0; i < arr_size_1; i++) {
                for (int k = 0; k < arr_size_3; k++) {
                    for (int l = 0; l < arr_size_4; l++) {
                        MatrixSU<Float, PrngClass> U(0);

                        for (int x = 0; x < arr_size_2; x++) {
                            U *= GetLink(i, x, k, l, direction);
                        }
                        sum += U.ReTrace();
                    }
                }
            }
            sum /= arr_size_1*arr_size_3*arr_size_4;
            break;

        case 3:
            for (int i = 0; i < arr_size_1; i++) {
                for (int j = 0; j < arr_size_2; j++) {
                    for (int l = 0; l < arr_size_4; l++) {
                        MatrixSU<Float, PrngClass> U(0);

                        for (int x = 0; x < arr_size_3; x++) {
                            U *= GetLink(i, j, x, l, direction);
                        }
                        sum += U.ReTrace();
                    }
                }
            }
            sum /= arr_size_1*arr_size_2*arr_size_4;
            break;



        case 4:
            for (int i = 0; i < arr_size_1; i++) {
                for (int j = 0; j < arr_size_2; j++) {
                    for (int k = 0; k < arr_size_3; k++) {
                        MatrixSU<Float, PrngClass> U(0);

                        for (int x = 0; x < arr_size_4; x++) {
                            U *= GetLink(i, j, k, x, direction);
                        }
                        sum += U.ReTrace();

                    }
                }
            }
            sum /= arr_size_1*arr_size_2*arr_size_3;
            break;
    }

    return sum;
}







template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class CycledArrayDim4,
                typename Float, class PrngClass>
double PeriodicGluodynamicsDim4_SU<MatrixSU, CycledArrayDim4, Float, PrngClass>
                    ::SinglePolyakovLoop(int direction, int i, int j, int k) const {
    if (direction < 1 || direction > 4) {
        throw invalid_argument("PeriodicAverageOrientedPolyakovLoop_WRONGPARAMETERS");
    }


    MatrixSU<Float, PrngClass> U(0);
    switch (direction) {
        case 1:
            for (int x = 0; x < arr_size_1; x++) {
                U *= GetLink(x, i, j, k, direction);
            }
            return U.ReTrace();

        case 2:
            for (int x = 0; x < arr_size_2; x++) {
                U *= GetLink(i, x, j, k, direction);
            }
            return U.ReTrace();

        case 3:
            for (int x = 0; x < arr_size_3; x++) {
                U *= GetLink(i, j, x, k, direction);
            }
            return U.ReTrace();

        case 4:
            for (int x = 0; x < arr_size_4; x++) {
                U *= GetLink(i, j, k, x, direction);
            }
            return U.ReTrace();
    }

    throw runtime_error("Unknown direction of Polyakov loop");
    return 0.0;
}























template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class SafeArrayDim4,
                typename Float, class PrngClass>
class xReflectionGluodynamicsDim4_SU;

template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class SafeArrayDim4,
                typename Float, class PrngClass>
ofstream &operator<< (ofstream &stream,
                            const xReflectionGluodynamicsDim4_SU<MatrixSU,
                                SafeArrayDim4, Float, PrngClass> &System);

template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class SafeArrayDim4,
                typename Float, class PrngClass>
ifstream &operator>> (ifstream &stream, xReflectionGluodynamicsDim4_SU<MatrixSU,
                            SafeArrayDim4, Float, PrngClass> &System);








template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class SafeArrayDim4,
                typename Float, class PrngClass>
class xReflectionGluodynamicsDim4_SU : public GluodynamicsWithoutBorderDim4_SU_Base
                                    <MatrixSU, SafeArrayDim4, Float, PrngClass>
{
    public:
        xReflectionGluodynamicsDim4_SU() {}
        xReflectionGluodynamicsDim4_SU(int _N1, int _N2, int _N3, int _N4,
                                int mode, float g0_input, unsigned int seed);

        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, SafeArrayDim4, Float, PrngClass>::m;

        MatrixSU<Float, PrngClass> GetLink(int x, int y, int z, int t,
                                                            int direction) const override;

        friend ofstream &operator<<
            <MatrixSU, SafeArrayDim4, Float, PrngClass>
            (ofstream &stream, const xReflectionGluodynamicsDim4_SU<MatrixSU,
                        SafeArrayDim4, Float, PrngClass> &System);
        friend ifstream &operator>>
            <MatrixSU, SafeArrayDim4, Float, PrngClass>
            (ifstream &stream, xReflectionGluodynamicsDim4_SU<MatrixSU,
                                SafeArrayDim4, Float, PrngClass> &System);






        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, SafeArrayDim4, Float, PrngClass>::rand_gen;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, SafeArrayDim4, Float, PrngClass>::arr_size_1;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, SafeArrayDim4, Float, PrngClass>::arr_size_2;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, SafeArrayDim4, Float, PrngClass>::arr_size_3;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, SafeArrayDim4, Float, PrngClass>::arr_size_4;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, SafeArrayDim4, Float, PrngClass>::g0;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, SafeArrayDim4, Float, PrngClass>::beta;
        using GluodynamicsWithoutBorderDim4_SU_Base
            <MatrixSU, SafeArrayDim4, Float, PrngClass>::t;
};




template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class SafeArrayDim4,
                typename Float, class PrngClass>
ofstream &operator<< (ofstream &stream, const xReflectionGluodynamicsDim4_SU
                    <MatrixSU, SafeArrayDim4, Float, PrngClass> &System) {
    System.Base_Output(stream);

    return stream;
}

template <template <typename _Float, class _PrngClass> class MatrixSU,
            template <typename Node> class SafeArrayDim4,
            typename Float, class PrngClass>
ifstream &operator>> (ifstream &stream, xReflectionGluodynamicsDim4_SU<MatrixSU,
                            SafeArrayDim4, Float, PrngClass> &System) {
    System.Base_Input(stream);

    return stream;
}






template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class SafeArrayDim4,
                typename Float, class PrngClass>
MatrixSU<Float, PrngClass> xReflectionGluodynamicsDim4_SU<MatrixSU,
                                    SafeArrayDim4, Float, PrngClass>
                        ::GetLink (int x, int y, int z, int t, int direction) const {
    if (direction < 0) {
        switch (direction) {
            case -1:
                x--;
                break;
            case -2:
                y--;
                break;
            case -3:
                z--;
                break;
            case -4:
                t--;
                break;
        }
    } else if (direction == 0) {
        throw runtime_error("GetLinkMobius");
    }


    t = (t >= 0 ? t%(2*arr_size_4) : 2*arr_size_4 + t%(2*arr_size_4))%
                                                                (2*arr_size_4);

    if (t/arr_size_4 == 1) {
        if (direction == 1) {
            x = arr_size_1 - x - 1;
            x = (x >= 0 ? x%arr_size_1 : arr_size_1 + x%arr_size_1)%arr_size_1;
            y = (y >= 0 ? y%arr_size_2 : arr_size_2 + y%arr_size_2)%arr_size_2;
            z = (z >= 0 ? z%arr_size_3 : arr_size_3 + z%arr_size_3)%arr_size_3;
            t = (t >= 0 ? t%arr_size_4 : arr_size_4 + t%arr_size_4)%arr_size_4;

            if (direction > 0) {
                return m[x][y][z][t].copy_of_up(direction).Inversed();
            } else {
                return m[x][y][z][t].copy_of_up(-direction);
            }
        } else {
            x = arr_size_1 - x;
            x = (x >= 0 ? x%arr_size_1 : arr_size_1 + x%arr_size_1)%arr_size_1;
            y = (y >= 0 ? y%arr_size_2 : arr_size_2 + y%arr_size_2)%arr_size_2;
            z = (z >= 0 ? z%arr_size_3 : arr_size_3 + z%arr_size_3)%arr_size_3;
            t = (t >= 0 ? t%arr_size_4 : arr_size_4 + t%arr_size_4)%arr_size_4;

            if (direction > 0) {
                return m[x][y][z][t].copy_of_up(direction);
            } else {
                return m[x][y][z][t].copy_of_up(-direction).Inversed();
            }
        }
    } else {
        x = (x >= 0 ? x%arr_size_1 : arr_size_1 + x%arr_size_1)%arr_size_1;
        y = (y >= 0 ? y%arr_size_2 : arr_size_2 + y%arr_size_2)%arr_size_2;
        z = (z >= 0 ? z%arr_size_3 : arr_size_3 + z%arr_size_3)%arr_size_3;

        if (direction > 0) {
            return m[x][y][z][t].copy_of_up(direction);
        } else {
            return m[x][y][z][t].copy_of_up(-direction).Inversed();
        }
    }
}




template <template <typename _Float, class _PrngClass> class MatrixSU,
                template <typename Node> class SafeArrayDim4,
                typename Float, class PrngClass>
xReflectionGluodynamicsDim4_SU<MatrixSU, SafeArrayDim4, Float, PrngClass>
            ::xReflectionGluodynamicsDim4_SU (int N1, int N2, int N3, int N4,
                            int mode, float g0_input, unsigned int _seed)
:   GluodynamicsWithoutBorderDim4_SU_Base<MatrixSU, SafeArrayDim4, Float, PrngClass>
                                        (N1, N2, N3, N4, mode, g0_input, _seed)
{
}




#endif // GLUODYNAMICSDIM4_SU_H
