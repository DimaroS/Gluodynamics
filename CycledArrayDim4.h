#ifndef ARRAYDIM4_H
#define ARRAYDIM4_H

#include "main.h"
#include <new>

using namespace std;

const int MAX_ARR_SIZE = MAX_SIZE_DIM4;



template <typename Node>
class CycledArrayDim1
{
    private:
        int arr_size;
    public:
        CycledArrayDim1(int s);
        CycledArrayDim1();
        Node m[MAX_ARR_SIZE];

        Node& operator[] (int k);
};

template <typename Node>
class CycledArrayDim2
{
    private:
        int arr_size;
    public:
        CycledArrayDim2(int s);
        CycledArrayDim2();
        CycledArrayDim1<Node> m[MAX_ARR_SIZE];

        CycledArrayDim1<Node>& operator[] (int k);
};

template <typename Node>
class CycledArrayDim3
{
    private:
        int arr_size;
    public:
        CycledArrayDim3(int s);
        CycledArrayDim3();
        CycledArrayDim2<Node> m[MAX_ARR_SIZE];

        CycledArrayDim2<Node>& operator[] (int k);
};



template <typename Node>
class CycledArrayDim4
{
    private:
        int arr_size;
    public:
        CycledArrayDim4(int s);
        CycledArrayDim4();
        CycledArrayDim3<Node> m[MAX_ARR_SIZE];

        CycledArrayDim3<Node>& operator[] (int k);
};






class CycledArrayDim1_InitError_UNSUPPORTEDSIZE
{
    int k = 0;
};

class CycledArrayDim2_InitError_UNSUPPORTEDSIZE
{
    int k = 0;
};

class CycledArrayDim3_InitError_UNSUPPORTEDSIZE
{
    int k = 0;
};

class CycledArrayDim4_InitError_UNSUPPORTEDSIZE
{
    int k = 0;
};




template <typename Node>
CycledArrayDim1<Node>::CycledArrayDim1()
:arr_size(0)
{
    return;
}

template <typename Node>
CycledArrayDim1<Node>::CycledArrayDim1(int s)
:arr_size(s)
{
    if (s <= 0 || s > MAX_ARR_SIZE) {
        throw(CycledArrayDim1_InitError_UNSUPPORTEDSIZE());
    }
}


template <typename Node>
Node& CycledArrayDim1<Node>::operator[] (int k) {
    return m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}





template <typename Node>
CycledArrayDim2<Node>::CycledArrayDim2()
:arr_size(0)
{
    return;
}

template <typename Node>
CycledArrayDim2<Node>::CycledArrayDim2(int s)
:arr_size(s)
{
    if (s <= 0 || s > MAX_ARR_SIZE) {
        throw(CycledArrayDim2_InitError_UNSUPPORTEDSIZE());
    }

    for (int i = 0; i < MAX_ARR_SIZE; i++) {
        m[i] = CycledArrayDim1<Node>(s);
    }
}


template <typename Node>
CycledArrayDim1<Node> &CycledArrayDim2<Node>::operator[] (int k) {
    return m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}





template <typename Node>
CycledArrayDim3<Node>::CycledArrayDim3()
:arr_size(0)
{
    return;
}

template <typename Node>
CycledArrayDim3<Node>::CycledArrayDim3(int s)
:arr_size(s)
{
    if (s <= 0 || s > MAX_ARR_SIZE) {
        throw(CycledArrayDim3_InitError_UNSUPPORTEDSIZE());
    }

    for (int i = 0; i < MAX_ARR_SIZE; i++) {
        m[i] = CycledArrayDim2<Node>(s);
    }
}


template <typename Node>
CycledArrayDim2<Node>& CycledArrayDim3<Node>::operator[] (int k) {
    return m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}






template <typename Node>
CycledArrayDim4<Node>::CycledArrayDim4()
:arr_size(0)
{
    return;
}

template <typename Node>
CycledArrayDim4<Node>::CycledArrayDim4(int s)
:arr_size(s)
{
    if (s <= 0 || s > MAX_ARR_SIZE) {
        throw(CycledArrayDim4_InitError_UNSUPPORTEDSIZE());
    }

    for (int i = 0; i < MAX_ARR_SIZE; i++) {
        m[i] = CycledArrayDim3<Node>(s);
    }
}


template <typename Node>
CycledArrayDim3<Node>& CycledArrayDim4<Node>::operator[] (int k) {
    return m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}

























template <typename Node>
class DynamicCycledArrayDim1
{
    private:
        int arr_size;
    public:
        DynamicCycledArrayDim1(int N);
        DynamicCycledArrayDim1();
        DynamicCycledArrayDim1(const DynamicCycledArrayDim1 &System);
        ~DynamicCycledArrayDim1();
        Node **m;

        Node& operator[] (int k);
        const Node& operator[] (int k) const;
};

template <typename Node>
class DynamicCycledArrayDim2
{
    private:
        int arr_size;
    public:
        DynamicCycledArrayDim2(int N1, int N2);
        DynamicCycledArrayDim2();
        DynamicCycledArrayDim2(const DynamicCycledArrayDim2 &System);
        ~DynamicCycledArrayDim2();
        DynamicCycledArrayDim1<Node> **m;

        DynamicCycledArrayDim1<Node>& operator[] (int k);
        const DynamicCycledArrayDim1<Node>& operator[] (int k) const;
};

template <typename Node>
class DynamicCycledArrayDim3
{
    private:
        int arr_size;
    public:
        DynamicCycledArrayDim3(int N1, int N2, int N3);
        DynamicCycledArrayDim3();
        DynamicCycledArrayDim3(const DynamicCycledArrayDim3 &System);
        ~DynamicCycledArrayDim3();
        DynamicCycledArrayDim2<Node> **m;


        DynamicCycledArrayDim2<Node>& operator[] (int k);
        const DynamicCycledArrayDim2<Node>& operator[] (int k) const;
};



template <typename Node>
class DynamicCycledArrayDim4
{
    private:
        int arr_size;
    public:
        DynamicCycledArrayDim4(int N1, int N2, int N3, int N4);
        DynamicCycledArrayDim4();
        DynamicCycledArrayDim4(const DynamicCycledArrayDim4 &System);
        ~DynamicCycledArrayDim4();
        DynamicCycledArrayDim3<Node> **m;

        DynamicCycledArrayDim3<Node>& operator[] (int k);
        const DynamicCycledArrayDim3<Node>& operator[] (int k) const;
};






class DynamicCycledArrayDim1_InitError_UNSUPPORTEDSIZE
{
    int k = 0;
};

class DynamicCycledArrayDim2_InitError_UNSUPPORTEDSIZE
{
    int k = 0;
};

class DynamicCycledArrayDim3_InitError_UNSUPPORTEDSIZE
{
    int k = 0;
};

class DynamicCycledArrayDim4_InitError_UNSUPPORTEDSIZE
{
    int k = 0;
};








template <typename Node>
DynamicCycledArrayDim1<Node>::DynamicCycledArrayDim1()
:arr_size(0)
{
    m = nullptr;
}


template <typename Node>
DynamicCycledArrayDim1<Node>::DynamicCycledArrayDim1
                            (const DynamicCycledArrayDim1 &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new Node*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new Node(*System.m[i]);
    }
}

template <typename Node>
DynamicCycledArrayDim1<Node>::DynamicCycledArrayDim1(int N)
:arr_size(N)
{
    if (arr_size <= 0) {
        throw(DynamicCycledArrayDim1_InitError_UNSUPPORTEDSIZE());
    }

    m = new Node*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new Node;
    }
}

template <typename Node>
DynamicCycledArrayDim1<Node>::~DynamicCycledArrayDim1() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node>
Node& DynamicCycledArrayDim1<Node>::operator[] (int k) {
    return *m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}

template <typename Node>
const Node& DynamicCycledArrayDim1<Node>::operator[] (int k) const {
    return *m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}





template <typename Node>
DynamicCycledArrayDim2<Node>::DynamicCycledArrayDim2() {
    arr_size = 0;
    m = nullptr;
}

template <typename Node>
DynamicCycledArrayDim2<Node>::DynamicCycledArrayDim2
                            (const DynamicCycledArrayDim2 &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new DynamicCycledArrayDim1<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicCycledArrayDim1<Node>(*System.m[i]);
    }
}

template <typename Node>
DynamicCycledArrayDim2<Node>::DynamicCycledArrayDim2(int N1, int N2) {
    if (N1 <= 0 || N2 <= 0) {
        throw(DynamicCycledArrayDim2_InitError_UNSUPPORTEDSIZE());
    }

    arr_size = N1;
    m = new DynamicCycledArrayDim1<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicCycledArrayDim1<Node>(N2);
    }
}

template <typename Node>
DynamicCycledArrayDim2<Node>::~DynamicCycledArrayDim2() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node>
DynamicCycledArrayDim1<Node>& DynamicCycledArrayDim2<Node>::operator[] (int k) {
    return *m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}

template <typename Node>
const DynamicCycledArrayDim1<Node>& DynamicCycledArrayDim2<Node>::operator[] (int k) const {
    return *m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}





template <typename Node>
DynamicCycledArrayDim3<Node>::DynamicCycledArrayDim3() {
    arr_size = 0;
    m = nullptr;
}

template <typename Node>
DynamicCycledArrayDim3<Node>::DynamicCycledArrayDim3
                            (const DynamicCycledArrayDim3 &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new DynamicCycledArrayDim2<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicCycledArrayDim2<Node>(*System.m[i]);
    }
}

template <typename Node>
DynamicCycledArrayDim3<Node>::DynamicCycledArrayDim3(int N1, int N2, int N3) {
    if (N1 <= 0 || N2 <= 0 || N3 <= 0) {
        throw(DynamicCycledArrayDim3_InitError_UNSUPPORTEDSIZE());
    }

    arr_size = N1;
    m = new DynamicCycledArrayDim2<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicCycledArrayDim2<Node>(N2, N3);
    }
}

template <typename Node>
DynamicCycledArrayDim3<Node>::~DynamicCycledArrayDim3() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node>
DynamicCycledArrayDim2<Node>& DynamicCycledArrayDim3<Node>::operator[] (int k) {
    return *m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}

template <typename Node>
const DynamicCycledArrayDim2<Node>& DynamicCycledArrayDim3<Node>::operator[] (int k) const {
    return *m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}






template <typename Node>
DynamicCycledArrayDim4<Node>::DynamicCycledArrayDim4() {
    arr_size = 0;
    m = nullptr;
}

template <typename Node>
DynamicCycledArrayDim4<Node>::DynamicCycledArrayDim4
                            (const DynamicCycledArrayDim4 &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new DynamicCycledArrayDim3<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicCycledArrayDim3<Node>(*System.m[i]);
    }
}

template <typename Node>
DynamicCycledArrayDim4<Node>::DynamicCycledArrayDim4(int N1, int N2, int N3, int N4) {
    if (N1 <= 0 || N2 <= 0 || N3 <= 0 || N4 <= 0) {
        throw(DynamicCycledArrayDim4_InitError_UNSUPPORTEDSIZE());
    }

    arr_size = N1;
    m = new DynamicCycledArrayDim3<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicCycledArrayDim3<Node>(N2, N3, N4);
    }
}

template <typename Node>
DynamicCycledArrayDim4<Node>::~DynamicCycledArrayDim4() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node>
DynamicCycledArrayDim3<Node> &DynamicCycledArrayDim4<Node>::operator[] (int k) {
    return *m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}

template <typename Node>
const DynamicCycledArrayDim3<Node> &DynamicCycledArrayDim4<Node>::operator[] (int k) const {
    return *m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}













template <typename Node>
class DynamicSafeArrayDim1
{
    private:
        int arr_size;
    public:
        DynamicSafeArrayDim1(int N);
        DynamicSafeArrayDim1();
        DynamicSafeArrayDim1(const DynamicSafeArrayDim1 &System);
        ~DynamicSafeArrayDim1();
        Node **m;

        Node& operator[] (int k);
};

template <typename Node>
class DynamicSafeArrayDim2
{
    private:
        int arr_size;
    public:
        DynamicSafeArrayDim2(int N1, int N2);
        DynamicSafeArrayDim2();
        DynamicSafeArrayDim2(const DynamicSafeArrayDim2 &System);
        ~DynamicSafeArrayDim2();
        DynamicSafeArrayDim1<Node> **m;

        DynamicSafeArrayDim1<Node>& operator[] (int k);
};

template <typename Node>
class DynamicSafeArrayDim3
{
    private:
        int arr_size;
    public:
        DynamicSafeArrayDim3(int N1, int N2, int N3);
        DynamicSafeArrayDim3();
        DynamicSafeArrayDim3(const DynamicSafeArrayDim3 &System);
        ~DynamicSafeArrayDim3();
        DynamicSafeArrayDim2<Node> **m;

        DynamicSafeArrayDim2<Node>& operator[] (int k);
};



template <typename Node>
class DynamicSafeArrayDim4
{
    private:
        int arr_size;
    public:
        DynamicSafeArrayDim4(int N1, int N2, int N3, int N4);
        DynamicSafeArrayDim4();
        DynamicSafeArrayDim4(const DynamicSafeArrayDim4 &System);
        ~DynamicSafeArrayDim4();
        DynamicSafeArrayDim3<Node> **m;

        DynamicSafeArrayDim3<Node>& operator[] (int k);
};






class DynamicSafeArrayDim1_WrongIndex
{
    public:
        DynamicSafeArrayDim1_WrongIndex() {}
    private:
        int k = 0;
};

class DynamicSafeArrayDim2_WrongIndex
{
    public:
        DynamicSafeArrayDim2_WrongIndex() {}
    private:
        int k = 0;
};

class DynamicSafeArrayDim3_WrongIndex
{
    public:
        DynamicSafeArrayDim3_WrongIndex() {}
    private:
        int k = 0;
};

class DynamicSafeArrayDim4_WrongIndex
{
    public:
        DynamicSafeArrayDim4_WrongIndex() {}
    private:
        int k = 0;
};



class DynamicSafeArrayDim1_InitError_UNSUPPORTEDSIZE
{
    public:
        DynamicSafeArrayDim1_InitError_UNSUPPORTEDSIZE() {}
    private:
        int k = 0;
};

class DynamicSafeArrayDim2_InitError_UNSUPPORTEDSIZE
{
    public:
        DynamicSafeArrayDim2_InitError_UNSUPPORTEDSIZE() {}
    private:
        int k = 0;
};

class DynamicSafeArrayDim3_InitError_UNSUPPORTEDSIZE
{
    public:
        DynamicSafeArrayDim3_InitError_UNSUPPORTEDSIZE() {}
    private:
        int k = 0;
};

class DynamicSafeArrayDim4_InitError_UNSUPPORTEDSIZE
{
    public:
        DynamicSafeArrayDim4_InitError_UNSUPPORTEDSIZE() {}
    private:
        int k = 0;
};










template <typename Node>
DynamicSafeArrayDim1<Node>::DynamicSafeArrayDim1()
:arr_size(0)
{
    m = nullptr;
}


template <typename Node>
DynamicSafeArrayDim1<Node>::DynamicSafeArrayDim1
                            (const DynamicSafeArrayDim1 &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new Node*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new Node(*System.m[i]);
    }
}

template <typename Node>
DynamicSafeArrayDim1<Node>::DynamicSafeArrayDim1(int N)
:arr_size(N)
{
    if (arr_size <= 0) {
        throw(DynamicSafeArrayDim1_InitError_UNSUPPORTEDSIZE());
    }

    m = new Node*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new Node;
    }
}

template <typename Node>
DynamicSafeArrayDim1<Node>::~DynamicSafeArrayDim1() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node>
Node& DynamicSafeArrayDim1<Node>::operator[] (int k) {
    if (k < 0 || k >= arr_size) {
        throw(DynamicSafeArrayDim1_WrongIndex());
    }
    return *m[k];
}






template <typename Node>
DynamicSafeArrayDim2<Node>::DynamicSafeArrayDim2() {
    arr_size = 0;
    m = nullptr;
}

template <typename Node>
DynamicSafeArrayDim2<Node>::DynamicSafeArrayDim2
                            (const DynamicSafeArrayDim2 &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new DynamicSafeArrayDim1<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicSafeArrayDim1<Node>(*System.m[i]);
    }
}

template <typename Node>
DynamicSafeArrayDim2<Node>::DynamicSafeArrayDim2(int N1, int N2) {
    if (N1 <= 0 || N2 <= 0) {
        throw(DynamicSafeArrayDim2_InitError_UNSUPPORTEDSIZE());
    }

    arr_size = N1;
    m = new DynamicSafeArrayDim1<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicSafeArrayDim1<Node>(N2);
    }
}

template <typename Node>
DynamicSafeArrayDim2<Node>::~DynamicSafeArrayDim2() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node>
DynamicSafeArrayDim1<Node>& DynamicSafeArrayDim2<Node>::operator[] (int k) {
    if (k < 0 || k >= arr_size) {
        throw(DynamicSafeArrayDim2_WrongIndex());
    }
    return *m[k];
}





template <typename Node>
DynamicSafeArrayDim3<Node>::DynamicSafeArrayDim3() {
    arr_size = 0;
    m = nullptr;
}

template <typename Node>
DynamicSafeArrayDim3<Node>::DynamicSafeArrayDim3
                            (const DynamicSafeArrayDim3 &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new DynamicSafeArrayDim2<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicSafeArrayDim2<Node>(*System.m[i]);
    }
}

template <typename Node>
DynamicSafeArrayDim3<Node>::DynamicSafeArrayDim3(int N1, int N2, int N3) {
    if (N1 <= 0 || N2 <= 0 || N3 <= 0) {
        throw(DynamicSafeArrayDim3_InitError_UNSUPPORTEDSIZE());
    }

    arr_size = N1;
    m = new DynamicSafeArrayDim2<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicSafeArrayDim2<Node>(N2, N3);
    }
}

template <typename Node>
DynamicSafeArrayDim3<Node>::~DynamicSafeArrayDim3() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node>
DynamicSafeArrayDim2<Node>& DynamicSafeArrayDim3<Node>::operator[] (int k) {
    if (k < 0 || k >= arr_size) {
        throw(DynamicSafeArrayDim3_WrongIndex());
    }
    return *m[k];
}






template <typename Node>
DynamicSafeArrayDim4<Node>::DynamicSafeArrayDim4() {
    arr_size = 0;
    m = nullptr;
}

template <typename Node>
DynamicSafeArrayDim4<Node>::DynamicSafeArrayDim4
                            (const DynamicSafeArrayDim4 &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new DynamicSafeArrayDim3<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicSafeArrayDim3<Node>(*System.m[i]);
    }
}

template <typename Node>
DynamicSafeArrayDim4<Node>::DynamicSafeArrayDim4(int N1, int N2, int N3, int N4) {
    if (N1 <= 0 || N2 <= 0 || N3 <= 0 || N4 <= 0) {
        throw(DynamicSafeArrayDim4_InitError_UNSUPPORTEDSIZE());
    }

    arr_size = N1;
    m = new DynamicSafeArrayDim3<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicSafeArrayDim3<Node>(N2, N3, N4);
    }
}

template <typename Node>
DynamicSafeArrayDim4<Node>::~DynamicSafeArrayDim4() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node>
DynamicSafeArrayDim3<Node> &DynamicSafeArrayDim4<Node>::operator[] (int k) {
    if (k < 0 || k >= arr_size) {
        throw(DynamicSafeArrayDim4_WrongIndex());
    }
    return *m[k];
}









































template <typename Node>
class DynamicUnsafeArrayDim1
{
    private:
        int arr_size;
    public:
        DynamicUnsafeArrayDim1(int N);
        DynamicUnsafeArrayDim1();
        DynamicUnsafeArrayDim1(const DynamicUnsafeArrayDim1 &System);
        ~DynamicUnsafeArrayDim1();
        Node **m;

        Node& operator[] (int k);
        const Node& operator[] (int k) const;
};

template <typename Node>
class DynamicUnsafeArrayDim2
{
    private:
        int arr_size;
    public:
        DynamicUnsafeArrayDim2(int N1, int N2);
        DynamicUnsafeArrayDim2();
        DynamicUnsafeArrayDim2(const DynamicUnsafeArrayDim2 &System);
        ~DynamicUnsafeArrayDim2();
        DynamicUnsafeArrayDim1<Node> **m;

        DynamicUnsafeArrayDim1<Node>& operator[] (int k);
        const DynamicUnsafeArrayDim1<Node>& operator[] (int k) const;
};

template <typename Node>
class DynamicUnsafeArrayDim3
{
    private:
        int arr_size;
    public:
        DynamicUnsafeArrayDim3(int N1, int N2, int N3);
        DynamicUnsafeArrayDim3();
        DynamicUnsafeArrayDim3(const DynamicUnsafeArrayDim3 &System);
        ~DynamicUnsafeArrayDim3();
        DynamicUnsafeArrayDim2<Node> **m;

        DynamicUnsafeArrayDim2<Node>& operator[] (int k);
        const DynamicUnsafeArrayDim2<Node>& operator[] (int k) const;
};



template <typename Node>
class DynamicUnsafeArrayDim4
{
    private:
        int arr_size;
    public:
        DynamicUnsafeArrayDim4(int N1, int N2, int N3, int N4);
        DynamicUnsafeArrayDim4();
        DynamicUnsafeArrayDim4(const DynamicUnsafeArrayDim4 &System);
        ~DynamicUnsafeArrayDim4();
        DynamicUnsafeArrayDim3<Node> **m;

        DynamicUnsafeArrayDim3<Node>& operator[] (int k);
        const DynamicUnsafeArrayDim3<Node>& operator[] (int k) const;
};




class DynamicUnsafeArrayDim1_InitError_UNSUPPORTEDSIZE
{
    public:
        DynamicUnsafeArrayDim1_InitError_UNSUPPORTEDSIZE() {}
    private:
        int k = 0;
};

class DynamicUnsafeArrayDim2_InitError_UNSUPPORTEDSIZE
{
    public:
        DynamicUnsafeArrayDim2_InitError_UNSUPPORTEDSIZE() {}
    private:
        int k = 0;
};

class DynamicUnsafeArrayDim3_InitError_UNSUPPORTEDSIZE
{
    public:
        DynamicUnsafeArrayDim3_InitError_UNSUPPORTEDSIZE() {}
    private:
        int k = 0;
};

class DynamicUnsafeArrayDim4_InitError_UNSUPPORTEDSIZE
{
    public:
        DynamicUnsafeArrayDim4_InitError_UNSUPPORTEDSIZE() {}
    private:
        int k = 0;
};










template <typename Node>
DynamicUnsafeArrayDim1<Node>::DynamicUnsafeArrayDim1()
:arr_size(0)
{
    m = nullptr;
}


template <typename Node>
DynamicUnsafeArrayDim1<Node>::DynamicUnsafeArrayDim1
                            (const DynamicUnsafeArrayDim1 &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new Node*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new Node(*System.m[i]);
    }
}

template <typename Node>
DynamicUnsafeArrayDim1<Node>::DynamicUnsafeArrayDim1(int N)
:arr_size(N)
{
    if (arr_size <= 0) {
        throw(DynamicUnsafeArrayDim1_InitError_UNSUPPORTEDSIZE());
    }

    m = new Node*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new Node;
    }
}

template <typename Node>
DynamicUnsafeArrayDim1<Node>::~DynamicUnsafeArrayDim1() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node>
Node& DynamicUnsafeArrayDim1<Node>::operator[] (int k) {
    return *m[k];
}

template <typename Node>
const Node& DynamicUnsafeArrayDim1<Node>::operator[] (int k) const {
    return *m[k];
}






template <typename Node>
DynamicUnsafeArrayDim2<Node>::DynamicUnsafeArrayDim2() {
    arr_size = 0;
    m = nullptr;
}

template <typename Node>
DynamicUnsafeArrayDim2<Node>::DynamicUnsafeArrayDim2
                            (const DynamicUnsafeArrayDim2 &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new DynamicUnsafeArrayDim1<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicUnsafeArrayDim1<Node>(*System.m[i]);
    }
}

template <typename Node>
DynamicUnsafeArrayDim2<Node>::DynamicUnsafeArrayDim2(int N1, int N2) {
    if (N1 <= 0 || N2 <= 0) {
        throw(DynamicUnsafeArrayDim2_InitError_UNSUPPORTEDSIZE());
    }

    arr_size = N1;
    m = new DynamicUnsafeArrayDim1<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicUnsafeArrayDim1<Node>(N2);
    }
}

template <typename Node>
DynamicUnsafeArrayDim2<Node>::~DynamicUnsafeArrayDim2() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node>
DynamicUnsafeArrayDim1<Node>& DynamicUnsafeArrayDim2<Node>::operator[] (int k) {
    return *m[k];
}

template <typename Node>
const DynamicUnsafeArrayDim1<Node>& DynamicUnsafeArrayDim2<Node>::operator[] (int k) const {
    return *m[k];
}





template <typename Node>
DynamicUnsafeArrayDim3<Node>::DynamicUnsafeArrayDim3() {
    arr_size = 0;
    m = nullptr;
}

template <typename Node>
DynamicUnsafeArrayDim3<Node>::DynamicUnsafeArrayDim3
                            (const DynamicUnsafeArrayDim3 &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new DynamicUnsafeArrayDim2<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicUnsafeArrayDim2<Node>(*System.m[i]);
    }
}

template <typename Node>
DynamicUnsafeArrayDim3<Node>::DynamicUnsafeArrayDim3(int N1, int N2, int N3) {
    if (N1 <= 0 || N2 <= 0 || N3 <= 0) {
        throw(DynamicUnsafeArrayDim3_InitError_UNSUPPORTEDSIZE());
    }

    arr_size = N1;
    m = new DynamicUnsafeArrayDim2<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicUnsafeArrayDim2<Node>(N2, N3);
    }
}

template <typename Node>
DynamicUnsafeArrayDim3<Node>::~DynamicUnsafeArrayDim3() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node>
DynamicUnsafeArrayDim2<Node>& DynamicUnsafeArrayDim3<Node>::operator[] (int k) {
    return *m[k];
}

template <typename Node>
const DynamicUnsafeArrayDim2<Node>& DynamicUnsafeArrayDim3<Node>::operator[] (int k) const {
    return *m[k];
}






template <typename Node>
DynamicUnsafeArrayDim4<Node>::DynamicUnsafeArrayDim4() {
    arr_size = 0;
    m = nullptr;
}

template <typename Node>
DynamicUnsafeArrayDim4<Node>::DynamicUnsafeArrayDim4
                            (const DynamicUnsafeArrayDim4 &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new DynamicUnsafeArrayDim3<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicUnsafeArrayDim3<Node>(*System.m[i]);
    }
}

template <typename Node>
DynamicUnsafeArrayDim4<Node>::DynamicUnsafeArrayDim4(int N1, int N2, int N3, int N4) {
    if (N1 <= 0 || N2 <= 0 || N3 <= 0 || N4 <= 0) {
        throw(DynamicUnsafeArrayDim4_InitError_UNSUPPORTEDSIZE());
    }

    arr_size = N1;
    m = new DynamicUnsafeArrayDim3<Node>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicUnsafeArrayDim3<Node>(N2, N3, N4);
    }
}

template <typename Node>
DynamicUnsafeArrayDim4<Node>::~DynamicUnsafeArrayDim4() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node>
DynamicUnsafeArrayDim3<Node> &DynamicUnsafeArrayDim4<Node>::operator[] (int k) {
    return *m[k];
}

template <typename Node>
const DynamicUnsafeArrayDim3<Node> &DynamicUnsafeArrayDim4<Node>::operator[] (int k) const {
    return *m[k];
}


































class DynamicCycledArray_InitError_UNSUPPORTEDSIZE
{
    public:
        int k = 0;
        DynamicCycledArray_InitError_UNSUPPORTEDSIZE(int Size) {k = Size;}
};


class DynamicCycledArray_ZERODIMENSION
{
    public:
        int k = 0;
        DynamicCycledArray_ZERODIMENSION() {}
};








template <typename Node, unsigned int Dimension>
class DynamicCycledArray
{
    private:
        int arr_size;
    public:
        DynamicCycledArray();
        DynamicCycledArray(int Size);
        DynamicCycledArray(const DynamicCycledArray<Node, Dimension> &System);
        ~DynamicCycledArray();

        DynamicCycledArray<Node, Dimension-1> **m = nullptr;

        DynamicCycledArray<Node, Dimension-1> &operator[] (int k);
};


template <typename Node, unsigned int Dimension>
DynamicCycledArray<Node, Dimension>::DynamicCycledArray() {
    m = nullptr;
    arr_size = 0;
}

template <typename Node, unsigned int Dimension>
DynamicCycledArray<Node, Dimension>::DynamicCycledArray
                            (const DynamicCycledArray<Node, Dimension> &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new DynamicCycledArray<Node, Dimension-1>*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new DynamicCycledArray<Node, Dimension-1>(System[i]);
    }
}

template <typename Node, unsigned int Dimension>
DynamicCycledArray<Node, Dimension>::DynamicCycledArray(int Size) {
    if (Size <= 0) {
        throw(DynamicCycledArray_InitError_UNSUPPORTEDSIZE(Size));
    }

    arr_size = Size;
    m = new DynamicCycledArray<Node, Dimension-1>*[Size];
    for (int i = 0; i < Size; i++) {
        m[i] = new DynamicCycledArray<Node, Dimension-1>(Size);
    }
}

template <typename Node, unsigned int Dimension>
DynamicCycledArray<Node, Dimension>::~DynamicCycledArray() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node, unsigned int Dimension>
DynamicCycledArray<Node, Dimension-1> &DynamicCycledArray<Node, Dimension>::operator[] (int k) {
    return *m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}








template <typename Node>
class DynamicCycledArray<Node, 1>
{
    private:
        int arr_size;
    public:
        DynamicCycledArray();
        DynamicCycledArray(int Size);
        DynamicCycledArray(const DynamicCycledArray<Node, 1> &System);
        ~DynamicCycledArray();

        Node **m = nullptr;

        Node &operator[] (int k);
};


template <typename Node>
DynamicCycledArray<Node, 1>::DynamicCycledArray() {
    m = nullptr;
    arr_size = 0;
}

template <typename Node>
DynamicCycledArray<Node, 1>::DynamicCycledArray
                            (const DynamicCycledArray<Node, 1> &System)
:arr_size(System.arr_size)
{
    if (System.m == nullptr) {
        m = nullptr;
        return;
    }

    m = new Node*[arr_size];
    for (int i = 0; i < arr_size; i++) {
        m[i] = new Node(System[i]);
    }
}

template <typename Node>
DynamicCycledArray<Node, 1>::DynamicCycledArray(int Size) {
    if (Size <= 0) {
        throw(DynamicCycledArray_InitError_UNSUPPORTEDSIZE(Size));
    }

    arr_size = Size;
    m = new Node*[Size];
    for (int i = 0; i < Size; i++) {
        m[i] = new Node;
    }
}

template <typename Node>
DynamicCycledArray<Node, 1>::~DynamicCycledArray() {
    for (int i = 0; i < arr_size; i++) {
        delete m[i];
    }
    delete [] m;
}


template <typename Node>
Node &DynamicCycledArray<Node, 1>::operator[] (int k) {
    return *m[(k >= 0 ? k%arr_size : arr_size + k%arr_size)%arr_size];
}





template <typename Node>
class DynamicCycledArray<Node, 0>
{
    public:
        DynamicCycledArray() {
            throw(DynamicCycledArray_ZERODIMENSION());
        }
};





#endif // ARRAYDIM4_H
