/*************************************************************************
	> File Name: conjugate_gradient.cpp
	> Author: CGN
	> Mail: 
	> Created Time: 二  6/12 15:43:51 2018
 ************************************************************************/

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <omp.h>
#include <random>
#include <functional>
#include "util.h"
#include "sparse_matrix_csr.h"


#ifdef DEBUG
const bool FLAG = true;
#else
const bool FLAG = false;
#endif

using namespace std;

template <typename T>
class Vec
{
public:
    bool is_row;
    int length;
    T* data;
    Vec(bool is_row, int length)
    {
        this->is_row = is_row;
        this->length = length;
        data = new T[length];
    }

    T& operator[](int pos)  { return data[pos]; }

    void add(Vec& vec)
    {
        for(int i = 0;i < length; i++)
            data[i] += vec[i];
    }

    void sub(Vec& vec)
    {
        for(int i = 0;i < length; i++)
            data[i] -= vec[i];
    }

    void print()
    {
        for(int i = 0; i < length; i ++)
            cout << data[i] << " ";
        cout << endl;
    }

    ~Vec()
    {
        if(data)
            delete [] data;
    }
};


template <typename T>
void spmv(SparseMatrix<T>& sm, Vec<T>& vec, Vec<T>& ret_vec)
{
    T* dv = vec.data;
    T* drv = ret_vec.data;
    int wb = (sm.nrow + THREADS - 1) / THREADS;
    int rb = sm.nrow - (THREADS-1)*wb;
#pragma omp parallel for  schedule(static) num_threads(THREADS)
    for(int t = 0;t < THREADS; t ++)
    {
        int wc = (t != THREADS-1) ? wb : rb;
        T tmp;
        int i, j;
        for(i = t*wb; i < t*wb+wc; i++)
        {
            tmp = T();
            for(j = sm.row_ptr[i]; j < sm.row_ptr[i+1]; j++)
                tmp += sm.values[j] * dv[sm.col_idx[j]];
            drv[i] = tmp;
        }
    }
}

template <typename T>
void sub_vec(Vec<T>& a, Vec<T>& b, Vec<T>& ret)
{
    T* da = a.data;
    T* db = b.data;
    T* dret = ret.data;
#pragma omp parallel num_threads(THREADS)
    for(int i = 0;i < ret.length; i++)
        dret[i] = da[i] - db[i];  
}

template <typename T>
void add_vec(Vec<T>& a, Vec<T>& b, Vec<T>& ret) 
{
    T* da = a.data;
    T* db = b.data;
    T* dret = ret.data;
#pragma omp parallel for num_threads(THREADS)
    for(int i = 0;i < ret.length; i++)
      dret[i] = da[i] + db[i];  
};

template <typename T>
bool is_small_enough(Vec<T>& v, T delta)
{
    T tmp = T();
    T* data = v.data;
    for(int i = 0; i < v.length; i++)
        tmp += fabs(data[i]);
    return (tmp < delta);
};

template <typename T>
T mul_vec(Vec<T>& a, Vec<T>& b)
{
    T tmp = T();
    T* da = a.data;
    T* db = b.data;
    for(int i = 0;i < a.length; i++)
        tmp += da[i] * db[i];
    return tmp;
};

template <typename T>

void vec_mul_scale(Vec<T>& a, T alpha, Vec<T>& b)
{
    T* da = a.data;
    T* db = b.data;
#pragma omp parallel for num_threads(THREADS)
    for(int i = 0;i < b.length; i ++)
        db[i] = da[i] * alpha;
};


/*********************************
 * note: ret_vec should be initialized with a init-value;
 * func: x = A_reverse * b
 * *******************************/
template <typename T>
void conjugate_sparsematrix(SparseMatrix<T>& A, Vec<T>& b, Vec<T>& x, T delta)
{

    CHECK(!b.is_row && !x.is_row);
    CHECK(A.ncol == b.length);
    CHECK(A.ncol == x.length);


    Vec<T> tb(false, b.length);
    Vec<T> r(false, x.length);
    Vec<T> p(false, r.length);
    Vec<T> tmp(false, p.length);


    spmv(A, x, tb);
    sub_vec(b, tb, r);

    for(int i = 0;i < p.length; i++)
        p[i] = r[i];

    int k = 0;

    T alpha = T();
    T beta = T();
    T rrk = T();


    do{
        rrk = mul_vec(r,r); // rk * rk
        spmv(A, p, tb); // tb = A*pk
        alpha = rrk / mul_vec(p, tb);

        vec_mul_scale(p, alpha, tmp); //tmp = alpha * pk
        x.add(tmp);
        
        vec_mul_scale(tb, alpha, tmp); //tmp = alpha * A * pk 
        r.sub(tmp); //rk+1

        if(is_small_enough(r, delta))
            break;

        beta = mul_vec(r, r) / rrk;
        vec_mul_scale(p, beta, tmp); //tmp = beta * pk;
        add_vec(r, tmp, p);
        k ++;
    }while(true);

    LOG("iteration times(k) : %d", k);

}

template <typename T>
void write_vec_to_file(string filename, Vec<T>& vec)
{
    int n = vec.length;
    //write to file
    std::ofstream fout(filename); 
    CHECK(fout.is_open());
    fout << n << ' ';
    for(int i = 0;i < n; i++)
        fout << vec[i] << ' ';
    fout.close();

}

template <typename T>
void read_vec_from_file(string filename, Vec<T>& vec)
{
	int n = vec.length;
	std::ifstream fin(filename);
	CHECK(fin.is_open());
	int pn;
	fin >> pn;
	CHECK(pn == n);
	for(int i = 0;i < pn; i ++)
		fin >> vec[i];
	fin.close();
}

int main(int argc, char** argv)
{



    if(argc == 2)
	{
    	int n = atoi(argv[1]);
    	SparseMatrix<double> sm;
    	Vec<double> x(false, n);
    	Vec<double> b(false, n);
    	Vec<double> x0(false, n);
        generate_conjugate_sparsematrix("matrix", n, n, 20);

    	sm.read_csr("matrix");

    	std::random_device rd;
    	std::minstd_rand gen(rd());
    	std::uniform_int_distribution<> dis(2000,18000);
    	auto dice = std::bind(dis, gen);
    	for(int i = 0;i < x.length; i++)
        	x[i] = 5.0 + (dice() - 10000)*1e-4;

    	spmv(sm, x, b);

    	for(int i = 0;i < x0.length; i++)
        	x0[i] = 5.0;

    	TIME_T ts, tt;
    	MARK_TIME(ts);
    	conjugate_sparsematrix<double>(sm, b, x0, 0.001 * n);
    	MARK_TIME(tt);
    	LOG("time %.5f", DIFF_TIME(tt, ts));

    	double diff = 0.0;
    	for(int i = 0;i < x0.length; i ++)
    	    diff += fabs(x0[i] - x[i]);
    	//write_vec_to_file("x_ret", x0);
    	LOG("difference between x and x0: %lf", diff);
	}
	else if(argc == 4)
	{
    	SparseMatrix<double> sm;

		sm.read_csr(argv[1]);

		int n = sm.nrow;

    	Vec<double> x(false, n);
    	Vec<double> b(false, n);
    	Vec<double> x0(false, n);

		read_vec_from_file<double>(argv[2], x);
		spmv(sm, x, b);

		read_vec_from_file<double>(argv[3], x0);

    	TIME_T ts, tt;
    	MARK_TIME(ts);
    	conjugate_sparsematrix<double>(sm, b, x0, 0.001 * n);
    	MARK_TIME(tt);
    	LOG("time %.5f", DIFF_TIME(tt, ts));

    	double diff = 0.0;
    	for(int i = 0;i < x0.length; i ++)
    	    diff += fabs(x0[i] - x[i]);
    	//write_vec_to_file("x_ret", x0);
    	LOG("difference between x and x0: %lf", diff);
	}
	else 
	{
		LOG("输入参数错误, 比如 ./main1 40000，或者 ./main1 matrix x x0");
		exit(0);;
	}

    //write_vec_to_file("x", vec);


    //write_vec_to_file("b", ret);

    
    //write_vec_to_file("x0", x0);


    return 0;
}


