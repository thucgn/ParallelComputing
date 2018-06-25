/*************************************************************************
	> File Name: sparse_matrix_csr.h
	> Author: CGN
	> Mail: 
	> Created Time: 二  6/12 16:43:50 2018
 ************************************************************************/

#ifndef _SPARSE_MATRIX_CSR_H
#define _SPARSE_MATRIX_CSR_H

#include <string>
#include <random>
#include <functional>
#include <fstream>
#include <iostream>
#include "util.h"


//sparse matrix in CSR format
template <typename T>
class SparseMatrix
{
public:
    int nrow; 
    int ncol; 
    int nnz; //number of non zeor

    T* values;
    int* col_idx;
    int* row_ptr;

    SparseMatrix():nrow(0), ncol(0), nnz(0), values(nullptr), col_idx(nullptr), row_ptr(nullptr) {}

    ~SparseMatrix();

    void read_csr(std::string filename);
    void print();

};

template <typename T>
void SparseMatrix<T>::print()
{
    LOG("row: %d, col: %d, nnz: %d", nrow, ncol, nnz);
    for(int i = 0;i < nnz; i ++)
        std::cout << values[i] << " ";
    std::cout << std::endl;

    for(int i = 0;i < nnz; i ++)
        std::cout << col_idx[i] << " ";
    std::cout << std::endl;

    for(int i = 0;i < nrow+1; i++)
        std::cout << row_ptr[i] << " ";
    std::cout << std::endl;

}

template <typename T>
void SparseMatrix<T>::read_csr(std::string filename)
{
    std::ifstream fin(filename);
    CHECK(fin.is_open());

    fin >> nrow;
    fin >> ncol;
    fin >> nnz;

    values = new T[nnz];
    col_idx = new int[nnz];
    row_ptr = new int[nrow+1];

    for(int i = 0;i < nnz; i++)
        fin >> values[i];
    for(int i = 0;i < nnz; i++)
        fin >> col_idx[i];
    for(int i = 0;i < nrow+1; i++)
        fin >> row_ptr[i];


    fin.close();
}

template <typename T>
SparseMatrix<T>::~SparseMatrix()
{
    if(values)
        delete [] values;
    if(col_idx)
        delete [] col_idx;
    if(row_ptr)
        delete [] row_ptr;
}

void generate_conjugate_sparsematrix(std::string filename, int row, int col, int sim_width)
{
    CHECK(row == col);

    int nnz = row;
    for(int k = row-1; k >= row-sim_width; k --)
        nnz += 2*k;

    double* values = new double[nnz];
    int* col_idx = new int[nnz];
    int* row_ptr = new int[row+1];

    int v_cnt = 0, c_cnt = 0, r_cnt = 0;

    //generate data
    std::random_device rd;
    std::minstd_rand gen(rd());
    std::uniform_int_distribution<> dis(5000,10000);
    auto dice = std::bind(dis, gen);
    double rand_value = dice() * 1e-4;

    for(int i = 0;i < row; i ++)
    {
        row_ptr[r_cnt ++] = v_cnt;
        for(int j = i-sim_width; j <= i+sim_width; j++)
        {
            if(j >= 0 && j < col)
            {
                if(i == j)
                {
                    values[v_cnt ++] = dice() * 1e-4 + 0.5;
                    col_idx[c_cnt ++] = j;
                }
                else 
                {
                    values[v_cnt ++] = rand_value + (i+j + 1.0)/(i*j + 1.0);
                    col_idx[c_cnt ++] = j;
                }
            }
        }
    }
    CHECK(v_cnt == nnz);
    row_ptr[r_cnt] = v_cnt;

    //write to file
    std::ofstream fout(filename); 
    CHECK(fout.is_open());
    fout << row << " " << col  << " "  << nnz << " ";
    for(int i = 0;i < nnz; i ++)
        fout << values[i] << ' ';
    for(int i = 0;i < nnz; i ++)
        fout << col_idx[i] << ' ';
    for(int i = 0;i < row+1; i++)
        fout << row_ptr[i] << ' ';

    fout.close();

    LOG("write matrix to %s, row: %d, col:%d, nnz: %d", filename.c_str(), row,col,nnz);

    delete [] values;
    delete [] col_idx;
    delete [] row_ptr;
}

#endif

