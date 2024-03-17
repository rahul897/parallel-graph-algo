#include <stdio.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
using namespace std;
//use initialise_matrix to create space in file for r rows and c columns with each element occupying s space
//use enter_element() to enter an element in matrix at rth row and cth column directly in file
//use enter_row() to enter row of elements or partial row of elements from start to end in matrix present in file
//use enter_column() to enter column of elements or partial column of elements from start to end in matrix present in file
//use fetch_element() to fetch an element in matrix at rth row and cth column directly from file
//use fetch_row() to fetch row of elements or partial row of elements from start to end in matrix present in file into a vector
//use fetch_column() to fetch column of elements or partial column of elements from start to end in matrix present in file into a vector
//use print_vector() to print a vector

//addition  of vectors and matrices is defined
//multiplication of vector and matrix is defined
//outer product of multiplication of two vectors is defined
//product of two matrices is defined

template <class T> class MATRIX
{
    public:

    int nr,nc,s;
    //fstream fin;
    MATRIX()
    {
    }
    MATRIX(int nr,int nc,int s)
    {
        this->nr=nr;
        this->nc=nc;
        this->s=s;
      //  this->fin=fin;
    }
    void initialise_matrix(fstream &fin)
    {
        for(int i=0;i<nr;i++)
        {
            for(int j=0;j<nc;j++)
            {
                for(int i=0;i<s;i++)
                fin<<" ";
            //fin<<"|";
            }
            fin<<"\n";
        }
    }
    void initialise_matrix(fstream &fin,int a)
    {
        for(int i=0;i<nr;i++)
        {
            for(int j=0;j<nc;j++)
            {
                fin<<a;
                for(int i=0;i<s-1;i++)
                fin<<" ";
            //fin<<"|";
            }
            fin<<"\n";
        }
    }

    T fetch_element(int r,int c,fstream &fin)
    {
        fin.clear();
        fin.seekg((r*nc+c)*s+r*1,ios::beg);
        T num;
        fin>>num;

       // cout<<num<<" num \n";
        return num;
    }
    void print_vector(vector<T> vect,fstream &fin)
    {
       // printf("\n");
        fin<<"\n";
        for(int i=0;i<vect.size();i++)
            fin<<vect[i]<<"  ";
        fin<<"\n";
       // printf("\n");
    }

    void print_matrix(vector<vector<T> > res,fstream &fin)
    {
        //printf("\n");
            fin<<"\n";
        for(int i=0;i<res.size();i++)
            print_vector(res[i],fin);
            fin<<"\n";// printf("\n");
    }

    void fetch_row(vector<T> &vec,int r,int start,int endc,fstream &fin)
    {
        vec.clear();
        for(int i=start;i<=endc;i++)
        {
            T num=fetch_element(r,i,fin);
            vec.push_back(num);
        }
    }
    void fetch_column(vector<T> &vec,int c,int start,int endc,fstream &fin)
    {
        vec.clear();
        for(int i=start;i<endc;i++)
        {
            T num=fetch_element(i,c,fin);
            vec.push_back(num);
        }
    }
    void fetch_diag(vector<T> &vec,fstream &fin)
    {
        unsigned long long int i,j;
        vec.clear();
        for(i=0;i<nr;i++)
        {
                T num;
                num=fetch_element(i,i,fin);
                vec.push_back(num);
        }
    }
    void enter_element(int r,int c,T num,fstream &fin)
    {
        fin.clear();
        fin.seekg((r*nc+c)*s+r*1,ios::beg);
        for(int i=0;i<s;i++)
            fin<<" ";
        fin.clear();
        fin.seekg((r*nc+c)*s+r*1,ios::beg);
        fin<<num;
    }
    int enter_row(vector<T> vec,int r,int start,int endc,fstream &fin)
    {
        if(vec.size()<endc-start+1)
            return 0;

        for(int i=start,j=0;i<=endc;i++,j++)
        {
            enter_element(r,i,vec[j],fin);
        }
        return 1;
    }
    int enter_column(vector<T> vec,int c,int start,int endc,fstream &fin)
    {
        if(vec.size()<endc-start+1)
            return 0;

        for(int i=start,j=0;i<endc;i++,j++)
        {
            enter_element(i,c,vec[j],fin);
        }
        return 1;
    }
    void fetch_matrix(vector<vector<T> > &res,int r1,int r2,int c1,int c2,fstream &fin)
    {
        res.clear();
        for(int i=r1;i<=r2;i++)
        {
            vector<T> get_t;

            fetch_row(get_t,i,c1,c2,fin);

            res.push_back(get_t);

        }
    }
    void enter_matrix(vector<vector<T> > res,int r1,int r2,int c1,int c2,fstream &fin)
    {
        int l=r2-r1+1;
       // cout<<"l= "<<l<<"\n";
        for(int i=r1,j=0;i<=r2&&j<l;i++,j++)
        {
            vector<T> get_t=res[j];

            enter_row(get_t,i,c1,c2,fin);

        }
    }
    void add_vector(vector<double> &res,vector<double> a,vector<double> b)
    {
        for(int i=0;i<a.size();i++)
            res.push_back(a[i]+b[i]);
    }
    void add_matrix(vector<vector<T> > &res,vector<vector<T> > a,vector<vector<T> > b)
    {
        for(int i=0;i<a.size();i++)
        {
            vector<T> temp;
            add_vector(temp,a[i],b[i]);
            res[i]=temp;
        }
    }

    void mult_2_mat(vector<T> A,vector<T> B,vector<vector<T> > &res,double val)
    {
        res.resize(A.size());
        for(int i=0;i<A.size();i++)
        {
            res[i].resize(B.size());
            for(int j=0;j<B.size();j++)
            {
                double sum=0;
                sum+=A[i]*B[j];
                res[i][j]=(val*sum);
            }
        }
    }

    void mult_vect_mat(fstream &fin,int nr,vector<T> B,vector<T> &res)
    {
        res.resize(nr);
        vector<T> A;
        for(int i=0;i<nr;i++)
        {
                A.clear();
                fetch_row(A,i,0,nr-1,fin);
                double sum=0;
                for(int k=0;k<A.size();k++)
                    sum+=A[k]*B[k];
                res[i]=(sum);
              //  cout<<"a alal"<<sum<<"  ";

        }
        //cout<<"\n";
    }

    void mult_matrices(vector<vector<T> > A,vector<vector<T> > B,vector<vector<T> > &res,double val)
    {
        res.resize(A.size());
        for(int i=0;i<A.size();i++)
        {
            res[i].resize(B[0].size());
            for(int j=0;j<B[0].size();j++)
            {
                double sum=0;
                for(int k=0;k<B.size();k++)
                    sum+=A[i][k]*B[k][j];
                res[i][j]=(val*sum);
            }
        }
    }
};
