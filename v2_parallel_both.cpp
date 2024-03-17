//version2 Parallel algorithm
//Implementation of Parallel NSD and PAUL using Pthreads

#include <iostream>
#include <math.h>
#include <fstream>
#include <set>
#include <vector>
#include <utility>
#include <unistd.h>
#include <stdlib.h>
#include "temp_matrix.h"
#include "gensvd.h"
using namespace std;

//Declaration of global buyers and global matching
set <pair<int,int> > s;                            
set <pair<int,int> >::iterator itr,itr1;
set <int> buyers;
set <int>::iterator itrbuyer;

#define MAT_SIZE 100
#define dummy_1d vector <double>
#define dummy_2d vector <vector<double > >

MATRIX<double> *a,*b,*atilde,*btilde,*u,*v,*sigmam,*w,*z,*Xi,*X;
int n,p,q,nr1,nr2,s1;
double factor,alp;
fstream finX,finoutX;
pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;

// Barrier variable
pthread_barrier_t barr;

struct param                                     // used for Nsd thread identification
{
    int r;
    int u;
    int s;
};

struct plocal                                    // used for updating local price
{
    double p;
    int buyer;
}*price;

/*
0 1 0 0 0 0 0 0 1 1
1 0 0 0 1 1 0 1 1 0
0 0 0 1 1 0 0 1 1 0
0 0 1 0 0 1 0 0 0 0
0 1 1 0 0 1 0 1 1 1
0 1 0 1 1 0 1 1 0 0
0 0 0 0 0 1 0 0 0 1
0 1 1 0 1 1 0 0 0 0
1 1 1 0 1 0 0 0 0 0
1 0 0 0 1 0 1 0 0 0

0 1 0 0 0 0 0 0 1 1
1 0 0 0 1 1 0 1 1 0
0 0 0 1 1 0 0 1 1 0
0 0 1 0 0 1 0 0 0 0
0 1 1 0 0 1 0 1 1 1
0 1 0 1 1 0 1 1 0 0
0 0 0 0 0 1 0 0 0 1
0 1 1 0 1 1 0 0 0 0
1 1 1 0 1 0 0 0 0 0
1 0 0 0 1 0 1 0 0 0
*/

//generate_matrix reads from the file the matrices into a 2-D Vector
void generate_matrix(dummy_2d &A,int a,int c,int b,fstream &fin)
{
    for(int i=0;i<a;i++)
        for(int j=0;j<c;j++)
        {
            double d;
            fin>>d;
            A[i].push_back(d);
        }
}

//computetilde computes the tilde of a given matrix and size
//The matrix is fetched row wise and row sum is computed
//The tilde matrix aij is aji/(sum_of_elements in row j)

void  computtilde(int nr1,int nc1,MATRIX<double> *m,MATRIX<double> *ati,fstream &fin,fstream &fin1)
{
    int rv[MAT_SIZE];
    for(int i=0;i<nr1;i++)
    {
        vector<double> row;
		m->fetch_row(row,i,0,nc1-1,fin);            // row is the vector which is used to fetch result   
        long long int sum=0;                        // ith row , 0 and nc1 represents starting and ending coulmn
        for(int j=0;j<nc1;j++)
            sum+=row[j];
        rv[i]=sum;
    }
    for(int i=0;i<nr1;i++)
        cout<<rv[i]<<" ";
    cout<<"\n\n";
    for(int i=0;i<nr1;i++)
    {
        for(int j=0;j<nc1;j++)
        {
            vector<double> vec;
			m->fetch_row(vec,j,0,nc1-1,fin);
            if(rv[j]!=0)
                ati->enter_element(i,j,vec[i]/rv[j],fin1);
            else
                ati->enter_element(i,j,0,fin1);
        }
    }
}

//generate_SVD generates the SVD of a given matrix
//It generates the left singular vector,right singular vector and singular values matrix
void generateSVD(dummy_2d &u,dummy_2d &v,dummy_2d &sigma,int nr1,int nr2)
{
    fstream finsvd;
    finsvd.open("svdf.txt");
    generate_matrix(u,nr1,nr1,0,finsvd);
    generate_matrix(v,nr2,nr2,0,finsvd);
    generate_matrix(sigma,nr1,nr2,0,finsvd);
}

//serial matching implementation of bipartite matching
void matching(dummy_2d X,int nr1,int nr2)
{
    for(int i=0;i<nr1;i++)
        buyers.insert(i);
    itrbuyer=buyers.begin();
    double *price=new double[nr2];
    for(int i=0;i<nr2;i++)
        price[i]=0;
    double elp=0.5;
    while(!buyers.empty())
    {
        int j=0,sj=0,k=0;double maxim=0;
        int i=*itrbuyer;
        maxim=0;
        for(k=0;k<nr2;k++)
        {
            if(X[i][k]-price[k]>=maxim)
            {
                maxim=X[i][k]-price[k];
                j=k;
            }
        }

        maxim=0;
        for(k=0;k<nr2;k++)
        {
            if(X[i][k]-price[k]>=maxim&&k!=j)
            {
                maxim=X[i][k]-price[k];
                sj=k;
            }
        }

        double u=X[i][j]-price[j];
        if(u<0)
        {
            buyers.erase(*itrbuyer);
            itrbuyer++;
            if(itrbuyer==buyers.end())
                itrbuyer=buyers.begin();
            continue;
        }
        double v=X[i][sj]-price[sj];
        price[j]=price[j]+u-v+elp;
        s.insert(make_pair(i,j));
        buyers.erase(*itrbuyer);
        for(int k=0;k<nr1;k++)
        if((itr1=s.find(make_pair(k,j)))!=s.end()&&k!=i)
        {
            s.erase(itr1);
            buyers.insert(k);
        }
        itrbuyer++;
        if(itrbuyer==buyers.end())
            itrbuyer=buyers.begin();
    }
}

//The function to print the global matching set

void print()
{
    itr1=s.begin();
    for(int i=0;i<s.size();i++)
    {
        cout<<"vertex "<<(*itr1).first+1<<" of first graph  is most similiar to  vertex "<<(*itr1).second+1<<" of second graph\n";
        itr1++;
    }
    cout<<"\n";
}

//The function to print matrix a having r rows and c columns
void printmatrix(dummy_2d a,int r,int c)
{
    cout<<"\n";
    for(int i=0;i<r;i++)
    {
        for(int j=0;j<c;j++)
            cout<<a[i][j]<<" ";
        cout<<"\n";
    }
    cout<<"\n";
}

//The function to multiply a matrix and vector and result is vector

void mult_vect_mat(dummy_2d A,dummy_1d B,dummy_1d &res)
{
    res.resize(A.size());
    for(int i=0;i<A.size();i++)
    {
            double sum=0;
            for(int k=0;k<A[0].size();k++)
                sum+=A[i][k]*B[k];
            res[i]=(sum);

    }
}

//The function to find the outer product of two vectors and result is a matrix
void mult_2_mat(dummy_1d A,dummy_1d B,dummy_2d &res,double val)
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
//The function to add two matrices
void mat_add(dummy_2d A,dummy_2d B,dummy_2d &res)
{
    for(int i=0;i<A.size();i++)
    {
        res[i].resize(A[0].size());
        for(int j=0;j<A[0].size();j++)
        {
                res[i][j]=A[i][j]+B[i][j];
		}
    }
}
//To fetch from file and add the partial matrix and save the result
void addmatrix(fstream &finx,vector<vector<double> > vec,int pos,int s1,int e1,int s2,int e2)
{
    int p=pos*nr1;
    int l,m;
    for(int i=s1,l=0;i<=e1;i++,l++)
    {
        for(int j=s2,m=0;j<=e2;j++,m++)
        {
            double a=Xi->fetch_element(p+i,j,finx);
            double l1=a+vec[l][m];
            Xi->enter_element(p+i,j,l1,finx);
        }
    }
}

//The Thread function to calculate the similarity matrix
//There are r*u threads each thread computes r*factor to (r+1)*factor-1 rows and  u*factor to (u+1)*factor-1 columns where factor is
//distribution to each thread
void *part_sim(void *arg)
{
    struct param * s1=(struct param *)arg;
    int r,u,s;
    double val;
    r=s1->r;
    u=s1->u;
    s=s1->s;
	int startu=r*factor;
    int endu=(r+1)*factor-1;
    if(endu>=nr1)
        endu=nr1-1;
    int startv=u*factor;
    int endv=(u+1)*factor-1;
    if(endv>=nr2)
        endv=nr2-1;

    //wmatrix and zmatrix contains vector iterates for all i's and k's within i's wik is a row in ith block kth row
    fstream finw,finz;
    finw.open("wmatrix.txt",ios::in|ios::out);
    finz.open("zmatrix.txt",ios::in|ios::out);

    vector<double> wi,zi;
    for(int i=0;i<s;i++)
    {
        for(int k=0;k<=n;k++)
        {
                val=pow(alp,k);
                dummy_2d temp;
                temp.resize(nr1);

                wi.clear();
                w->fetch_row(wi,i*(n+1)+k,startu,endu,finw);
                zi.clear();
                z->fetch_row(zi,i*(n+1)+k,startv,endv,finz);
                mult_2_mat(wi,zi,temp,val);
                if(k<=n-1)
                {
                    if(k==0)
                    {
                        pthread_mutex_lock(&mutex);
                        Xi->enter_matrix(temp,i*nr1+startu,(i)*nr1+endu,startv,endv,finX);
                        pthread_mutex_unlock(&mutex);
                    }
                    else
                    {

                        pthread_mutex_lock(&mutex);
                        addmatrix(finX,temp,i,startu,endu,startv,endv);
                        pthread_mutex_unlock(&mutex);

                    }
                }
                else
                if(k==n)
                {
                        pthread_mutex_lock(&mutex);
                        dummy_2d temp1;
                        temp1.resize(nr1);

                        Xi->fetch_matrix(temp1,i*nr1+startu,(i)*nr1+endu,startv,endv,finX);

                        for(int i1=0;i1<temp1.size();i1++)
                            for(int j1=0;j1<temp1[0].size();j1++)
                                temp1[i1][j1]=(1-alp)*temp1[i1][j1];
                        Xi->enter_matrix(temp1,i*nr1+startu,(i)*nr1+endu,startv,endv,finX);
                        addmatrix(finX,temp,i,startu,endu,startv,endv);

                        pthread_mutex_unlock(&mutex);
                }
        }
    }
    dummy_2d res;
    res.resize(nr1);
    for(int i=0;i<s;i++)
    {
        dummy_2d temp1;

        pthread_mutex_lock(&mutex);
        Xi->fetch_matrix(temp1,i*nr1+startu,(i)*nr1+endu,startv,endv,finX);
        pthread_mutex_unlock(&mutex);

        if(i==0)
            res=temp1;
        else
        {
                dummy_2d t;t.resize(temp1.size());
               Xi->add_matrix(t,temp1,res);
               res=t;
        }
    }
    pthread_mutex_lock(&mutex);
    X->enter_matrix(res,startu,endu,startv,endv,finoutX);
    pthread_mutex_unlock(&mutex);
    finw.close();
    finz.close();
}
//This is to read from normal file and enter into dynamic file

void set_matrix(int nr1,int nc1,fstream &fin,fstream &finm,MATRIX<double> *m)
{
    for(int i=0;i<nr1;i++)
    {
        for(int j=0;j<nc1;j++)
        {
            double a1;
            fin>>a1;
            m->enter_element(i,j,a1,finm);
        }
    }
}
//The thread function for bipartite matching algorithm
void *match(void *arg)
{
    int id=*(int *)arg;
    id=id-1;
    int start=id*(factor);
    int endd=(id+1)*factor;
    if(endd>=nr1)
        endd=nr1;


    set <pair<int,int> > s_local;
    set <pair<int,int> >::iterator itr_local,itr1_local;
    set <int> buyers_local;
    set <int>::iterator itrbuyer_local;
    set <pair<int,int> >::iterator itr,itr1;

    for(int i=start;i<endd;i++)
        buyers_local.insert(i);
    itrbuyer_local=buyers_local.begin();
    struct plocal price_local[nr2];
    for(int i=0;i<nr2;i++)
    {
        price_local[i].p=0;
        price_local[i].buyer=0;
    }
    vector<vector<double> > X;
    MATRIX<double> *MX;
    MX=new MATRIX<double>(nr1,nr2,s1);
    pthread_mutex_lock(&mutex);
    MX->fetch_matrix(X,start,endd-1,0,nr2-1,finoutX);
    pthread_mutex_unlock(&mutex);
    double elp_local;
    elp_local=0.1;
  // printmatrix(X,endd-start,nr2);
    int l=1;
    while(!buyers.empty())
    {
        struct plocal inprice_local[nr2];
        for(int i=0;i<nr2;i++)
        {
            inprice_local[i].p=price_local[i].p;
            inprice_local[i].buyer=price_local[i].buyer;
        }
        int j=0,sj=0,k=0;double maxim=0;
        int sindex=start;
        if(buyers_local.size()>0)
          {
          itrbuyer_local=buyers_local.begin();

        }
        for(;itrbuyer_local!=buyers_local.end()&&buyers_local.size()>0;itrbuyer_local++)
        {
           int i=*itrbuyer_local;
            i=i-sindex;
            maxim=0;

            for(k=0;k<nr2;k++)
            {
                if((X[i][k]-price_local[k].p)>=maxim)
                {
                    maxim=X[i][k]-price_local[k].p;
                    j=k;
                }
            }
            double u=X[i][j]-price_local[j].p;
            if(u>0)
            {
                maxim=0;
                for(k=0;k<nr2;k++)
                {
                    if((X[i][k]-price_local[k].p)>=maxim&&k!=j)
                    {
                        maxim=X[i][k]-price_local[k].p;
                        sj=k;
                    }
                }
                double v=X[i][sj]-price_local[sj].p;

                double cp=price_local[j].p+u-v+elp_local;
                if(inprice_local[j].p<cp)
                {
                    inprice_local[j].p=cp;
                    inprice_local[j].buyer=sindex+i;
                }
            }
            else{
                //buyers_local.erase(i);
                //buyers.erase(i);
                //continue;
            }
        }
       for(int i=0;i<nr2;i++)
        {
            price_local[i].p=inprice_local[i].p;
            price_local[i].buyer=inprice_local[i].buyer;
        }

        pthread_mutex_lock(&mutex);
        for(int i=0;i<nr2;i++)
        {
            if(price[i].p<price_local[i].p)
            {
                price[i].p=price_local[i].p;
                price[i].buyer=price_local[i].buyer;
            }
        }

        pthread_mutex_unlock(&mutex);


         // Synchronization point
        int rc = pthread_barrier_wait(&barr);
        if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
        {
            printf("Could not wait on barrier\n");
            exit(-1);
        }


        //itrbuyer_local=buyers_local.begin();
        //for(;itrbuyer_local!=buyers_local.end();itrbuyer_local++)
        for(int i=start;i<endd;i++)
        {
            //int i=*itrbuyer_local;

            for(int j=0;j<nr2;j++)
            {
                pthread_mutex_lock(&mutex);
                int con=price_local[j].buyer==i&&price_local[j].p>=price[j].p&&!((price_local[j].p==price[j].p)&&price[j].p==0);
                int con1=price_local[j].buyer==i&&price_local[j].p<price[j].p;
  //              cout<<" i = "<<i<<" j =  "<<j<<"\n";
                pthread_mutex_unlock(&mutex);

                if(con)
                {
                    /*pthread_mutex_lock(&mutex);
                    cout<<"iteration number << "<<l<<"\n";
                    cout<<"buyer i = "<<i<<" prod j = "<<j<<" local price = "<<price_local[j].p<<" gl price = "<<price[j].p<<"\n\n";
                    pthread_mutex_unlock(&mutex);*/
                    for(int i1=start;i1<endd;i1++)
                    if((itr1=s_local.find(make_pair(i1,j)))!=s_local.end())
                    {

                        s_local.erase(itr1);
                        buyers_local.insert(i1);
                        if(s_local.size()==0)
                            break;

                   }
                    s_local.insert(make_pair(i,j));
                   if((itrbuyer_local=buyers_local.find(i))!=buyers_local.end())
                        buyers_local.erase(itrbuyer_local);
                   break;
                }
                else
                if(con1)
                {
                    //pthread_mutex_lock(&mutex);
                    //cout<<"iteration number << "<<l<<"\n";
                    //cout<<"check :: buyer i = "<<i<<" prod j = "<<j<<" local price = "<<price_local[j].p<<" gl price = "<<price[j].p<<"\n\n";
                  // pthread_mutex_unlock(&mutex);
                    if((itr1=s_local.find(make_pair(i,j)))!=s_local.end())
                    {


                        s_local.erase(itr1);
                        buyers_local.insert(i);
                        //price_local[j].p=price[j].p;
                         break;
                    }
                }
            }
        }
        rc = pthread_barrier_wait(&barr);
        if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
        {
            printf("Could not wait on barrier\n");
            exit(-1);
        }
        l++;
        pthread_mutex_lock(&mutex);
        for(int i=start;i<endd;i++)
        {
            if((itrbuyer_local=buyers_local.find(i))==buyers_local.end())
            {
                itrbuyer=buyers.find(i);
                if(itrbuyer!=buyers.end())
                    buyers.erase(itrbuyer);

            }
            else
            {
                buyers.insert(i);

            }
        }
        pthread_mutex_unlock(&mutex);
        rc = pthread_barrier_wait(&barr);
        if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
        {
            printf("Could not wait on barrier\n");
            exit(-1);
        }
    }
    pthread_mutex_lock(&mutex);
    itr1_local=s_local.begin();
    for(int i=0;i<s_local.size();i++,itr1_local++)
    {
        s.insert(*itr1_local);
    }
    pthread_mutex_unlock(&mutex);

}


void set_global_buyers()
{
    for(int i=0;i<nr1;i++)
        buyers.insert(i);
    itrbuyer=buyers.begin();
    price=new plocal[nr2];
    for(int i=0;i<nr2;i++)
    {
        price[i].p=0;
        price[i].buyer=0;
    }
}

int main(int argc,char ** argv)
{
    //similarity matrix is always of size of B (nB) * A(nA)
    // U matrix is always of size nB * nB
    // V matrix is always of size nA * nA
    // sigma matrix is always of size nB * nA

    //B matrix is nr1 * nc1 && A matrix is nr2 * nc2
    //nr1<=nr2
    int nc1,nc2;
    if(argc<3)
    {
        cout<<"enter ./a.out nr1 nr2\n";
        exit(1);
    }

    nr1=nc1=atoi(argv[1]);nr2=nc2=atoi(argv[2]);
    getsvd(nr1,nr2);

    cout<<"enter the value of s \n";
    cin>>s1;
    dummy_1d dum1;
    cout<<"\nenter the value of n : \n";
    cin>>n;
    cout<<"enter the value of alpha :\n";
    cin>>alp;
    fstream fin;fin.open("input.txt");
    fstream finsv;finsv.open("svdf.txt");

    fstream fina;fina.open("Amatrix.txt",ios::in|ios::out|ios::trunc);
    fstream finb;finb.open("Bmatrix.txt",ios::in|ios::out|ios::trunc);
    fstream finat;finat.open("atilde.txt",ios::in|ios::out|ios::trunc);
    fstream finbt;finbt.open("btilde.txt",ios::in|ios::out|ios::trunc);
    fstream finu;finu.open("umatrix.txt",ios::in|ios::out|ios::trunc);
    fstream finv;finv.open("vmatrix.txt",ios::in|ios::out|ios::trunc);
    fstream fins;fins.open("smatrix.txt",ios::in|ios::out|ios::trunc);
    finoutX.open("outmatrix.txt",ios::in|ios::out|ios::trunc);

    finX.open("Xmatrix.txt",ios::in|ios::out|ios::trunc);

    a=new MATRIX<double>(nr2,nc2,2);
    b=new MATRIX<double>(nr1,nc1,2);
    atilde=new MATRIX<double>(nr2,nc2,s1);
    btilde=new MATRIX<double>(nr1,nc1,s1);
    u=new MATRIX<double>(nr1,nr1,s1);
    v=new MATRIX<double>(nr2,nr2,s1);
    sigmam=new MATRIX<double>(nr1,nr2,s1);
    int s=min(nr1,nr2);
    Xi=new MATRIX<double>(s*nr1,nr2,s1);
    X=new MATRIX<double>(nr1,nr2,s1);

    a->initialise_matrix(fina);
    b->initialise_matrix(finb);
    atilde->initialise_matrix(finat);
    btilde->initialise_matrix(finbt);
    u->initialise_matrix(finu);
    v->initialise_matrix(finv);
    Xi->initialise_matrix(finX);
    X->initialise_matrix(finoutX);
    sigmam->initialise_matrix(fins);
    set_matrix(nr2,nc2,fin,fina,a);
    set_matrix(nr1,nc1,fin,finb,b);
    set_matrix(nr1,nr1,finsv,finu,u);
    set_matrix(nr2,nr2,finsv,finv,v);
    set_matrix(nr1,nr2,finsv,fins,sigmam);

    computtilde(nr2,nc2,a,atilde,fina,finat);
    computtilde(nr1,nc1,b,btilde,finb,finbt);
    vector<double> sigv;
    sigmam->fetch_diag(sigv,fins);
    //Initialization of iterates
    double val;dummy_2d temp;
    temp.resize(nr2);
    dummy_2d res;res.resize(nr2);
    fstream finw;finw.open("wmatrix.txt",ios::in|ios::out|ios::trunc);
    fstream finz;finz.open("zmatrix.txt",ios::in|ios::out|ios::trunc);
    w=new MATRIX<double>(s*(n+1),nr1,s1);
    z=new MATRIX<double>(s*(n+1),nr2,s1);
    w->initialise_matrix(finw);
    z->initialise_matrix(finz);

    for(int i=0;i<s;i++)
    {
        //generation of iterates
        vector<double> wi,zi;
        u->fetch_column(wi,i,0,nr1,finu);
        v->fetch_column(zi,i,0,nr2,finv);
        for(int j=0;j<nr1;j++)
        {
            wi[j]=wi[j]*sigv[i];
        }
        w->enter_row(wi,i*(n+1),0,nr1-1,finw);
        z->enter_row(zi,i*(n+1),0,nr2-1,finz);
        for(int k=1;k<=n;k++)
        {
            wi.clear();
            w->fetch_row(wi,i*(n+1)+k-1,0,nr1-1,finw);
            zi.clear();
            z->fetch_row(zi,i*(n+1)+k-1,0,nr2-1,finz);
            vector<double> res;
            w->mult_vect_mat(finbt,nr1,wi,res);

            w->enter_row(res,i*(n+1)+k,0,nr1-1,finw);
            res.clear();
            z->mult_vect_mat(finat,nr2,zi,res);
            z->enter_row(res,i*(n+1)+k,0,nr2-1,finz);
        }
    }
        finw.close();
        finz.close();
        //Finding X[i]

        factor=2;
        p=ceil(nr1/factor);q=ceil(nr2/factor);
        int pq=p*q;
        pthread_t threads[p*q+1];
        for(int i=0;i<p;i++)
        {
            for(int j=0;j<q;j++)
            {
                struct param *s1=new param;
                s1->r=i;
                s1->u=j;
                s1->s=s;
                pthread_create(&threads[i*p+j],NULL,part_sim,(void *)s1);
            }
        }
        for(int i=0;i<p;i++)
        {
            for(int j=0;j<q;j++)

            {
                pthread_join(threads[i*p+j],NULL);
            }
        }



        cout<<"simi matrix\n\n";
        dummy_2d Xlastres;
        X->fetch_matrix(Xlastres,0,nr1-1,0,nr2-1,finoutX);
		for(int i=0;i<Xlastres.size();i++)
			Xlastres[i][i]=0;
        X->enter_matrix(Xlastres,0,nr1-1,0,nr2-1,finoutX);
       // sleep(1);
		printmatrix(Xlastres,nr1,nr2);
		set_global_buyers();

        int NP=ceil(nr1/factor);
        // Barrier initialization
        if(pthread_barrier_init(&barr, NULL, NP))
        {
            printf("Could not create a barrier\n");
            return -1;
        }
        for(int i=0;i<NP;i++)
        {
            int id=i+1;
            pthread_create(&threads[i],NULL,match,(int *)&id);
            sleep(1);
        }
        for(int i=0;i<NP;i++)
        {
            pthread_join(threads[i],NULL);
        }
        print();
        fina.close();
        finb.close();
        finat.close();
        finbt.close();
        finu.close();
        finv.close();
        fins.close();
        finsv.close();
        finX.close();
        finoutX.close();
        fin.close();
        free(a);free(b);free(atilde);free(btilde);free(u);free(v);free(sigmam);
        free(w);free(z);free(Xi);free(X);
}

