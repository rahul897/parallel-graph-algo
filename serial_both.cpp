#include <iostream>
#include <math.h>
#include <fstream>
#include <set>
#include <vector>
#include <utility>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

set <pair<int,int> > s;
set <pair<int,int> >::iterator itr,itr1;
set <int> buyers;
set <int>::iterator itrbuyer;
#define MAT_SIZE 100
#define dummy_1d vector <double>
#define dummy_2d vector <vector<double > >

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

//generate_matrix reads from the file the matrices into a Matrix object
void generate_matrix(dummy_2d &A,int a,int c,int b,ifstream &fin)
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
void  computtilde(dummy_2d &Atilde,int nr1,int nc1,dummy_2d A)
{
    int rv[MAT_SIZE];
    for(int i=0;i<nr1;i++)
    {
        long long int sum=0;
        for(int j=0;j<nc1;j++)
            sum+=A[i][j];
        rv[i]=sum;
    }
    for(int i=0;i<nr1;i++)
        cout<<rv[i]<<" ";
    cout<<"\n\n";

    for(int i=0;i<nr1;i++)
    {
        for(int j=0;j<nc1;j++)
        {
            if(rv[j]!=0)
            Atilde[i].push_back(A[j][i]/rv[j]);
            else
            Atilde[i].push_back(0);
        }
    }
}
//generate_SVD generates the SVD of a given matrix
//It generates the left singular vector,right singular vector and singular values matrix
void generateSVD(dummy_2d &u,dummy_2d &v,dummy_2d &sigma,int nr1,int nr2)
{
    ifstream finsvd;
    finsvd.open("svdf.txt");
    generate_matrix(u,nr1,nr1,0,finsvd);
    generate_matrix(v,nr2,nr2,0,finsvd);
    generate_matrix(sigma,nr1,nr2,0,finsvd);
}

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
void compute_transpose(dummy_2d a,int r,int c)
{
     for(int i=0;i<r;i++)
        for(int j=0;j<c;j++)
            a[i][j]=a[j][i];
}

void printmatrix(dummy_2d a,int r,int c)
{
    for(int i=0;i<r;i++)
    {
        for(int j=0;j<c;j++)
            cout<<a[i][j]<<" ";
        cout<<"\n";
    }
}

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
int main(int argc,char ** argv)
{
    int nr1,nc1,nr2,nc2,n;double alp;
    /*if(argc < 3)
    {
        cout << "usage :  ./a.out num_vertices1 num_vertices2 \n";
        exit(1);
    }
    nr1=atoi(argv[1]);
    nc1=nr1;
    nr2=atoi(argv[2]);
    nc2=nr2;*/
    nr1=nr2=nc1=nc2=8;
    dummy_1d dum1;
    cout<<"\nenter the value of n : \n";
    cin>>n;
    cout<<"enter the value of alpha :\n";
    cin>>alp;

    ifstream fin;fin.open("input.txt");
    dummy_2d B;
    for(int i=0;i<nr1;i++)
        B.push_back(dum1);
    generate_matrix(B,nr1,nc1,2,fin);
    dummy_2d A;
    for(int i=0;i<nr2;i++)
        A.push_back(dum1);
    generate_matrix(A,nr2,nc2,2,fin);
    cout<<"\n----------Adjacency Matrix B--------------\n ";
    printmatrix(B,nr1,nc1);
    cout<<"\n";
    cout<<"\n----------Adjacency Matrix A--------------\n ";
    printmatrix(A,nr2,nc2);
    cout<<"\n";
    dummy_2d Btilde,Atilde;
    for(int i=0;i<nr1;i++)
        Btilde.push_back(dum1);
    for(int i=0;i<nr2;i++)
        Atilde.push_back(dum1);


    computtilde(Btilde,nr1,nc1,B);
    computtilde(Atilde,nr2,nc2,A);
    //compute_transpose(Atilde,nr2,nc2);
    cout<<"\n----------Tilde Matrix B--------------\n ";
    printmatrix(Btilde,nr1,nc1);
    cout<<"\n----------Tilde Matrix A--------------\n ";
    printmatrix(Atilde,nr1,nc1);

    dummy_2d u,v,sigma;
    u.resize(MAT_SIZE);v.resize(MAT_SIZE);sigma.resize(MAT_SIZE);
    generateSVD(u,v,sigma,nr1,nr2);
    cout<<"\n----------SVD DECOMPOSITION--------------\n ";

   // u=u*-1;v=v*-1;
    //to verify if SVD has taken place properly
    //cout<<u * sigma * v.transpose();

    cout<<"\n\nU Matrix-----\n";
    printmatrix(u,nr1,nc1);
    cout<<"\n\nV Matrix -----\n";
    printmatrix(v,nr2,nc2);
    cout<<"\n\nSigma Matrix -----------\n";
    printmatrix(sigma,nr1,nr2);


    int s=min(nr1,nr2);
    vector<double> w[s],z[s];
    //Initialization of iterates
    for(int i=0;i<s;i++)
    {
        vector<double> vect;
        for(int j=0;j<nr1;j++)
            vect.push_back(sigma[i][i]*u[j][i]);
        w[i]=vect;
        vect.clear();
        for(int j=0;j<nr2;j++)
            vect.push_back(v[j][i]);
        z[i]=vect;
    }

    dummy_2d X[s];
    for(int i=0;i<s;i++)
    {
        X[i].resize(nr1);
    }
    dummy_2d Xlast;Xlast.resize(nr1);
    for(int i=0;i<nr1;i++)
        for(int j=0;j<nr2;j++)
            Xlast[i].push_back(0);
    dummy_2d Xlastres;Xlastres.resize(nr1);
        double val;dummy_2d temp;
            temp.resize(nr1);
    dummy_2d res;res.resize(nr1);


    for(int i=0;i<s;i++)
    {
        //generation of iterates
        for(int l=0;l<nr1;l++)
            for(int m=0;m<nr2;m++)
                X[i][l].push_back(0);
        vector<double> itrw[n+1],itrz[n+1];
        itrw[0]=w[i];
        itrz[0]=z[i];
        for(int k=1;k<=n;k++)
        {
           mult_vect_mat(Btilde , itrw[k-1],itrw[k]);
            mult_vect_mat(Atilde , itrz[k-1],itrz[k]);
        }
        //Finding X[i]
        for(int k=0;k<=n-1;k++)
        {
            val=pow(alp,k);
            mult_2_mat(itrw[k],itrz[k],temp,val);
            mat_add(temp,X[i],res);
            X[i]=res;
        }
        val=pow(alp,n);
        mult_2_mat(itrw[n],itrz[n],temp,val);
        for(int p=0;p<X[i].size();p++)
            for(int q=0;q<X[i][0].size();q++)
                X[i][p][q]=(1-alp)*X[i][p][q];
        mat_add(temp,X[i],res);
        X[i]=res;
    }
    //Taking the sum of all X[i]'s
   for(int i=0;i<s;i++)
    {
        mat_add(Xlast,X[i],Xlastres);
        Xlast=Xlastres;

    }
    cout<<"simi matrix\n\n";
    printmatrix(Xlastres,nr1,nr2);
    matching(Xlast,nr1,nr2);
    print();
}
