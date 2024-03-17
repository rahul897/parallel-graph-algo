#include <iostream>
#include <octave/oct.h>
#include <ov.h>
#include <fstream>
#include<set>
#include<utility>
#include <octave/builtin-defun-decls.h>
using namespace std;
set <pair<int,int> > s;
set <pair<int,int> >::iterator itr,itr1;
set <int> buyers;
set <int>::iterator itrbuyer;
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
Matrix generate_matrix(int a,int c,int b,ifstream &fin)
{
    Matrix A=Matrix(a,c);
    for(int i=0;i<a;i++)
            for(int j=0;j<c;j++)
                fin>>A(i,j);
    return A;
}
//computetilde computes the tilde of a given matrix and size
Matrix computtilde(int nr1,int nc1,Matrix A)
{
    RowVector rv=RowVector(nr1);
    for(int i=0;i<nr1;i++)
    {
        Matrix c=A.row(i);
        c=c.sum();
        rv(i)=c(0);
    }
    cout<<rv<<"\n\n";
    Matrix Atilde=Matrix(nr1,nc1);
    for(int i=0;i<nr1;i++)
    {
        for(int j=0;j<nc1;j++)
        {
            if(rv(j)!=0)
            Atilde(i,j)=A(j,i)/rv(j);
            else
            Atilde(i,j)=0;
        }
    }
    return Atilde;
}
//generate_SVD generates the SVD of a given matrix
//It generates the left singular vector,right singular vector and singular values matrix
void generateSVD(Matrix &u,Matrix &v,DiagMatrix &sigma,int nr1,int nr2)
{
    octave_value_list in,out;
    in(0)=nr1;
    in(1)=nr2;
    out=Frand(in,2);
    Matrix H=out(0).matrix_value();
    cout<<"\n----------Initial  Matrix H--------------\n ";
    cout<<"\n\n"<<H;
    SVD s=SVD(H);
    u=s.left_singular_matrix();
    v=s.right_singular_matrix();
    sigma=s.singular_values();
}

void matching(Matrix X,int nr1,int nr2)
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
            if(X(i,k)-price[k]>=maxim)
            {
                maxim=X(i,k)-price[k];
                j=k;
            }
        }

        maxim=0;
        for(k=0;k<nr2;k++)
        {
            if(X(i,k)-price[k]>=maxim&&k!=j)
            {
                maxim=X(i,k)-price[k];
                sj=k;
            }
        }

        double u=X(i,j)-price[j];
        double v=X(i,sj)-price[sj];

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
    //cout<<s.size()<<"\n";
    //cout<<buyers.size()<<"\n";
    itr1=s.begin();
    for(int i=0;i<s.size();i++)
    {
        cout<<"vertex "<<(*itr1).first+1<<" of first graph  is most similiar to  vertex "<<(*itr1).second+1<<" of second graph\n";
        itr1++;
    }
    cout<<"\n";

}

int main(int argc,char ** argv)
{
    int nr1,nc1,nr2,nc2,n;double alp;
    if(argc < 3)
    {
        cout << "usage :  ./a.out num_vertices1 num_vertices2 \n";
        exit(1);
    }
    nr1=atoi(argv[1]);
    nc1=nr1;
    nr2=atoi(argv[2]);
    nc2=nr2;

    cout<<"\nenter the value of n : \n";
    cin>>n;
    cout<<"enter the value of alpha :\n";
    cin>>alp;

    ifstream fin;fin.open("input.txt");
    Matrix B=generate_matrix(nr1,nc1,2,fin);
    Matrix A=generate_matrix(nr2,nc2,2,fin);
    cout<<"\n----------Adjacency Matrix B--------------\n ";
    cout<<B<<"\n\n";
    cout<<"\n----------Adjacency Matrix A--------------\n ";
    cout<<A<<"\n\n";
    Matrix Btilde=computtilde(nr1,nc1,B);
    Matrix Atilde=computtilde(nr2,nc2,A);
    Matrix Atildetr=Atilde.transpose();
    cout<<"\n----------Tilde Matrix B--------------\n ";
    cout<<Btilde<<"\n\n";
    cout<<"\n----------Transpose Tilde Matrix A--------------\n ";
    cout<<Atildetr<<"\n\n";
    Matrix u,v;DiagMatrix sigma;

    generateSVD(u,v,sigma,nr1,nr2);
    cout<<"\n----------SVD DECOMPOSITION--------------\n ";

    u=u*-1;v=v*-1;
    //to verify if SVD has taken place properly
    cout<<u * sigma * v.transpose();
    cout<<"\n\nU Matrix-----\n"<<u<<"\n\nV Matrix -----\n"<<v<<"\n\nSigma Matrix -----------\n"<<sigma<<"\n\n";

    int s=min(nr1,nr2);
    ColumnVector w[s],z[s];
    //Initialization of iterates
    for(int i=0;i<s;i++)
    {
        w[i]=sigma(i,i)*u.column(i);
        z[i]=v.column(i);
        //cout<<w[i]<<" ";
    }

    Matrix X[s];
    Matrix Xlast=Matrix(nr1,nr2,0);
    for(int i=0;i<s;i++)
        X[i]=Matrix(nr1,nr2,0);

    for(int i=0;i<s;i++)
    {
        //generation of iterates
        ColumnVector itrw[n+1],itrz[n+1];
        itrw[0]=w[i];
        itrz[0]=z[i];
        for(int k=1;k<=n;k++)
        {
            itrw[k]=Btilde * itrw[k-1];
            itrz[k]=Atilde * itrz[k-1];
          //  cout<<itrw[k]<<"hah\n";

        }

    //Finding X[i]
        for(int k=0;k<=n-1;k++)
        {
            X[i]=X[i]+pow(alp,k)*itrw[k]*itrz[k].transpose();
        }
        X[i]=(1-alp)*X[i]+pow(alp,n)*itrw[n]*itrz[n].transpose();
    }
    //cout<<Xlast<<"\n\n";
    //Taking the sum of all X[i]'s
    for(int i=0;i<s;i++)
    {
        //cout<<"\n"<<X[i]<<"\n\n";
        Xlast=Xlast+X[i];
    }
    cout<<"simi matrix\n\n"<<Xlast<<"\n";
    ofstream fout;fout.open("outmat.txt");
    fout<<Xlast;
    matching(Xlast,nr1,nr2);
    print();
}
