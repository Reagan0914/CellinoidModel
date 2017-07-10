#include "mex.h"
#include <math.h>
#include <stdio.h>

int eejcb(double a[],double v[]);
int eejcb_luxp(double a[],double v[]);
void mexFunction(int hlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *V;
    double *A;
    
    double MatrixA[9];  //for input matrix C and store eigen values in diagnal
    double MatrixV[9];  //for eigen values
    int i,j;
    
    A = mxGetPr(prhs[0]);
    for(i=0;i<9;i++)
    {
        MatrixA[i] = A[i];
    }
    
    if(eejcb_luxp(MatrixA, MatrixV)!=1)
        mexErrMsgTxt("Error in Eig.c Function!\n");
    
    //output
    plhs[0] = mxCreateDoubleMatrix(3,3, mxREAL);
    V = mxGetPr(plhs[0]);
    for(i=0;i<9;i++)
        V[i] = MatrixV[i];
    
    plhs[1] = mxCreateDoubleMatrix(3,3,mxREAL);
    A = mxGetPr(plhs[1]);
    for(i=0;i<9;i++)
    {
        if(i%4 == 0)
            A[i] = MatrixA[i];
        else
            A[i] = MatrixA[i];//0.0
    }
}
/*  Calculate Eigen Pair by Jacaobi Method (Modified by LUXP)*/
int eejcb_luxp(double a[],double Q[])
{
	int i,j,p,q,u,w,t,s,l;
	const int n=3;
	double eps=0.00000000001;
	int jt=100;
    double fm,cn,sn,omega,x,y,d;
    int Order[3];
    double v[9];
    l=1;
	
    for (i=0; i<=n-1; i++)
	{
		v[i*n+i]=1.0;
        for (j=0; j<=n-1; j++)
		{
			if (i!=j)
			{
				v[i*n+j]=0.0;
			}
		}
	}
    while (1==1)
	{
		fm=0.0;
        for (i=0; i<=n-1; i++)
		{
			for (j=0; j<=n-1; j++)
			{
				d=fabs(a[i*n+j]);
				if ((i!=j)&&(d>fm))
				{
					fm=d;
					p=i;
					q=j;
				}
			}
		}
        if (fm<eps)
		{
                
            //modified by LUXP, adjusted the eigen vector sequence to Right-Hand-System
            //Order: Get the order of 3 eigen values from smallest to largest
            if(a[0] <= a[4])
            {
                if (a[4] <= a[8])
                {
                    Order[0] = 0;
                    Order[1] = 1;
                    Order[2] = 2;
                }else if(a[8] >= a[0])
                {
                    Order[0] = 0;
                    Order[1] = 2;
                    Order[2] = 1;
                }else{
                    Order[0] = 2;
                    Order[1] = 0;
                    Order[2] = 1;
                }

            }else{
                if (a[4] >= a[8])
                {
                    Order[0] = 2;
                    Order[1] = 1;
                    Order[2] = 0;
                }else if(a[8] >= a[0])
                {
                    Order[0] = 1;
                    Order[1] = 0;
                    Order[2] = 2;
                }else{
                    Order[0] = 1;
                    Order[1] = 2;
                    Order[2] = 0;
                }        
            }
            //adjust eigen vectors sequence according to the eigen values order
            Q[0]= v[Order[0]];    Q[3]= v[Order[1]];     Q[6]=v[Order[2]];
            Q[1]= v[Order[0]+3];  Q[4]= v[Order[1]+3];   Q[7]=v[Order[2]+3];
            Q[2]= v[Order[0]+6];  Q[5]= v[Order[1]+6];   Q[8]=v[Order[2]+6];
            //store 3 eigenvalues to a[1,2,3] from small to large,original 0.0, then go back!
            a[1]= a[4*Order[0]];    a[2]= a[4*Order[1]];    a[3]= a[4*Order[2]]; 
            a[0]= a[1]; a[4]= a[2];  a[8]= a[3];
            a[1]=0.0;   a[2]= 0.0;  a[3]= 0.0;
            //adjust to Right Hand System, adjust x,z axis, then cross for y axis
            if (Q[0] < 0)
            {
                Q[0] *= -1;
                Q[1] *= -1;
                Q[2] *= -1;
            }
            if (Q[8] < 0)
            {
                Q[6] *= -1;
                Q[7] *= -1;
                Q[8] *= -1;
            }
            Q[3]=Q[7]*Q[2]-Q[8]*Q[1];
            Q[4]=Q[8]*Q[0]-Q[6]*Q[2];
            Q[5]=Q[6]*Q[1]-Q[7]*Q[0];
            
            //end modification
			return(1);
		}
        if (l>jt)
		{
			return(-1);
		}
        l=l+1;
        u=p*n+q;
		w=p*n+p;
		t=q*n+p;
		s=q*n+q;
        x=-a[u];
		y=(a[s]-a[w])/2.0;
        omega=x/sqrt(x*x+y*y);
        if (y<0.0)
		{
			omega=-omega;
		}
        sn=1.0+sqrt(1.0-omega*omega);
        sn=omega/sqrt(2.0*sn);
        cn=sqrt(1.0-sn*sn);
        fm=a[w];
        a[w]=fm*cn*cn+a[s]*sn*sn+a[u]*omega;
        a[s]=fm*sn*sn+a[s]*cn*cn-a[u]*omega;
        a[u]=0.0;
		a[t]=0.0;
        for (j=0; j<=n-1; j++)
		{
			if ((j!=p)&&(j!=q))
			{
				u=p*n+j;
				w=q*n+j;
				fm=a[u];
				a[u]=fm*cn+a[w]*sn;
				a[w]=-fm*sn+a[w]*cn;
			}
		}
        for (i=0; i<=n-1; i++)
		{
			if ((i!=p)&&(i!=q))
            {
				u=i*n+p;
				w=i*n+q;
				fm=a[u];
				a[u]=fm*cn+a[w]*sn;
				a[w]=-fm*sn+a[w]*cn;
            }
		}
        for (i=0; i<=n-1; i++)
		{
			u=i*n+p;
			w=i*n+q;
            fm=v[u];
            v[u]=fm*cn+v[w]*sn;
            v[w]=-fm*sn+v[w]*cn;
		}
	}
    
	return(1);
}

/*  Calculate Eigen Pair by Jacaobi Method */
int eejcb(double a[],double v[])
{
	int i,j,p,q,u,w,t,s,l;
	int n=3;
	double eps=0.00000000001;
	int jt=100;
    double fm,cn,sn,omega,x,y,d;
    l=1;
	
    for (i=0; i<=n-1; i++)
	{
		v[i*n+i]=1.0;
        for (j=0; j<=n-1; j++)
		{
			if (i!=j)
			{
				v[i*n+j]=0.0;
			}
		}
	}
    while (1==1)
	{
		fm=0.0;
        for (i=0; i<=n-1; i++)
		{
			for (j=0; j<=n-1; j++)
			{
				d=fabs(a[i*n+j]);
				if ((i!=j)&&(d>fm))
				{
					fm=d;
					p=i;
					q=j;
				}
			}
		}
        if (fm<eps)
		{
			return(1);
		}
        if (l>jt)
		{
			return(-1);
		}
        l=l+1;
        u=p*n+q;
		w=p*n+p;
		t=q*n+p;
		s=q*n+q;
        x=-a[u];
		y=(a[s]-a[w])/2.0;
        omega=x/sqrt(x*x+y*y);
        if (y<0.0)
		{
			omega=-omega;
		}
        sn=1.0+sqrt(1.0-omega*omega);
        sn=omega/sqrt(2.0*sn);
        cn=sqrt(1.0-sn*sn);
        fm=a[w];
        a[w]=fm*cn*cn+a[s]*sn*sn+a[u]*omega;
        a[s]=fm*sn*sn+a[s]*cn*cn-a[u]*omega;
        a[u]=0.0;
		a[t]=0.0;
        for (j=0; j<=n-1; j++)
		{
			if ((j!=p)&&(j!=q))
			{
				u=p*n+j;
				w=q*n+j;
				fm=a[u];
				a[u]=fm*cn+a[w]*sn;
				a[w]=-fm*sn+a[w]*cn;
			}
		}
        for (i=0; i<=n-1; i++)
		{
			if ((i!=p)&&(i!=q))
            {
				u=i*n+p;
				w=i*n+q;
				fm=a[u];
				a[u]=fm*cn+a[w]*sn;
				a[w]=-fm*sn+a[w]*cn;
            }
		}
        for (i=0; i<=n-1; i++)
		{
			u=i*n+p;
			w=i*n+q;
            fm=v[u];
            v[u]=fm*cn+v[w]*sn;
            v[w]=-fm*sn+v[w]*cn;
		}
	}    
	return(1);
}