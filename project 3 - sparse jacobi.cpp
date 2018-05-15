#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;


// diag(i)=A(i,i)
// lx(i)=A(i-1,i)
// ly(i)=A(i-ngrid,i)
// ux(i)=A(i,i+1)
// uy(i)=A(i,i+ngrid)
// n=ngrid*ngrid

void BminusAX(int n, int ngrid, double* ly, double* lx, double* diag, 
			  double* ux, double* uy, double* x, double* b, double* r) 
{
	for (int i=0;i<n;i++) 
	{
		r[i]=b[i]-diag[i]*x[i];
		if (i>0)
			r[i]-=lx[i]*x[i-1];
		if (i+1<n)
			r[i]-=ux[i]*x[i+1];
		if (i>ngrid-1)
			r[i]-=ly[i]*x[i-ngrid];
		if (i+ngrid<n)
			r[i]-=uy[i]*x[i+ngrid];
	}
}

void Atimesd(int n, int ngrid, double* ly, double* lx, double* diag, 
			 double* ux, double* uy, double* d, double* Ad) 
{
	for (int i=0;i<n;i++) 
	{
		Ad[i]=diag[i]*d[i];
		if (i>0)
			Ad[i]+=lx[i]*d[i-1];
		if (i+1<n)
			Ad[i]+=ux[i]*d[i+1];
		if (i>ngrid-1)
			Ad[i]+=ly[i]*d[i-ngrid];
		if (i+ngrid<n)
			Ad[i]+=uy[i]*d[i+ngrid];
	} 
}

void PolyPre(int p, int n, int ngrid, double* ly, double* lx, double* diag, 
			  double* ux, double* uy, double* r, double* z)
{
	double* res=new double[n];

	for (int i=0;i<n;i++)
		z[i] = 0;		//setting z[0] = 0 vector
	
	if (p==0)
	{
		for(int i=0; i<n; i++)
			z[i] = r[i];
	}

	if(p==1)
	{
		for(int i=0; i<n; i++)
			z[i] = r[i] / diag[i];
	}

	if (p>1)
	{
		for (int j=0;j<=p;j++)
		{
			for (int i=0;i<n;i++)
			{
				BminusAX(n, ngrid, ly, lx, diag, ux, uy, z, r, res); //res = r-Az
				z[i] += res[i]/diag[i];
			}
		}
	}
}


int main() 
{
	int n,ngrid,N,p;
	double h,eps, alpha, beta;

	eps=1.0e-10;
	ngrid=2;
	h=1.0/ngrid;
	n=ngrid*ngrid;
	srand(0);

	N = 1000;
	
	cout << "Enter value for p: ";
	cin >> p;
 
	// A_{i,i}
	double* diag=new double[n];
   // A_{i,i-1}
	double* lx=new double[n];
   // A_{i,i+1}
	double* ux=new double[n];
   // A_{i,i-ngrid}
	double* ly=new double[n];
   // A_{i,i+ngrid}
	double* uy=new double[n];

	double* b=new double[n];
	double* r=new double[n];
	double* x=new double[n];
	double* Ad=new double[n];
	double* d=new double[n];
	double* s=new double[n];
	
	for (int i=0;i<n;i++) 
	{
		int jgrid = i/ngrid;
		int igrid = i%ngrid;
		lx[i]=0.0;
		ly[i]=0.0;
		ux[i]=0.0;
		uy[i]=0.0;
		diag[i]=4.0;
		if (igrid>0)
			lx[i]=-1.0;
		if (igrid<ngrid-1)
			ux[i]=-1.0;
		if (jgrid>0)
			ly[i]=-1.0;
		if (jgrid<ngrid-1)
			uy[i]=-1.0;	
		b[i]=static_cast<double> (rand() % 100);  
		b[i]*=(h*h);
	}
 
	//start of painless algorithm
	int k=0;
	for(int i=0; i<n; i++)
		x[i] = 0;		//setting x_0 = {0,...,0}

	BminusAX(n, ngrid, ly, lx, diag, ux, uy, x, b, r); //solving r=b-Ax

	PolyPre(p,n,ngrid,ly,lx,diag,ux,uy,r,d);
	
	
	double DeltaOld;
	double DeltaNew = 0;
	for (int i=0; i<n; i++)
		DeltaNew += r[i] * d[i];
	DeltaOld = DeltaNew;

	while ((k<N) && (DeltaNew > (eps*eps*DeltaOld)))
		{
			Atimesd(n, ngrid, ly, lx, diag, ux, uy, d, Ad);	//q=Ad

			alpha = 0;
			for(int i=0; i<n; i++)
				alpha += d[i] * Ad[i];
			alpha = DeltaNew / alpha;
			
			for(int i=0; i<n; i++)
				x[i] += alpha * d[i];	
			
			BminusAX(n, ngrid, ly, lx, diag, ux, uy, x, b, r); //solving r=b-Ax
			
			PolyPre(p,n,ngrid,ly,lx,diag,ux,uy,r,s);


			DeltaOld = DeltaNew;
			DeltaNew = 0;
			for(int i=0; i<n; i++)
				DeltaNew += r[i] * s[i];
			
			beta = DeltaNew / DeltaOld;
			
			for(int i=0; i<n; i++)
				d[i] = s[i] + (beta * d[i]);

				k++;
		}

	cout << "great success!" << endl;
	cout << " x = ( ";
	for(int i=0; i<n; i++)
		cout << x[i] << " ";
	cout << ")" << endl;
	cout << "in " << k << " iterations" << endl;
}

