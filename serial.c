/*  Author: Kareem Abdelshafy. Date: 04/14/2016
		 High Performance computing Final Project
serial version of navier stoke solver for a fluid flow inside a channel   */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <time.h>

int main ( int argc, char **argv)
{

	// variables declaration and initalization
	double rho = 1; double Re = 10; double Nx = 128; double Ny = 64; 
	double dt = 0.0005; double tend = 0.5; double L = 2; double h = 1;
	double dx = L/Nx;  double dy = h/Ny; double B=dx/dy;
	double w=1; int flag1=1 ; int flag=1; double R=0; double R1=0;
	double R2=0; double R3=0; int count=0; int count1=0;

	// array declaration  
	int ry= Ny+2;  int cx = Nx+2;

	double **u = (double **)malloc(ry * sizeof(double *));  // u
	for (int j=0; j<ry; j++)  
		{u[j] = (double *)malloc(cx * sizeof(double));}

	double **v = (double **)malloc(ry * sizeof(double *));  // v
	for (int j=0; j<ry; j++)  
		{v[j] = (double *)malloc(cx * sizeof(double));}

	double **P = (double **)malloc(ry * sizeof(double *));  // P
	for (int j=0; j<ry; j++)  
		{P[j] = (double *)malloc(cx * sizeof(double));}

	double **uo = (double **)malloc(ry * sizeof(double *));  // uo
	for (int j=0; j<ry; j++)  
		{uo[j] = (double *)malloc(cx * sizeof(double));}

	double **vo = (double **)malloc(ry * sizeof(double *));  // vo
	for (int j=0; j<ry; j++)  
		{vo[j] = (double *)malloc(cx * sizeof(double));}

	double **C = (double **)malloc(ry * sizeof(double *));  // C
	for (int j=0; j<ry; j++)  
		{C[j] = (double *)malloc(cx * sizeof(double));}

	double **Pa = (double **)malloc(ry * sizeof(double *));  // Pa
	for (int j=0; j<ry; j++)  
		{Pa[j] = (double *)malloc(cx * sizeof(double));}

	double **Pb = (double **)malloc(ry * sizeof(double *));  // Pb
	for (int j=0; j<ry; j++)  
		{Pb[j] = (double *)malloc(cx * sizeof(double));}

	double **Po = (double **)malloc(ry * sizeof(double *));  // Po
	for (int j=0; j<ry; j++)  
		{Po[j] = (double *)malloc(cx * sizeof(double));}


	clock_t begin, end;  //  Start time counting
	double time_spent;
	begin = clock();


	// initial condition

	for (int j=0 ; j < ry ;j++)  {	//u
	    for (int i=0 ; i < cx ;i++) {
		u[j][i]=0 ; } }

	for (int j=0 ; j < ry ;j++) {	//u
		u[j][0]=1 ; } 

	for (int j=0 ; j < ry ;j++)  {	//v
	    for (int i=0 ; i < cx ;i++) {
		v[j][i]=0 ; } }

	for (int j=0 ; j < ry ;j++)  {	//P
	    for (int i=0 ; i < cx ;i++) {
		P[j][i]=0 ; } }

	for (int j=0 ; j < ry ;j++)  {	//Pa
	    for (int i=0 ; i < cx ;i++) {
		Pa[j][i]=0 ; } }

	for (int j=0 ; j < ry ;j++)  {	//Pb
	    for (int i=0 ; i < cx ;i++) {
		Pb[j][i]=0 ; } }

	for (int j=0 ; j < ry ;j++)  {	//Po
	    for (int i=0 ; i < cx ;i++) {
		Po[j][i]=0 ; } }

	for (int j=0 ; j < ry ;j++)  {	//uo
	    for (int i=0 ; i < cx ;i++) {
		uo[j][i]=u[j][i] ; } }

	for (int j=0 ; j < ry ;j++)  {	//C
	    for (int i=0 ; i < cx ;i++) {
		C[j][i]=0 ; } }

 	
//////////////////////////// start of time undependant solver //////////////////

 for (double t=0 ; t < tend ; t+=dt ) {

printf("t = %f count1= %d  \n", t , count1);

flag1=1;  flag=1; 

    while (flag1 == 1) {

	count1++ ;

	for (int j=1 ; j < ry-1 ;j++)  {	//uo
	    for (int i=1 ; i < cx-1 ;i++) {

		uo[j][i]=u[j][i] - dt/dx/4* ( pow(u[j][i]+u[j][i+1],2) - pow(u[j][i]+u[j][i-1],2) )

                -dt/dy/4* ( (u[j][i]+u[j+1][i]) * (v[j][i]+v[j][i+1]) - (u[j][i]+u[j-1][i]) * (v[j-1][i]+v[j-1][i+1]) )

                -( P[j][i+1]-P[j][i] )*dt/dx +

		 dt/Re*( (u[j][i+1]-2*u[j][i]+u[j][i-1])/(dx*dx) + (u[j+1][i]-2*u[j][i]+u[j-1][i])/(dy*dy) ); } }






	for (int j=0 ; j < ry ;j++) {	//uo left boundary
		uo[j][0]=1 ; } 

	for (int i=0 ; i < cx ;i++) {	//uo lower boundary
		uo[0][i]=-1*uo[1][i] ; } 

	for (int i=0 ; i < cx ;i++) {	//uo  upper boundary
		uo[ry-1][i]=-1*uo[ry-2][i] ; } 

	for (int j=0 ; j < ry ;j++) {	//uo  right boundary
		uo[j][cx-1]=uo[j][cx-2] ; } 






	for (int j=1 ; j < ry-1 ;j++)  {	//vo
	    for (int i=1 ; i < cx-1 ;i++) {

            vo[j][i]=v[j][i] -dt/dy/4 * ( pow(v[j][i]+v[j+1][i],2) - pow(v[j][i]+v[j-1][i],2) )

                -dt/dx/4* ( (u[j][i]+u[j+1][i])*(v[j][i]+v[j][i+1]) - (u[j][i-1]+u[j+1][i-1])*(v[j][i]+v[j][i-1]) )

                -(P[j+1][i]-P[j][i])*dt/dy

		+dt/Re*( (v[j][i+1]-2*v[j][i]+v[j][i-1])/(dx*dx) +(v[j+1][i]-2*v[j][i]+v[j-1][i])/(dy*dy) ) ; } }

	for (int j=0 ; j < ry ;j++) {	//vo left boundary
		vo[j][0]=-1*vo[j][1] ; } 

	for (int i=0 ; i < cx ;i++) {	//vo lower boundary
		vo[0][i]=0 ; } 

	for (int i=0 ; i < cx ;i++) {	//vo  upper boundary  tricky !!
		vo[ry-2][i]=0 ; } 

	for (int j=0 ; j < ry ;j++) {	//vo  right boundary
		vo[j][cx-1]=vo[j][cx-2] ; } 






	for (int j=1 ; j < ry-1 ;j++)  {	//C
	    for (int i=1 ; i < cx-1 ;i++) {
		C[j][i]= (uo[j][i]-uo[j][i-1])*dx/dt+(vo[j][i]-vo[j-1][i])*B*B*dy/dt ; } }


//////////////////////// LSOR  Poission Equation ////////////////////

	while (flag == 1) {
	
	count++ ;

	for (int j=0 ; j < ry ;j++)  {	// Pt = Pa
	    for (int i=0 ; i < cx ;i++) {
		Pb[j][i]=Pa[j][i];  } }



	for (int j=1 ; j < ry-1 ;j++)  {	//C
	    for (int i=1 ; i < cx-1 ;i++) {
                Pa[j][i]=0.25*(Pb[j][i+1]+Pb[j][i-1]+B*B*(Pb[j+1][i]+Pb[j-1][i])-C[j][i]);  } }

	 R=0;
	for (int j=0 ; j < ry ;j++)  {	// Tolerance check
	    for (int i=0 ; i < cx ;i++) {
		double Rtemp= fabs((Pa[j][i]-Pb[j][i]) / (Pb[j][i]+0.0000001));

		R=std::max(R,Rtemp);	} }

	if ( R > 0.001) { flag = 1;}   else  {flag = 2;}

	}

	flag = 1;
///////////////////////  end of LSOR   //////////////////////////////

	for (int j=0 ; j < ry ;j++)  {	// Pt = P
	    for (int i=0 ; i < cx ;i++) {
		Po[j][i]=P[j][i];  } }

	for (int j=1 ; j < ry-1 ;j++)  {  // P = P + Pa
	    for (int i=1 ; i < cx-1 ;i++) {
		P[j][i]=P[j][i] + Pa[j][i];  } }

	R1=0;
	for (int j=0 ; j < ry ;j++)  {	// Tolerance check
	    for (int i=0 ; i < cx ;i++) {
		double R1temp= fabs((P[j][i]-Po[j][i]) / (Po[j][i]+0.0000001));
		R1=std::max(R1,R1temp);	} }

	if (R1 > 0.00001) { flag1 = 1;}   else  {flag1 = 2;}

//printf("flag1 %d \n", flag1);
	
    }  // end of while time independant loop



  	for (int j=1 ; j < ry-1 ;j++)  {	// Tolerance check
	    for (int i=1 ; i < cx-1 ;i++) {
		v[j][i]=vo[j][i]+(Pa[j][i]-Pa[j+1][i])*dt/dy;
		u[j][i]=uo[j][i]+(Pa[j][i]-Pa[j][i+1])*dt/dx;  }  }

  } // end for dt loop
	

end = clock();
time_spent = (double)(end - begin)/ CLOCKS_PER_SEC;



	FILE *f = fopen("uoS", "w");
	fprintf(f, "run time = %f " , time_spent);
	fprintf(f, "uo=[ ");
	for (int j = 0; j < ry; j++) {
	   for (int i=0 ; i < cx ;i++) {
		fprintf(f, "%.5f ", uo[j][i] ); } 
		fprintf(f, "; " ) ; }
		fprintf(f, "]; \n ");
		fprintf(f, "surf(uo) \n ");
	  
	fclose(f);


printf("count = %d count1= %d and time spent = %f  \n", count , count1, time_spent);
return 0;
}

