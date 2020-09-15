/*  Author: Kareem Abdelshafy. Date: 04/14/2016
		 High Performance computing Final Project
Parallel version of navier stoke solver for a fluid flow inside a channel   */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <algorithm>
#include "mpi.h"

//using namespace std;

int main ( int argc, char **argv)
{
	
	double rho = 1; double Re = 10; double Nx = 128-2; double Ny = 64-2; 
	double dt = 0.0005; double tend = 0.5; double L = 2; double h = 1;
	double dx = L/Nx;  double dy = h/Ny; double B=dx/dy;
	double w=1; 
	int ry= Ny+2;  int cx = Nx+2;

	double **out = (double **)malloc(ry * sizeof(double *));  // uo
	for (int j=0; j<ry; j++)  
		{out[j] = (double *)malloc(cx * sizeof(double));}

	for (int j=0 ; j < ry ;j++)  {	//C
	    for (int i=0 ; i < cx ;i++) {
		out[j][i]=0.5 ; } }

/////////////////////////////MPI vector and task preparation////////////////////////
	double start, stop;
	int size, myrank;

	MPI_Init (&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


///////////////////////////////
   int flag1=1 ; int flag=1; double R=0; double R1=0;
	 int count=0; int count1=0; double Rtot = 0;
	double R1tot = 0; double R1temp = 0; double Rtemp = 0; 
//////////////////////////////// Algorithm to split the balance the work load /////////////
	
	int mypart;
	int othersparts[size] ;

	for (int i=0 ; i < size ; i++ ) {
	  othersparts[i] = 0 ;
	}
	
   	int reminder= (ry % size);
	if (myrank < reminder)
	 { mypart =  ceil (float(ry) / float(size))  ;} 	
	else
	 { mypart = floor( float(ry) / float(size) ) ;} 

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Allgather ( &mypart, 1 , MPI_INT , othersparts, 1 , MPI_INT , MPI_COMM_WORLD) ;
	
	MPI_Barrier(MPI_COMM_WORLD);


	int leftColInd =0;

	for (int i=0 ; i < myrank ; i++ ) {
	  leftColInd += othersparts[i] ;
	}
	
	int rightColInd = leftColInd + mypart -1;   
	
	int mysize ;
	
	if (myrank == 0 || myrank == size-1) { mysize = mypart + 1; }
	else 
	{ mysize = mypart+2 ; }

	int rightNeighbor = ( myrank )+1;
	int leftNeighbor = ( myrank ) - 1 ;

	 

////////////////////////////// Memory allocation //////////////////////////////
 

	double **u = (double **)malloc(mysize * sizeof(double *));  // u
	for (int j=0; j<mysize; j++)  
		{u[j] = (double *)malloc(cx * sizeof(double));}

	double **v = (double **)malloc(mysize * sizeof(double *));  // v
	for (int j=0; j<mysize; j++)  
		{v[j] = (double *)malloc(cx * sizeof(double));}

	double **P = (double **)malloc(mysize * sizeof(double *));  // P
	for (int j=0; j<mysize; j++)  
		{P[j] = (double *)malloc(cx * sizeof(double));}

	double **uo = (double **)malloc(mysize * sizeof(double *));  // uo
	for (int j=0; j<mysize; j++)  
		{uo[j] = (double *)malloc(cx * sizeof(double));}

	double **vo = (double **)malloc(mysize * sizeof(double *));  // vo
	for (int j=0; j<mysize; j++)  
		{vo[j] = (double *)malloc(cx * sizeof(double));}

	double **C = (double **)malloc(mysize * sizeof(double *));  // C
	for (int j=0; j<mysize; j++)  
		{C[j] = (double *)malloc(cx * sizeof(double));}

	double **Pa = (double **)malloc(mysize * sizeof(double *));  // Pa
	for (int j=0; j<mysize; j++)  
		{Pa[j] = (double *)malloc(cx * sizeof(double));}

	double **Pb = (double **)malloc(mysize * sizeof(double *));  // Pb
	for (int j=0; j<mysize; j++)  
		{Pb[j] = (double *)malloc(cx * sizeof(double));}

	double **Po = (double **)malloc(mysize * sizeof(double *));  // Po
	for (int j=0; j<mysize; j++)  
		{Po[j] = (double *)malloc(cx * sizeof(double));}  

//////////////////// variable initialization /////////////////////////////////////////

	for (int j=0 ; j < mysize ;j++)  {	//u
	    for (int i=0 ; i < cx ;i++) {
		u[j][i]=0 ; } }

	for (int j=0 ; j < mysize ;j++) {	//u
		u[j][0]=1 ; } 

	for (int j=0 ; j < mysize ;j++)  {	//v
	    for (int i=0 ; i < cx ;i++) {
		v[j][i]=0 ; } }

	for (int j= 0 ; j < mysize ;j++)  {	//P
	    for (int i=0 ; i < cx ;i++) {
		P[j][i]=0 ; } }

	for (int j=0 ; j < mysize ;j++) {	//Pa
	    for (int i=0 ; i < cx ;i++) {
		Pa[j][i]=0 ; } }

	for (int j=0 ; j < mysize ;j++)  {	//Pb
	    for (int i=0 ; i < cx ;i++) {
		Pb[j][i]=0 ; } }

	for (int j=0 ; j < mysize ;j++)  {	//Po
	    for (int i=0 ; i < cx ;i++) {
		Po[j][i]=0 ; } }

	for (int j=0 ; j < mysize ;j++)  {	//uo
	    for (int i=0 ; i < cx ;i++) {
		uo[j][i]=u[j][i] ; } }

	for (int j=0 ; j < mysize ;j++)  {	//C
	    for (int i=0 ; i < cx ;i++) {
		C[j][i]=0 ; } }




MPI_Barrier(MPI_COMM_WORLD);
start = MPI_Wtime();
/////////////////////////////////  start of computation /////////////////////////////////////

	MPI_Request request;  MPI_Status status; 

for ( double t=0 ; t < tend ; t+=dt )  { 

	flag1=1 ; flag=1;

	while (flag1 == 1) {

	 count1++ ;


///////////////////////////////  update uo  /////////////////
	
	for (int j= 1 ; j < mysize-1 ;j++)  {	//uo
	    for (int i=1 ; i < cx-1 ;i++) {

		uo[j][i]=u[j][i] - dt/dx/4* ( pow(u[j][i]+u[j][i+1],2) - pow(u[j][i]+u[j][i-1],2) )

                -dt/dy/4* ( (u[j][i]+u[j+1][i]) * (v[j][i]+v[j][i+1]) - (u[j][i]+u[j-1][i]) * (v[j-1][i]+v[j-1][i+1]) )

                -( P[j][i+1]-P[j][i] )*dt/dx +

		 dt/Re*( (u[j][i+1]-2*u[j][i]+u[j][i-1])/(dx*dx) + (u[j+1][i]-2*u[j][i]+u[j-1][i])/(dy*dy) );   } }


	for (int j= 0 ; j < mysize ;j++) {	//uo left boundary
		uo[j][0]=1 ;  } 
	
	for (int j=0 ; j < mysize ;j++) {	//uo  right boundary
		uo[j][cx-1]=uo[j][cx-2] ; } 

	if (myrank == 0 ){
	for (int i=0 ; i < cx ;i++) {	//uo lower boundary
		uo[0][i]=-1*uo[1][i] ; }  }

	if (myrank == size-1 ){
	for (int i=0 ; i < cx ;i++) {	//uo  upper boundary
		uo[mysize-1][i]=-1*uo[mysize-2][i] ; }  }


	if (myrank < size-1)
	{
	MPI_Isend( & (uo[mysize-2][0]) , cx , MPI_DOUBLE, rightNeighbor, 0, MPI_COMM_WORLD, & request);
	MPI_Recv( & (uo[mysize-1][0])  , cx , MPI_DOUBLE , rightNeighbor , 1 , MPI_COMM_WORLD, & status) ;
	}

	if (myrank > 0)
	{
	MPI_Isend( & (uo[1][0]) , cx , MPI_DOUBLE, leftNeighbor, 1, MPI_COMM_WORLD, & request);
	MPI_Recv( & (uo[0][0]) , cx , MPI_DOUBLE , leftNeighbor , 0 , MPI_COMM_WORLD, & status) ;
	}


///////////////////////////////  update vo //////////

	for (int j=1 ; j < mysize-1 ;j++)  {	//vo
	    for (int i=1 ; i < cx-1 ;i++) {

            vo[j][i]=v[j][i] -dt/dy/4 * ( pow(v[j][i]+v[j+1][i],2) - pow(v[j][i]+v[j-1][i],2) )

                -dt/dx/4* ( (u[j][i]+u[j+1][i])*(v[j][i]+v[j][i+1]) - (u[j][i-1]+u[j+1][i-1])*(v[j][i]+v[j][i-1]) )

                -(P[j+1][i]-P[j][i])*dt/dy

		+dt/Re*( (v[j][i+1]-2*v[j][i]+v[j][i-1])/(dx*dx) +(v[j+1][i]-2*v[j][i]+v[j-1][i])/(dy*dy) ) ; } }

	for (int j= 0 ; j < mysize ;j++) {	//vo left boundary
		vo[j][0]=-1*vo[j][1] ; } 

	if (myrank == 0 ){
	for (int i=0 ; i < cx ;i++) {	//vo lower boundary
		vo[0][i]=0 ; }   }

	if (myrank == size-1 ){
	for (int i=0 ; i < cx ;i++) {	//vo  upper boundary  tricky !!
		vo[mysize-2][i]=0 ; }  }

	for (int j=0 ; j < mysize ;j++) {	//vo  right boundary
		vo[j][cx-1]=vo[j][cx-2] ; } 


	if (myrank < size-1)
	{
	MPI_Isend( & (vo[mysize-2][0]) , cx , MPI_DOUBLE, rightNeighbor, 2, MPI_COMM_WORLD, & request);
	MPI_Recv( & (vo[mysize-1][0])  , cx , MPI_DOUBLE , rightNeighbor , 3 , MPI_COMM_WORLD, & status) ;
	}

	if (myrank > 0)
	{
	MPI_Isend( & (vo[1][0]) , cx , MPI_DOUBLE, leftNeighbor, 3, MPI_COMM_WORLD, & request);
	MPI_Recv( & (vo[0][0]) , cx , MPI_DOUBLE , leftNeighbor , 2 , MPI_COMM_WORLD, & status) ;
	}

///////////////////////// update C ////////////

	for (int j=1 ; j < mysize-1 ;j++)  {	//C
	    for (int i=1 ; i < cx-1 ;i++) {
		C[j][i]= (uo[j][i]-uo[j][i-1])*dx/dt+(vo[j][i]-vo[j-1][i])*B*B*dy/dt ; } }

	if (myrank < size-1)
	{
	MPI_Isend( & (C[mysize-2][0]) , cx , MPI_DOUBLE, rightNeighbor, 4, MPI_COMM_WORLD, & request);
	MPI_Recv( & (C[mysize-1][0])  , cx , MPI_DOUBLE , rightNeighbor , 5 , MPI_COMM_WORLD, & status) ;
	}

	if (myrank > 0)
	{
	MPI_Isend( & (C[1][0]) , cx , MPI_DOUBLE, leftNeighbor, 5, MPI_COMM_WORLD, & request);
	MPI_Recv( & (C[0][0]) , cx , MPI_DOUBLE , leftNeighbor , 4 , MPI_COMM_WORLD, & status) ;
	}


//////////////////////// LSOR  Poission Equation ////////////////////
////////////////////////////////////////////////////////////////////

	while (flag == 1) {
	
	//count++ ;

	for (int j=0 ; j < mysize ;j++)  {	// Pt = Pa
	    for (int i=0 ; i < cx ;i++) {
		Pb[j][i]=Pa[j][i];  } }

	for (int j=1 ; j < mysize-1 ;j++)  {	//C
	    for (int i=1 ; i < cx-1 ;i++) {
                Pa[j][i]=0.25*(Pb[j][i+1]+Pb[j][i-1]+B*B*(Pb[j+1][i]+Pb[j-1][i])-C[j][i]);  } }

	if (myrank < size-1)
	{
	MPI_Isend( & (Pa[mysize-2][0]) , cx , MPI_DOUBLE, rightNeighbor, 6, MPI_COMM_WORLD, & request);
	MPI_Recv( & (Pa[mysize-1][0])  , cx , MPI_DOUBLE , rightNeighbor , 7 , MPI_COMM_WORLD, & status) ;
	}

	if (myrank > 0)
	{
	MPI_Isend( & (Pa[1][0]) , cx , MPI_DOUBLE, leftNeighbor, 7, MPI_COMM_WORLD, & request);
	MPI_Recv( & (Pa[0][0]) , cx , MPI_DOUBLE , leftNeighbor , 6 , MPI_COMM_WORLD, & status) ;
	}
	
	 R = 0;   Rtot =0;   Rtemp=0;

	for (int j=0 ; j < myrank ;j++)  {	// Tolerance check
	    for (int i=0 ; i < cx ;i++) {
		 Rtemp= fabs((Pa[j][i]-Pb[j][i]) / (Pb[j][i]+0.0000001));
		R= std::max(R,Rtemp);	 } }


	MPI_Allreduce ( & R, & Rtot, 1 ,  MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD) ;
	if ( Rtot > 0.001) { flag = 1;}   else  {flag = 2;}
	
	MPI_Barrier(MPI_COMM_WORLD);

	}

	flag = 1;

//////////////////////////////////////////////////////////////////////
///////////////////////  end of LSOR   //////////////////////////////
	

	for (int j=0 ; j < mysize ;j++)  {	// Pt = P
	    for (int i=0 ; i < cx ;i++) {
		Po[j][i]=P[j][i];  } }


	for (int j= 1 ; j < mysize-1 ;j++)  {  // P = P + Pa
	    for (int i=1 ; i < cx-1 ;i++) {
		P[j][i]=P[j][i] + Pa[j][i];  } }


	if (myrank < size-1)
	{
	MPI_Isend( & (P[mysize-2][0]) , cx , MPI_DOUBLE, rightNeighbor, 8, MPI_COMM_WORLD, & request);
	MPI_Recv( & (P[mysize-1][0])  , cx , MPI_DOUBLE , rightNeighbor , 9 , MPI_COMM_WORLD, & status) ;
	}

	if (myrank > 0)
	{
	MPI_Isend( & (P[1][0]) , cx , MPI_DOUBLE, leftNeighbor, 9, MPI_COMM_WORLD, & request);
	MPI_Recv( & (P[0][0]) , cx , MPI_DOUBLE , leftNeighbor , 8 , MPI_COMM_WORLD, & status) ;
	}

	R1=0;  R1tot=0;  R1temp=0; 


	for (int j=0 ; j < mysize ;j++)  {	// Tolerance check
	    for (int i=0 ; i < cx ;i++) {
		 R1temp= fabs((P[j][i]-Po[j][i]) / (Po[j][i]+0.0000001));
		R1= std::max(R1,R1temp);			} }


	MPI_Allreduce ( & R1, & R1tot, 1 ,  MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD) ;
	if ( R1tot > 0.00001) { flag1 = 1;}   else  {flag1 = 2;}
	
	MPI_Barrier(MPI_COMM_WORLD);

	
    }  // end of while time independant loop


	for (int j=1 ; j < mysize-1 ;j++)  {	
	    for (int i=1 ; i < cx-1 ;i++) {
		v[j][i]=vo[j][i]+(Pa[j][i]-Pa[j+1][i])*dt/dy;
		u[j][i]=uo[j][i]+(Pa[j][i]-Pa[j][i+1])*dt/dx;  }  }

	if (myrank < size-1)
	{
	MPI_Isend( & (u[mysize-2][0]) , cx , MPI_DOUBLE, rightNeighbor, 10, MPI_COMM_WORLD, & request);
	MPI_Recv( & (u[mysize-1][0])  , cx , MPI_DOUBLE , rightNeighbor , 11 , MPI_COMM_WORLD, & status) ;
	MPI_Isend( & (v[mysize-2][0]) , cx , MPI_DOUBLE, rightNeighbor, 12, MPI_COMM_WORLD, & request);
	MPI_Recv( & (v[mysize-1][0])  , cx , MPI_DOUBLE , rightNeighbor , 13 , MPI_COMM_WORLD, & status) ;
	}

	if (myrank > 0)
	{
	MPI_Isend( & (u[1][0]) , cx , MPI_DOUBLE, leftNeighbor, 11, MPI_COMM_WORLD, & request);
	MPI_Recv( & (u[0][0]) , cx , MPI_DOUBLE , leftNeighbor , 10 , MPI_COMM_WORLD, & status) ;
	MPI_Isend( & (v[1][0]) , cx , MPI_DOUBLE, leftNeighbor, 13, MPI_COMM_WORLD, & request);
	MPI_Recv( & (v[0][0]) , cx , MPI_DOUBLE , leftNeighbor , 12 , MPI_COMM_WORLD, & status) ;
	}

printf ("myrank = %d , time is %f \n" , myrank, t);


} // end  of time dependant loop for (t=0 -> tend)

MPI_Barrier(MPI_COMM_WORLD);
stop = MPI_Wtime();
double timespent= stop - start;
/////////////////////  Writing the data to output file  :( /////////////////////////////
MPI_Status status99, status98;


if (myrank !=0 )
{

	MPI_Ssend( & mysize , 1 , MPI_DOUBLE, 0, 98, MPI_COMM_WORLD);
	MPI_Ssend( & (uo[0][0]) , (mysize)*(cx+2) , MPI_DOUBLE, 0, 99 , MPI_COMM_WORLD);

}
else
{
 	FILE *f = fopen("uo_MPI", "w");
	fprintf (f, "time_spent = %f ;\n", timespent);
	fprintf (f, "no_of_processors = %d ; \n", size);
	fprintf(f, "uo=[ ");
	for (int j = 0; j < mysize; j++) 	{
	   for (int i=0 ; i < cx ;i++) 		 {
		fprintf(f , "%.3f ", uo[j][i] );     } 
		fprintf(f, " ; "); 		 }



	for ( int n= 1 ; n < size ; n++ ) {
		
	MPI_Recv( & mysize , 1 , MPI_DOUBLE, n, 98, MPI_COMM_WORLD, &status98);
	MPI_Recv( & (out[0][0]) , (mysize)*(cx+2) , MPI_DOUBLE, n, 99 , MPI_COMM_WORLD, &status99 );

	for (int j = 0; j < (mysize); j++) 	{
	   for (int i=0 ; i < cx ;i++) 		 {
		fprintf(f ,"%.3f ", out[j][i] );     } 
		fprintf(f ," ; "); 		 }

					 }
 		fprintf(f, "]; \n ");
		fprintf(f, "surf(uo) \n ");
		fclose(f);
}	



//////////////////////////////////////////////////////////////////////////////////////////

MPI_Finalize();

return 0;
}

