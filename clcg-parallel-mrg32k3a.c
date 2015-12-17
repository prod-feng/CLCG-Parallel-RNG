#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

/*
Parallel Combined Linear Congruential MRG32K3A random number generator. OpenMP safe. 
*/
#define Aa  1403580
#define Bb  -810728
#define Cc  527612
#define Dd  -1370589


static unsigned long long PMOD1        = 4294967087;
static unsigned long long PMOD2        = 4294944443;

long double  matrix1[3][3] = { {0,Aa,Bb}, {1,0,0}, {0,1,0} };

long double  matrix2[3][3] = {{ Cc,0,Dd}, {1,0,0}, {0,1,0} };

long double  m1[3][3] = { {0,Aa,Bb}, {1,0,0}, {0,1,0} };

long double  m2[3][3] = {{ Cc,0,Dd}, {1,0,0}, {0,1,0} };

/*default seeds for X and Y*/

long double  xseeds[3]={1,0,0};
long double  yseeds[3]={1,0,0};

#define MAX_THREADS 128

#define DIMENSION 1 //3?

long double rseedX[3][DIMENSION];
long double rseedY[3][DIMENSION];

static int first[MAX_THREADS][DIMENSION][4];

//static  unsigned long long tmpseedX[MAX_THREADS][3][DIMENSION];

//static  unsigned long long tmpseedY[MAX_THREADS][3][DIMENSION];

long double tmpseedX[MAX_THREADS][3][DIMENSION];
long double tmpseedY[MAX_THREADS][3][DIMENSION];


long double xyz[DIMENSION];

#pragma omp threadprivate(rseedX,rseedY,xyz)

int matrixMult(int np,long double *mm, long double *matrix,int col,unsigned long long PMOD, long double *retmatrix)
{
  int i,j,k,steps;
  int64_t mod;

  long double tmp1,mult_tmp;
  long double sign,ret1[3][3];
  for(steps=0; steps <np; steps++)
  { 
    for(i=0;i<3;i++){
      for(j=0; j<col;j++){
        tmp1 = 0;
        for(k=0; k<3;k++){
           /*calculate matrix*/
           mult_tmp = (mm[i*3+k]*matrix[k*col+j]);
           mod = mult_tmp/PMOD;
           mult_tmp = mult_tmp -(long double)mod * (long double)PMOD;
           if (mult_tmp < 0.0 && col <=1) mult_tmp += (long double)PMOD;

           tmp1 += mult_tmp;
        }
        mod =(tmp1)/PMOD;
        tmp1 = tmp1 - (long double)mod * (long double)PMOD;
        if (tmp1 < 0.0 && col <=1) tmp1 += (long double)PMOD;
        ret1[i][j] = tmp1;
      }
    }
    /*update the state matrix*/
    for(i=0;i<3;i++){
      for(j=0; j<col;j++){
        retmatrix[i*col+j] = ret1[i][j];
      }
    }
  }
// loop to jump ahead np steps to calculate the steps to jump.
}


int clcgrandom(int num, long double * ret)
{
   unsigned long long random_next;
   int i,j,sign;

   long double  tmpvarX[3],tmpvarY[3];
   long double tmpsum;
   int64_t mod;

   if(num > DIMENSION) num = DIMENSION;
   int id=omp_get_thread_num();

   for(j=0; j< num; j++){
     //re-use the first random numbers generated for different threads.
     if(id > 0 && first[id][j][0] == 1){
       first[id][j][0]=0;
     }
     else
     {

         tmpvarX[0] = rseedX[0][j];
         tmpvarX[1] = rseedX[1][j];
         tmpvarX[2] = rseedX[2][j];

         tmpvarY[0] = rseedY[0][j];
         tmpvarY[1] = rseedY[1][j];
         tmpvarY[2] = rseedY[2][j];

         matrixMult(1,&m1[0][0],&tmpvarX[0],1,PMOD1,&tmpvarX[0]);
         matrixMult(1,&m2[0][0],&tmpvarY[0],1,PMOD2,&tmpvarY[0]);


         rseedX[0][j] = tmpvarX[0];
         rseedX[1][j] = tmpvarX[1];
         rseedX[2][j] = tmpvarX[2];

         rseedY[0][j] = tmpvarY[0];
         rseedY[1][j] = tmpvarY[1];
         rseedY[2][j] = tmpvarY[2];

     }

     tmpsum = rseedX[0][j] + rseedY[0][j];
  
     mod =(tmpsum)/PMOD1;
     tmpsum = tmpsum - (long double)mod * (long double)PMOD1;
     if (tmpsum < 0.0) tmpsum += (long double)PMOD1;
  
      ret[j]=tmpsum/(long double)PMOD1;

      if(ret[j] < 0) ret[j] = ret[j]+1;
   }
   return num;
}

void setseed(int64_t *iseed)
{
   int i, j, id, nthreads;
   int k;

   long double tmpvarX[3],tmpvarY[3];


   id = omp_get_thread_num();
 
   #pragma omp single
   {

      nthreads = omp_get_num_threads();
      //initilize the state matrix.
      matrixMult(nthreads-1,&m1[0][0],&matrix1[0][0],3,PMOD1,&m1[0][0]);

      matrixMult(nthreads-1,&m2[0][0],&matrix2[0][0],3,PMOD2,&m2[0][0]);
      // Init the flags for all threads.
      for(j=0; j< DIMENSION; j++){
        for(i=0; i< MAX_THREADS; i++){
          first[i][j][0]=1;
        }
        //

        tmpvarX[0] = xseeds[0];
        tmpvarX[1] = xseeds[1];
        tmpvarX[2] = xseeds[2];

        tmpvarY[0] = yseeds[0];
        tmpvarY[1] = yseeds[1];
        tmpvarY[2] = yseeds[2];


        tmpseedX[0][0][j] = xseeds[0];
        tmpseedX[0][1][j] = xseeds[1];
        tmpseedX[0][2][j] = xseeds[2];

        tmpseedY[0][0][j] = yseeds[0];
        tmpseedY[0][1][j] = yseeds[1];
        tmpseedY[0][2][j] = yseeds[2];

        for (i = 1; i < nthreads; ++i)
        {
          matrixMult(1,&matrix1[0][0],&tmpvarX[0],1,PMOD1,&tmpvarX[0]);
          matrixMult(1,&matrix2[0][0],&tmpvarY[0],1,PMOD2,&tmpvarY[0]);

          /*store the seeds here*/
          tmpseedX[i][0][j]  = tmpvarX[0];
          tmpseedX[i][1][j]  = tmpvarX[1];
          tmpseedX[i][2][j]  = tmpvarX[2];

          tmpseedY[i][0][j]  = tmpvarY[0];
          tmpseedY[i][1][j]  = tmpvarY[1];
          tmpseedY[i][2][j]  = tmpvarY[2];

        }
      }
  }
  for(j=0; j< DIMENSION; j++){
     rseedX[0][j] =  tmpseedX[id][0][j];
     rseedX[1][j] =  tmpseedX[id][1][j];
     rseedX[2][j] =  tmpseedX[id][2][j];

     rseedY[0][j] =  tmpseedY[id][0][j];
     rseedY[1][j] =  tmpseedY[id][1][j];
     rseedY[2][j] =  tmpseedY[id][2][j];
  }
}


static long num_trials = 10000000;
int main ()
{
 long i,j; long Ncirc = 0; 
 double pi, x, y;

 int64_t myseed[3] = {1,0,0};

 int ndim=0;


#pragma omp parallel private (x, y) reduction (+:Ncirc)
  {
     setseed(myseed);
     #pragma omp for
     for(i=0;i<num_trials; i++)
     {
       ndim = clcgrandom(1,xyz);
       x=xyz[0];
       ndim = clcgrandom(1,xyz);
       y=xyz[0];
       //printf("XYZ= %llf   %d\n",xyz[0],omp_get_thread_num());  
       if (( x*x + y*y) <= 1.0) Ncirc++;
     }
  }
  pi = 4.0 * ((double)Ncirc/(double)num_trials);
  printf("\n %d trials, pi = %f \n",num_trials, pi);
}
