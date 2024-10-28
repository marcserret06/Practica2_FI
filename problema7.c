#define N 512
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float V1[N], V2[N], V3[N];
float Mat[N][N], MatDD[N][N];
void InitData(){
int i,j;
srand(334411);
for( i = 0; i < N; i++ )
 for( j = 0; j < N; j++ ){
 Mat[i][j]=(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
 if ( (abs(i - j) <= 3) && (i != j))
 MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0*RAND_MAX));
 else if ( i == j )
 MatDD[i][j]=(((i*j)%3)?-1:1)*(10000.0*(rand()/(1.0*RAND_MAX)));
 else MatDD[i][j] = 0.0;
 }
for( i = 0; i < N; i++ ){
 V1[i]=(i<N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
 V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
 V3[i]=(((i*j)%5)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
}
}
void PrintVect( float vect[N], int from, int numel ) {
int to = from + numel;
to++;
for (int i = from; i < to; i++) {
printf("%f ", vect[i]);
}
}
void PrintRow( float mat[N][N], int row, int from, int numel ) {
int to = from + numel;
to++;
for (int i=from;i<to;i++) {
printf("%f ", mat[row][i]);
}
}
void MultEscalar( float vect[N], float vectres[N], float alfa ) {
for (int i = 0; i < N; i++) {
vectres[i] = vect[i] * alfa;
}
}

float Scalar( float vect1[N], float vect2[N] ) {
float resul=0;
float nvec;
for (int i = 0; i < N; i++) {
nvec = vect1[i] * vect2[i];
resul = resul + nvec;
}
return resul;
}

float Magnitude( float vect[N] ) {
float res, suma;
suma = Scalar(vect, vect);
res=sqrt(suma);
return res;
}

int Ortogonal( float vect1[N], float vect2[N] ) {
int res;
float Ortog;
Ortog = Scalar(vect1, vect2);
if (Ortog == 0) {
res = 1;
}
else {
res = 0;
}
return res;
}

void  Projection( float vect1[N], float vect2[N], float vectres[N] ) {
float magv2, sca, div;
magv2=Magnitude(vect2);
sca=Scalar(vect1, vect2);
div=sca/magv2;
MultEscalar(vect2, vectres, div);
}

int main() {
float Vres[N];
InitData();
Projection(V2, V3, Vres);
PrintVect(Vres, 0, 10);
}
