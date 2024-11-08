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
printf("%f\n", mat[row][i]);
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

float Infininorm( float M[N][N] ) {
float suma1=0;
float sfila[N];
for (int i=0;i<N;i++) {
for (int n=0;n<N;n++) {
suma1 = suma1 + fabs(M[i][n]);
}
sfila[i] = suma1;
suma1=0;
}
float max=0;
for (int i=0;i<N;i++) {
if (sfila[i] > max) {
max = sfila[i];
}
}
return max;
}

float Onenorm(float M[N][N]) {
float suma2=0;
float scolumna[N];
for (int i=0;i<N;i++) {
for (int n=0;n<N;n++) {
suma2 = suma2 + fabs(M[n][i]);
}
scolumna[i] = suma2;
suma2=0;
}
float max=0;
for (int i=0;i<N;i++) {
if (scolumna[i] > max) {
max = scolumna[i];
}
}
return max;
}
float NormFrobenius( float M[N][N] ) {
float sumatot = 0;
float resul;
for (int i = 0; i < N; i++) {
for (int n = 0; n < N; n++) {
sumatot = sumatot + M[i][n]*M[i][n];
}
}
resul = sqrt(sumatot);
return resul;
}

int DiagonalDom( float M[N][N] ) {
float dom;
int res;
int dd = 1;
float sum = 0;
for (int i=0;i<N;i++) {
for (int j=0;j<N;j++) {
if (i == j) {
dom = abs(M[i][j]);
}
else {
sum = sum + abs(M[i][j]);
}
}
if (dom < sum) {
dd = 0;
break;
}
sum = 0;
}
if (dd == 0) {
res = 0;
}
else {
res = 1;
}
return res;
}

void Matriu_x_Vector( float M[N][N], float vect[N], float vectres[N] ) {
for (int i=0; i<N; i++) {
int sum = 0;
int prod = 0;
for (int j=0; j<N; j++) {
prod = M[i][j] * vect[j];
sum = sum + prod;
vectres[i] = sum;
}
}
}


int main() {
InitData();
float vectr[N];
Matriu_x_Vector(Mat, V2, vectr);
PrintVect(vectr, 0, 9);
}
