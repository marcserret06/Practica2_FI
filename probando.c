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
int DiagonalDom( float M[N][N] ) {
float dom;
int res;
int dd = 1;
float sum = 0;
for (int i=0;i<N;i++) {
for (int j=0;j<N;j++) {
if (i == j) {
dom = fabs(M[i][j]);
}
else {
sum = sum + fabs(M[i][j]);
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
void PrintVect( float vect[N], int from, int numel ) {
int to = from + numel;
to++;
for (int i = from; i < to; i++) {
printf("%f ", vect[i]);
}
}
int Jacobi(float M[N][N], float vect[N], float vectres[N], unsigned iter) {
int dd = DiagonalDom(M);
if (dd == 0) {
return 0;
}
for (int m = 0; m<N; m++) {
vectres[m] = 0;
}
float v_aux[N];
for (int n = 0; n < iter; n++) {
float sum;
float diag;
float prod;
for (int i = 0; i<N; i++) {
sum = 0;
diag = 0;
prod = 0;
for (int j = 0; j<N; j++) {
if (i != j) {
prod = M[i][j] * vectres[j];
sum += prod;
}
else {
diag = M[i][i];
}
}
v_aux[i] = ((vect[i] - sum)/diag);
}
}
for (int k = 0;k<N;k++) {
vectres[k] = v_aux[k];
}
return 1;
}

int main() {
InitData();
float Jc1[N], Jc2[N], Jc3[N];
int x1, x2, x3;
printf("Empezamos\n");
x1 = Jacobi(MatDD, V3, Jc1, 1);
printf("x1 completado");
x2 = Jacobi(MatDD, V3, Jc2, 1000);
printf("x2 completado");
x3 = Jacobi(Mat, V3, Jc3, 1000);

if (x1 == 0) {
printf("La matriz MatDD no es Diagonal Dominante, no se puede aplicar\n");
}
else {
PrintVect(Jc1, 0, 9);
}
printf("\n");

if (x2 == 0) {
printf("La matriz MatDD no es Diagonal Dominante, no se puede aplicar\n");
}
else {
PrintVect(Jc2, 0, 9);
}
printf("\n");

if (x3 == 0) {
printf("La matriz Mat no es Diagonal Dominante, no se puede aplicar\n");
}
else {
PrintVect(Jc3, 0, 9);
}
printf("\n");
}
