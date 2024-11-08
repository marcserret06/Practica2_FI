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
resul = sqrt(sumatot); return resul;
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
for (int k = 0;k<N;k++) {
vectres[k] = v_aux[k];
}
}
return 1;
}

int main() {
InitData();
// Resultados A
printf("Los valores entre 0 y 9, y entre 256 y 265 del V1 son:\n");
PrintVect(V1, 0, 9);
printf("\n");
PrintVect(V1, 256, 9);
printf("\n");
printf("Los valores entre 0 y 9, y entre 256 y 265 del V2 son:\n");
PrintVect(V2, 0, 9);
printf("\n");
PrintVect(V2, 256, 9);
printf("\n");
printf("Los valores entre 0 y 9, y entre 256 y 265 del V2 son:\n");
PrintVect(V3, 0, 9);
printf("\n");
PrintVect(V3, 256, 9);
printf("\n\n");
// Resultados B
printf("Los valores entre 0 y 9, de las filas 0 y 100 de Mat son:\n");
PrintRow(Mat,0,0,9);
printf("\n");
PrintRow(Mat,100,0,9);
printf("\n\n");
//Resultados C
printf("Los valores entre 0 y 9 de la fila 0 y desde el 95 al 104 de la fila 100 de MatDD son:\n");
PrintRow(MatDD,0,0,9);
printf("\n");
PrintRow(MatDD,100,95,9);
printf("\n\n");
//Resultado D
float a1, b1, a2, b2, a3, b3;
int a4, b4;
a1 = Infininorm(Mat);
a2 = Onenorm(Mat);
a3 = NormFrobenius(Mat);
a4 = DiagonalDom(Mat);
b1 = Infininorm(MatDD);
b2 = Onenorm(MatDD);
b3 = NormFrobenius(MatDD);
b4 = DiagonalDom(MatDD);
printf("La infininorma de la matriz Mat = %f\n", a1);
printf("La norma uno de la matriz Mat = %f\n", a2);
printf("La nomra de Frobenius de la matriz Mat = %f\n", a3);
if (a4==1) {
printf("La matriz Mat es Diagonal Dominante\n");
}
else {
printf("La matriz Mat no es Diagonal Dominante\n");
}
printf("\n");
printf("La infininorma de la matriz MatDD = %f\n", b1);
printf("La norma uno de la matriz MatDD = %f\n", b2);
printf("La nomra de Frobenius de la matriz MatDD = %f\n", b3);
if (b4==1) {
printf("La matriz MatDD es Diagonal Dominante\n");
}
else {
printf("La matriz MatDD no es Diagonal Dominante\n");
}
printf("\n\n");
//Resultado E
float v1v2, v2v3, v1v3;
v1v2 = Scalar(V1, V2);
v1v3 = Scalar(V1, V3);
v2v3 = Scalar(V2, V3);
printf("El producto escalar entre V1 y V2 es: %f\nEl producto escalar entre V1 y V3 es: %f\nEl producto escalar entre V2 y V3 es: %f\n", v1v2, v1v3, v2v3);
printf("\n");
//Resultado F
float mgv1, mgv2, mgv3;
mgv1 = Magnitude(V1);
mgv2 = Magnitude(V2);
mgv3 = Magnitude(V3);
printf("Las magnitudes de los vectores V1, V2 y V3 son: %f, %f, %f\n", mgv1, mgv2, mgv3);
printf("\n");
//Resultado G
int Ortv12;
Ortv12 = Ortogonal(V1, V2);
if (Ortv12 == 1) {
printf("Los vectores V1 y V2 son ortogonales\n");
} else {
printf("Los vectores V1 y V2 no son ortogonales\n");
}
printf("\n");
//Resultado H
float vecMS3[N];
MultEscalar(V3, vecMS3, 2);
printf("Los valores entre 0 y 9, y 256 y 265 del V3 al multiplicarlo por el escalar 2 es:\n");
PrintVect(vecMS3, 0, 9);
printf("\n");
PrintVect(vecMS3, 256, 9);
printf("\n\n");

//Resultado I

float Projv2v3[N], Projv1v2[N];
Projection(V2, V3, Projv2v3);
Projection(V1, V2, Projv1v2);
printf("Los elementos del 0 al 9 de la proyección del V2 sobre V3 son:\n");
PrintVect(Projv2v3, 0, 9);
printf("\n");
printf("Los elementos del 0 al 9 de la proyección del V1 sobre V2 son:\n");
PrintVect(Projv1v2, 0, 9);
printf("\n\n");

//Resultado J

float Mmatv2[N];
Matriu_x_Vector(Mat, V2, Mmatv2);
printf("Los valores entre el 0 y el 9 resultantes de multiplicar la Matriz Mat por el vector V2 son:\n");
PrintVect(Mmatv2, 0, 9);
printf("\n\n");
// Resultado K
float Jc1[N], Jc2[N], Jc3[N];
int x1, x2, x3;
x1 = Jacobi(MatDD, V3, Jc1, 1);
x2 = Jacobi(MatDD, V3, Jc2, 1000);
x3 = Jacobi(Mat, V3, Jc3, 1000);
if (x1 == 0) {
printf("La matriz MatDD no es Diagonal Dominante, no se puede aplicar\n");
}
else {
printf("Los primeros 10 valores del vector resultante ´x´ de la Matriz MatDD por el vector V3 tras una iteración son:\n");
PrintVect(Jc1, 0, 9);
}
printf("\n");

if (x2 == 0) {
printf("La matriz MatDD no es Diagonal Dominante, no se puede aplicar\n");
}
else {
printf("Los primeros 10 valores del vector resultante ´x´ de la Matriz MatDD por el vector V3 tras mil iteraciones son:\n");
PrintVect(Jc2, 0, 9);
}
printf("\n");

if (x3 == 0) {
printf("La matriz Mat no es Diagonal Dominante, no se puede aplicar\n");
}
else {
printf("Los primeros 10 valores del vector resultante ´x´ de la Matriz Mat por el vector V3 tras una iteración son:\n");
PrintVect(Jc3, 0, 9);
}
printf("\n");
}
