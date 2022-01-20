#define ARRAY(MAT) &(MAT[0][0])

int** buildMatrix(int N);
void destroyMatrix(int** mat);
int** readMatrix(int N);
void fillMatrix(int** mat, int value, int N);
void copyMatrix(int** source, int** dest, int N);
void printMatrix(int** mat, int N);
