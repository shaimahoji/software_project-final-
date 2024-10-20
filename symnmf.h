#ifndef SYMNMF_H
#define SYMNMF_H
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION


/*we want to represent the matrix as pointer of int pointres*/
typedef struct {
    int rows;
    int cols;
    float** mat;
}Matrix;

Matrix create_matrix(int rows, int cols);
//Matrix read_matrix(char filename[]);
//double euclid_distance(float *point1, float *point2, int n);
void print_matrix(Matrix matrix);
//void free_matrix(Matrix matrix);
Matrix sym(Matrix matrix);
//float sum(float *arr, int size);
Matrix ddg(Matrix matrix);
//Matrix multiply_matrices(Matrix A, Matrix B);
//Matrix norm(Matrix matrix);

#endif 