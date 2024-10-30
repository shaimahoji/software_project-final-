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
void print_matrix(Matrix matrix);
Matrix sym(Matrix matrix);
Matrix ddg(Matrix matrix);
Matrix norm(Matrix matrix);
void free_matrix(Matrix matrix);

#endif 