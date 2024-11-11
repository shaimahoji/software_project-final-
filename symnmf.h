#ifndef SYMNMF_H
#define SYMNMF_H
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/* Header file for Symmetric Non-negative Matrix Factorization (SymNMF) functions
 * Defines the Matrix structure and declares functions for matrix operations
 */


/* Matrix struct to represent a 2D matrix.
 * - rows: Number of rows in the matrix
 * - cols: Number of columns in the matrix
 * - mat: Pointer to an array of double pointers, representing a 2D array of doubles
 */
typedef struct {
    int rows;
    int cols;
    double** mat;
}Matrix;

Matrix createMatrix(int rows, int cols);

void printMatrix(Matrix matrix);

Matrix sym(Matrix matrix);

Matrix ddg(Matrix matrix);

Matrix norm(Matrix matrix);

void freeMatrix(Matrix matrix);

Matrix converge(Matrix init_H,Matrix W);

#endif

// old comment: /* Representing the matrix as a pointer to int pointres*/
/*
Matrix create_matrix(int rows, int cols);
void print_matrix(Matrix matrix);
Matrix sym(Matrix matrix);
Matrix ddg(Matrix matrix);
Matrix norm(Matrix matrix);
void free_matrix(Matrix matrix);
Matrix converge(Matrix init_H,Matrix W);
*/
