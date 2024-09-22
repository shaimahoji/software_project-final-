#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//we want to represent the matrix as pointer of int pointres

typedef struct {
    int rows;
    int cols;
    float** mat;
}Matrix;

//using the file create the matrix
Matrix read_matrix(char filename[]){
    FILE *file;
    float **matrix;
    int rows = 0, cols = 0;
    char line[256];
    Matrix mat;

    // Open the file for reading
    file = fopen(filename, "r");
    if (file == NULL) {
        return mat;
    }

    // First, determine the number of rows and columns
    while (fgets(line, sizeof(line), file)) {
        rows++;
        int current_cols = 0;
        char *token = strtok(line, ",");
        while (token) {
            current_cols++;
            token = strtok(NULL, ",");
        }
        if (current_cols > cols) {
            cols = current_cols;  // Update the maximum columns found
        }
    }

    // Allocate memory for the matrix
    matrix = (float **)malloc(rows * sizeof(float *));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (float *)malloc(cols * sizeof(float));
    }

    mat.cols = cols;
    mat.rows = rows;

    // Reset file pointer to the beginning of the file
    rewind(file);

    // Read the matrix values
    int row = 0;
    while (fgets(line, sizeof(line), file) && row < rows) {
        int col = 0;
        char *token = strtok(line, ",");
        while (token!=NULL) {
            matrix[row][col++] = atof(token);  // Convert string to int
            token = strtok(NULL, ",");
        }
        row++;
    }

    // Close the file
    mat.mat = matrix;
    fclose(file);

    return mat;
    }

double euclid_distance(float *point1, float *point2, int n) {
    double sum = 0.0;
    for (int i = 0 ; i<n ;i++) {
        double diff = point1[i]- point2[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

// Print the matrix
void print_matrix(Matrix matrix){
    for (int i = 0; i < matrix.rows; i++) {
        for (int j = 0; j < matrix.cols; j++) {
            if (j+1 == matrix.cols){printf("%.4f", matrix.mat[i][j]);}
            else{printf("%.4f,", matrix.mat[i][j]);}
        }
        printf("\n");
    }
}
// Free allocated memory
void free_matrix(Matrix matrix){
    for (int i = 0; i < matrix.rows; i++) {
        free(matrix.mat[i]);
    }
    free(matrix.mat);
}
// Function to create an empty matrix of size rows x cols
Matrix create_matrix(int rows, int cols) {
    Matrix matrix;
    matrix.rows = rows;
    matrix.cols = cols;
    matrix.mat = (float **)malloc(rows * sizeof(float *));
    for (int i = 0; i < rows; i++) {
        matrix.mat[i] = (float *)malloc(cols * sizeof(float));
    }
    return matrix;
}
//Compute the Similarity Matrix
Matrix sym(Matrix matrix){
    int n  = matrix.rows; // the number of data points
    Matrix sym = create_matrix(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if(i!=j){
                float dist = euclid_distance(matrix.mat[i],matrix.mat[j],matrix.cols);
                sym.mat[i][j]= exp(-(dist*dist)/2);
                }
            else{sym.mat[i][j]=0;}
        }
    }
    return sym;  
}
// the sum of two data points (a row in the matrix) 
float sum(float *arr, int size) {
    float total = 0.0;
    for (int i = 0; i < size; i++) {
        total += arr[i];  // Dereference the pointer to get the value
    }
    return total;
}
//Compute the diagonal degree Matrix
Matrix ddg(Matrix matrix){
    Matrix sym_mat = sym(matrix);
    int n  = sym_mat.rows; // the number of data points
    Matrix ddg = create_matrix(n,n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if(i==j){
                ddg.mat[i][j]= sum(sym_mat.mat[i],n);
                }
            else{ddg.mat[i][j]=0;}
        }
    }  
    return ddg;  
}
// Function to multiply matrices A and B, resulting in matrix C
Matrix multiply_matrices(Matrix A, Matrix B) {
    // Check if multiplication is possible (A.cols must equal B.rows)
    if (A.cols != B.rows) {
        printf("Matrix multiplication is not possible due to incompatible dimensions.\n");
        exit(EXIT_FAILURE);
    }

    // Create result matrix C with dimensions A.rows x B.cols
    Matrix C = create_matrix(A.rows, B.cols);

    // Perform the multiplication
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < B.cols; j++) {
            C.mat[i][j] = 0;  // Initialize to 0
            for (int k = 0; k < A.cols; k++) {
                C.mat[i][j] += A.mat[i][k] * B.mat[k][j];
            }
        }
    }

    return C;
}

//Compute the normalized similarity matrix
Matrix norm(Matrix matrix){
    Matrix D = ddg(matrix);
    Matrix inverse_sqrt_D = create_matrix(D.rows,D.cols);
    Matrix A = sym(matrix);
    Matrix norm; 

    for (int i=0; i<D.rows;i++){
        for (int j=0; j<D.cols; j++){
            
            if(i==j){
                inverse_sqrt_D.mat[i][j] = 1 / sqrt(D.mat[i][j]);}  
        }
    }

    norm = multiply_matrices(multiply_matrices(inverse_sqrt_D,A),inverse_sqrt_D);
    return norm;
}
int main(int argc, char *argv[]) {   
    char goal[256];
    char filename[256];

    strcpy(goal, argv[1]);
    strcpy(filename, argv[2]);

    Matrix matrix = read_matrix(filename);

    if (strcmp(goal, "sym") == 0){print_matrix(sym(matrix));}
    if (strcmp(goal, "ddg") == 0){print_matrix(ddg(matrix));}
    if (strcmp(goal, "norm") == 0){print_matrix(norm(matrix));}
}