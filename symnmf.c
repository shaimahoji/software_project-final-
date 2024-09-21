#include <stdio.h>
#include <stdlib.h>

//we want to represent the matrix as pointer of int pointres
//using the file create the matrix
float** read_matrix(char filename[]){
    FILE *file;
    double **matrix;
    int rows = 0, cols = 0;
    char line[256];

    // Open the file for reading
    file = fopen(filename, "r");
    if (file == NULL) {
        return EXIT_FAILURE;
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

    // Reset file pointer to the beginning of the file
    rewind(file);

    // Read the matrix values
    int row = 0;
    while (fgets(line, sizeof(line), file) && row < rows) {
        int col = 0;
        char *token = strtok(line, ",");
        while (token) {
            matrix[row][col++] = atoi(token);  // Convert string to int
            token = strtok(NULL, ",");
        }
        row++;
    }

    // Close the file
    fclose(file);
}

void print_matrix(float** matrix){
        // Print the matrix
    printf("The matrix is:\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.4f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void free_matrix(float** matrix){
        // Free allocated memory
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);

    return matrix; 
}

//Compute the Similarity Matrix
float** sym(float** matrix){
    sym_mat;
}

//Compute the diagonal degree Matrix
float** ddg(float** matrix){}

//Compute the normalized similarity matrix
float** norm(float** matrix){}

int main(int argc, char *argv[]) {
    char goal[256];
    char filename[256];
    float** matrix = read_matrix(filename);

    if (strcmp(goal, 'sym') == 0){sym(matrix);}
    if (strcmp(goal, 'ddg') == 0){ddg(matrix);}
    if (strcmp(goal, 'norm') == 0){norm(matrix);}
}