import pandas as pd
import numpy as np
import sys
import mysymnmf as symnmf

np.random.seed(1234)

# Get X (the N datapoints) from the input file
def read_data(filename):
    data_array =  np.loadtxt(filename, delimiter=',')
    return data_array.tolist() # python list

# Initialize H : 1.4.1
def initialize_H(W, k):

    # Calculating the average of all entries of W
    m = np.mean(W)

    # The upper bound for values in H
    upper_bound = 2 * np.sqrt(m / k)

    # Initializing H with random values from the interval [0, upper_bound]
    H = np.random.uniform(0, upper_bound, (int(W.shape[0]), int(k)))

    return H.tolist() # Return H as a Python list



def handle_symnmf(X, k):

    # Find the normalized version of X to obtain W
    W = handle_norm(X)
    
    # Convert W to a numpy array for further processing
    W = np.array(W)

    # Initialize the H matrix
    init_H = initialize_H(W,k)

    # Call C-implemented SymNMF function
    H = symnmf.Csymnmf(init_H,W.tolist())

    return H

# Call C-implemented function that computes the similarity matrix from X
def handle_sym(X):
    A = symnmf.Csym(X)
    return A # python list

# Call C-implemented function that computes the Diagonal Degree Matrix
def handle_ddg(X):
    D = symnmf.Cddg(X)
    return D # python list

# Call C-implemented function that computes the normalized similarity Matrix
def handle_norm(X):
    W =  symnmf.Cnorm(X)
    return W # python list

# Dictionary to map goal strings to corresponding function handlers
switch = {
    "symnmf": handle_symnmf,
    "sym": handle_sym,
    "ddg": handle_ddg,
    "norm": handle_norm
}

# Calling the appropriate function based on the provided goal
def switch_case_operation(X, k, operation):
    if operation in switch:
        if operation == "symnmf":
            return switch[operation](X,k)
        else:
            return switch[operation](X)
    else:
        return None

# Helper function to print matrix values with 4 decimal places
def print_matrix(matrix):
    for row in matrix:
        formatted_row = ",".join(f"{val:.4f}" for val in row)
        print(formatted_row)


def main(): 
    # Parse command line arguments - K, goal, file_name (based on project file: we can assume they're valid)
    K = int(sys.argv[1])
    goal = sys.argv[2]
    input_file = sys.argv[3]

    # Load data points from the input file
    X = read_data(input_file)

    # Executing the operation based on the goal
    result = switch_case_operation(X, K, goal)

    # Checking if the operation was successful
    if result is None:
        print("An Error Has Occurred")
        sys.exit()
    
    print_matrix(result)


if __name__ == "__main__":
    main()
