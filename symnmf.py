import pandas as pd
import numpy as np
import sys
import mysymnmf as symnmf

np.random.seed(1234)

# Get X (the N datapoints) from the input file.
def read_data(filename):
    data_array =  np.loadtxt(filename, delimiter=',')
    return data_array.tolist() # python list

# Initialize H : 1.4.1
def initialize_H(W, k):
    # Calculating the average of all entries of W.
    m = np.mean(W)

    # helper variable
    upper_bound = 2 * np.sqrt(m / k)

    # Initializing H with values from the interval [0, upper_bound]
    H = np.random.uniform(0, upper_bound, (int(W.shape[0]), int(k)))

    return H.tolist()



def handle_symnmf(X, k):

    # A and D are found within finding W.
    W = handle_norm(X)
    
    # convert python list to numpy array 
    W = np.array(W)

    init_H = initialize_H(W,k)

    H = symnmf.Csymnmf(init_H,W.tolist())

    return H


def handle_sym(X):
    A = symnmf.Csym(X)
    return A # python list


def handle_ddg(X):
    D = symnmf.Cddg(X)
    return D # python list


def handle_norm(X):
    W =  symnmf.Cnorm(X)
    return W # python list

# Dictionary mapping goals to functions.
switch = {
    "symnmf": handle_symnmf,
    "sym": handle_sym,
    "ddg": handle_ddg,
    "norm": handle_norm
}


def switch_case_operation(X, k, operation):
    if operation in switch:
        if operation == "symnmf":
            return switch[operation](X,k)
        else:
            return switch[operation](X)
    else:
        return None


def print_matrix(matrix):
    for row in matrix:
        formatted_row = ",".join(f"{val:.4f}" for val in row)
        print(formatted_row)


def main(): 
    # Parse command line arguments - K, goal, file_name (based on project file: we can assume they're valid)

    # int, < N 
    K = int(sys.argv[1])

    # One of the following: symnmf, sym, ddg or norm.
    goal = sys.argv[2]

    # Contains N datapoints for all above possible goals.
    input_file = sys.argv[3]

    # datapoints
    X = read_data(input_file)

    result = switch_case_operation(X, K, goal)

    if result is None:
        print("An Error Has Occurred")
        sys.exit()
    
    print_matrix(result)


if __name__ == "__main__":
    main()
