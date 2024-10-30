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

    return H

# Update H : 1.4.2
def H_next_t(H_prev_t, W, beta=0.5):
    # The Numerator (top part) of the fraction in the rule (in the project file)
    W_mult_prevH = np.dot(W,H_prev_t)

    # Previous H (Hi,j to the power of t) multiplied (matrix multiplication) by its transpose
    pH_pHT_pH = np.dot(H_prev_t,H_prev_t.T)

    # The Denominator (bottom part) of the fraction
    pH_pHT_pH = np.dot(pH_pHT_pH,H_prev_t)

    # The fraction itself
    fract = beta * (W_mult_prevH / pH_pHT_pH)

    # The desrired result (H sub i,j to the power of t+1)
    next_H = H_prev_t * ((1-beta) + fract)

    return next_H

# Convergence : 1.4.3
def converge(W, k, max_iter=300, epsilon=1e-4, beta=0.5):
    init_H = initialize_H(W,k)
    prev_H = init_H

    for i in range(max_iter):

        current_H = H_next_t(prev_H, W, beta)

        # Calculating squared ||.||F
        frobe_squared = np.sum((current_H - prev_H)**2)

        if i == max_iter or frobe_squared < epsilon:
            return current_H
        
        prev_H = current_H


def handle_symnmf(X, k):

    # A and D are found within finding W.
    W = handle_norm(X)
    
    # convert python list to numpy array 
    W = np.array(W)

    H = converge(W, k)

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
