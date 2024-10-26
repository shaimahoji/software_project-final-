import pandas as pd
import numpy as np
import sys
import mysymnmf as symnmf

np.random.seed(1234)

# Get X (the N datapoints) from the input file.
def read_data(filename):
    data_array =  np.loadtxt(filename, delimiter=',')
    return data_array.tolist() # python list



def similarity_matrix(X):

    n = X.shape[0]

    # Initialization (guarantees that for all i=j: Aij = 0)
    A = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i + 1, n):  # Only calculate for j > i
            # Compute the squared Euclidean distance
            dist_squared = np.sum((X[i] - X[j]) ** 2)
            
            # Apply the similarity formula
            similarity = np.exp(-dist_squared / 2)
            
            # Assign the value symmetrically
            A[i, j] = similarity
            A[j, i] = similarity  # A[i, j] == A[j, i]

    return A


def diagonal_degree_matrix(X):
    A = similarity_matrix(X)
    # Computing the degrees. Example: d1 = the sum of the elements in the first row in A.
    degrees = np.sum(A, axis=1)  # Sum of each row
    D = np.diag(degrees)
    return D


def normalized_similarity_matrix(A, D):
    # Computing D^{-1/2}
    D_inv_sqrt = np.diag(1 / np.sqrt(np.diag(D)))

    # Computing the normalized similarity matrix W. (@ is used for matrix multiplication)
    W = D_inv_sqrt @ A @ D_inv_sqrt
    return W

##### check #####
def deriving_hard_clustering(H):
    # Deriving hard clustering
    return np.argmax(H, axis=1)

# Initialize H : 1.4.1
def initialize_H(W, k):
    # Calculating the average of all entries of W.
    m = np.mean(W)

    # helper variable
    upper_bound = 2 * np.sqrt(m / k)


    # Initializing H with values from the interval [0, upper_bound]
    #H = np.random.uniform(0, upper_bound, (m, k))
    H = np.random.uniform(0, upper_bound, (int(W.shape[0]), int(k)))

    return H

# Update H : 1.4.2
def H_next_t(H_prev_t, W, beta=0.5):

    W_mult_prevH = np.dot(W,H_prev_t)

    # H_prev_t * H_prev_t transposed * H_prev_t
    pH_pHT_pH = np.dot(H_prev_t,H_prev_t.T)
    pH_pHT_pH = np.dot(pH_pHT_pH,H_prev_t)

    # The fraction in the rule
    fract = beta * (W_mult_prevH / pH_pHT_pH)

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
    # Parse command line arguments - K, goal, file_name
    # Given in assignment, we can assume that inputs are valid:

    # int, < N 
    K = int(sys.argv[1])

    # One of the following: symnmf, sym, ddg or norm.
    goal = sys.argv[2]

    #contains N data points for all above goal.
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
