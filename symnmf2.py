import pandas as pd
import numpy as np
import sys
import mysymnmf as symnmf

np.random.seed(1234)

# Get X (the N datapoints) from the input file.
def read_data(filename):
    data_array =  np.loadtxt(filename, delimiter=',')
    # Exclude the first column
    filtered_array = data_array[:, 1:]
    
    return filtered_array


def similarity_matrix(X):

    n = X.shape[0]

    # Initialization (guarantees that for all i=j: Aij = 0)
    A = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            if i != j:
                # Computing the squared Euclidean distance.
                dist_squared = np.sum((X[i] - X[j]) ** 2)

                # Applying the similarity formula.
                A[i, j] = np.exp(-dist_squared / 2)
    
    return A


def diagonal_degree_matrix(X):
    # Computing the degrees. Example: d1 = the sum of the elements in the first row in X.
    degrees = np.sum(X, axis=1)  # Sum of each row
    D = np.diag(degrees)
    return D


def normalized_similarity_matrix(A, D):
    # Computing D^{-1/2}
    D_inv_sqrt = np.diag(1 / np.sqrt(np.diag(D)))

    # Computing the normalized similarity matrix W. (@ is used for matrix multiplication)
    W = D_inv_sqrt @ A @ D_inv_sqrt
    return W


def Initialize_H(k,W):
    # Calculating the average of all entries of W.
    m = np.mean(W)

    # helper variable
    upper_bound = 2 * np.sqrt(m / k)

    # Initializing H with values from the interval [0, upper_bound]
    H = np.random.uniform(0, upper_bound, (m, k))


def convergence_condition(H_old, H_new, epsilon):
    """
    Check the Frobenius norm convergence condition:
    ||H(t+1) - H(t)||_F < epsilon
    """
    diff = np.linalg.norm(H_new - H_old, 'fro')
    return diff < epsilon


def optimize_H_until_convergence(W, k, max_iter=300, epsilon=1e-4, beta=0.5):
    """
    Optimize H until either the convergence condition is met or the maximum number of iterations is reached.
    """
    H = initialize_H(W, k)
    for t in range(max_iter):
        H_new = update_H(W, H, beta)
        frobenius_diff = np.linalg.norm(H_new - H, 'fro')
        #print(f"Iteration {t+1}: Frobenius norm difference = {frobenius_diff}")
        
        # Check for convergence
        if convergence_condition(H, H_new, epsilon):
            break

        H = H_new
    
    return H


def deriving_hard_clustering(H):
    # Deriving hard clustering
    return np.argmax(H, axis=1)


def handle_symnmf(X, k):
    #print("You selected symnmf")
    A = handle_sym(X)
    D = handle_ddg(X)
    W = handle_norm(X)
    tmp_H = optimize_H_until_convergence(W,k)
    H = deriving_hard_clustering(tmp_H)
    return H


def handle_sym(X):
    #print("You selected sym")
    A = symnmf.Csym(X)
    #A = similarity_matrix(X)
    return A


def handle_ddg(X):
    #print("You selected ddg")
    D = symnmf.Cddg(X)
    # D = diagonal_degree_matrix(X)
    return D


def handle_norm(X):
    #print("You selected norm")
    A = handle_sym(X)
    D = handle_ddg(X)
    W = normalized_similarity_matrix(A, D)
    return W

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
            return handle_symnmf(X,k)
        else:
            return switch[operation](X)
    else:
        #return "An Error Has Occurred"
        return None


def update_H(W, H, beta=0.5):

    # Initialize H_new as a copy of H
    H_new = np.copy(H)
    
    #The @ operator is used for matrix multiplication. It behaves like np.dot when working with arrays\matrices in NumPy.
    WH = W @ H
    HHT_H = H @ H.T @ H
    
    # Update rule: H_new_ij = H_ij * (0.5 + 0.5 * (WH_ij / HHT_H_ij))
    for i in range(H.shape[0]):
        for j in range(H.shape[1]):
            if HHT_H[i, j] != 0:
                H_new[i, j] = H[i, j] * (1 - beta + beta * (WH[i, j] / HHT_H[i, j]))
    
    return H_new


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

    ########### TRY 1 ##########

    # datapoints
    X = read_data(input_file)
    
    # Number of datapoints
    #N = len(data_array)

    result = switch_case_operation(X, K, goal)

    if result == None:
        print("An Error Has Occurred")
        sys.exit()
    
    print_matrix(result)


if __name__ == "__main__":
    main()