import pandas as pd
import numpy as np
import sys
import math
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

np.random.seed(1234)

# Get X (the N datapoints) from the input file.
def read_data(filename):
    data_array =  np.loadtxt(filename, delimiter=',')
    return data_array

##### SymNMF section (APPARENTLY ALSO NEED THE ENTIRE FILE, AND IT NEEDS TO BE WIH C MODULE, FOR NOW I USE MY PYTHON VERSION) #####

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
    A = similarity_matrix(X)
    D = diagonal_degree_matrix(X)
    W = normalized_similarity_matrix(A, D)
    H = converge(W, k)

    return H
########################## END of SymNMF section ##########################


########################## START of HW1 section ##########################
def euclDist(vector1, vector2):

    my_len = len(vector1)
    sum = 0

    for i in range(my_len):
        inner_sum = (float(vector1[i]) - float(vector2[i]))**2
        sum += inner_sum
    
    return math.sqrt(sum)


def initializeCentroids(lines, K):

    init_cent = []

    #Initialize centroids as first K datapoints
    for i in range(K):

        init_cent.append(lines[i])
    
    return init_cent


def closestCluster(data_point, centroids):

    distances = [euclDist(data_point, centroid) for centroid in centroids]
    return distances.index(min(distances))


def updateCentroids(data_points, points_clusters, curr_cents, K):
    
    first_line = data_points[0]
    point_len = len(first_line)

    updated = [[0] * point_len for _ in range(K)]
    counts = [0] * K
    
    for point, centroid in zip(data_points, points_clusters):
        for i in range(point_len):
            updated[centroid][i] += float(point[i])
        counts[centroid] += 1
    
    for i in range(K):
        if counts[i] != 0:
            updated[i] = [x_i / counts[i] for x_i in updated[i]]
    
    return updated


def HW1Convergence(updated_cents, prev_cents):
    epsilon = 0.001
    my_len = len(updated_cents)
    for i in range(my_len):
        if euclDist(updated_cents[i],prev_cents[i]) >= epsilon:
            return False
    
    return True

# Main func
def Kmeans(lines, K, max_iter = 300):
   
    curr_cents = initializeCentroids(lines, K)
    new_cents = curr_cents
    iter_num = 0

    # Implement do while loop to print list items 
    while(True): 
        curr_cents = new_cents
        points_clusters = []

        # Assign every Xi to its closest cluster
        for line in lines:
            # for each Xi : the element i in this list represents the index (between 0 and K-1) of the cluster that Xi belongs to.
            points_clusters.append(closestCluster(line, curr_cents))

        new_cents = updateCentroids(lines, points_clusters, curr_cents, K)
        iter_num += 1
        if(HW1Convergence(new_cents, curr_cents) == False and iter_num < max_iter): 
            continue
        else: 
            break

    return points_clusters
#### end of HW1 section #####


def silhouetteScoreH(H, data_points):
    cluster_labels = np.argmax(H, axis=1)
    #return silhouette_score(H, cluster_labels)
    return silhouette_score(data_points, cluster_labels, metric='euclidean')

# MAIN
def main(): 
    # Parse command line arguments - K, goal, file_name
    # Given in assignment, we can assume that inputs are valid:

    # int, < N 
    K = int(sys.argv[1])

    #contains N data points for all above goal.
    input_file = sys.argv[2]

    # datapoints
    X = read_data(input_file)

    #the SymNMF clusters are in here, they're accessed as explained in 1.5 in the PDF
    H = handle_symnmf(X,K) #according to previous tests: CORRECT

    #the kmeans HW1 clusters
    hw1_labels = Kmeans(X, K)
                    
    symnmf_score = silhouetteScoreH(H, X)
    print("nmf: %.4f" % symnmf_score) 

    hw1_score = silhouette_score(X, hw1_labels, metric='euclidean')
    print("kmeans: %.4f" % hw1_score)

if __name__ == "__main__":
    main()