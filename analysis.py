import pandas as pd
import numpy as np
import sys
import symnmf as NMF 
import hw1
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


np.random.seed(1234)

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
    X = NMF.read_data(input_file)
    #the SymNMF clusters are in here, they're accessed as explained in 1.5 in the PDF
    H = NMF.handle_symnmf(X,K) #according to previous tests: CORRECT

    #the kmeans HW1 clusters
    hw1_labels = hw1.Kmeans(X, K)
                    
    symnmf_score = silhouetteScoreH(H, X)
    print("nmf: %.4f" % symnmf_score) 

    hw1_score = silhouette_score(X, hw1_labels, metric='euclidean')
    print("kmeans: %.4f" % hw1_score)

if __name__ == "__main__":
    main()