import numpy as np
import sys
import symnmf as NMF 
import hw1
from sklearn.metrics import silhouette_score


np.random.seed(1234)

# Calculating the silhouette score for SymNMF clusters (in matrix H)
def silhouetteScoreH(H, data_points):

    # Determining the cluster labels by taking the index of the maximum value in each row of H
    cluster_labels = np.argmax(H, axis=1)

    return silhouette_score(data_points, cluster_labels, metric='euclidean')


def main(): 
    # Parse command line arguments - K, file_name (based on project requirements: we can assume they're valid)
    K = int(sys.argv[1])
    input_file = sys.argv[2]

    # Reading datapoints from the input file
    X = NMF.read_data(input_file)

    # Obtaining matrix H, which contains the SymNMF clusters(they're accessed as explained in 1.5 in the project file)
    H = NMF.handle_symnmf(X,K)

    # Obtaining the kmeans HW1 cluster labels
    hw1_labels = hw1.Kmeans(X, K)
                    
    symnmf_score = silhouetteScoreH(H, X)
    print("nmf: %.4f" % symnmf_score) 

    hw1_score = silhouette_score(X, hw1_labels, metric='euclidean')
    print("kmeans: %.4f" % hw1_score)


if __name__ == "__main__":
    main()
