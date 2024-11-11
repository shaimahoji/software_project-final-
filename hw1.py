import math 

# Calculating the Euclidean distance between two vectors
def euclDist(vector1, vector2):

    my_len = len(vector1)
    sum = 0
    
    # Sum of squared differences between corresponding elements
    for i in range(my_len):
        inner_sum = (float(vector1[i]) - float(vector2[i]))**2
        sum += inner_sum
    
    return math.sqrt(sum)

# Initializing centroids as first K datapoints
def initializeCentroids(lines, K):

    init_cent = []

    for i in range(K):

        init_cent.append(lines[i])
    
    return init_cent

# Finding the index of the closest cluster centroid for a given data point
def closestCluster(data_point, centroids):

    # Calculating the distances from the data point to each centroid
    distances = [euclDist(data_point, centroid) for centroid in centroids]

    return distances.index(min(distances)) #index of the minimum distance (closest centroid)


def updateCentroids(data_points, points_clusters, curr_cents, K):
    
    first_line = data_points[0]
    point_len = len(first_line)

    updated = [[0] * point_len for _ in range(K)] # Initializing updated centroids
    counts = [0] * K # Initializing count of the number of points per cluster
    
    # Sum coordinates of points in each cluster
    for point, centroid in zip(data_points, points_clusters):
        for i in range(point_len):
            updated[centroid][i] += float(point[i])
        counts[centroid] += 1

    # Calculating the means (by dividing each sum by the number of points in the cluster)
    for i in range(K):
        if counts[i] != 0:
            updated[i] = [x_i / counts[i] for x_i in updated[i]]
    
    return updated

# Checking for convergence by comparing previous and updated(current) centroids
def HW1Convergence(updated_cents, prev_cents):

    epsilon = 0.001 # The threshold for convergence
    my_len = len(updated_cents)

    for i in range(my_len):

        # If the Euclidean distance between any pair of old and new centroids is greater than epsilon, then: not converged
        if euclDist(updated_cents[i],prev_cents[i]) >= epsilon:
            return False
    
    return True

def Kmeans(lines, K, max_iter = 300):
   
    curr_cents = initializeCentroids(lines, K)
    new_cents = curr_cents
    iter_num = 0

    # Implement do while loop to print list items 
    while(True): 
        curr_cents = new_cents
        points_clusters = []

        # Assigning each datapoint (Xi) to its closest cluster (using its index, between 0 and K-1)
        for line in lines:
            points_clusters.append(closestCluster(line, curr_cents))

        # Updating centroids based on the assigned clusters
        new_cents = updateCentroids(lines, points_clusters, curr_cents, K)
        iter_num += 1

        # Checking for convergence or maximum iteration condition
        if(HW1Convergence(new_cents, curr_cents) == False and iter_num < max_iter): 
            continue
        else: 
            break

    return points_clusters
