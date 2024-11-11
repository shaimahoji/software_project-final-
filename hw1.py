import math


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