import sys
from functools import lru_cache
import numpy as np

# Function to calculate Levenshtein (edit) distance between two strings
def levenshtein(s1, s2):
    len_s1, len_s2 = len(s1), len(s2)
    dp = np.zeros((len_s1 + 1, len_s2 + 1), dtype=int)

    for i in range(len_s1 + 1):
        dp[i][0] = i
    for j in range(len_s2 + 1):
        dp[0][j] = j

    for i in range(1, len_s1 + 1):
        for j in range(1, len_s2 + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            dp[i][j] = min(dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + cost)

    return dp[len_s1][len_s2]

# Function to convert strings into a distance matrix (based on Levenshtein distance)
def create_distance_matrix(cities):
    n = len(cities)
    dist = np.zeros((n, n), dtype=int)

    for i in range(n):
        for j in range(i + 1, n):
            distance = levenshtein(cities[i], cities[j])
            dist[i][j] = dist[j][i] = distance

    return dist

# Dynamic Programming approach to solve TSP using Bitmasking
def tsp(dist, n):
    # dp[mask][i] represents the minimum cost to visit all cities in mask and end at city i
    dp = [[sys.maxsize] * n for _ in range(1 << n)]
    dp[1][0] = 0  # Start from city 0 with only city 0 visited

    # Fill the dp table
    for mask in range(1 << n):
        for u in range(n):
            if (mask & (1 << u)) == 0:  # City u not in the mask
                continue
            # Try to extend the tour from city u to city v
            for v in range(n):
                if (mask & (1 << v)) == 0:  # City v not visited
                    new_mask = mask | (1 << v)
                    dp[new_mask][v] = min(dp[new_mask][v], dp[mask][u] + dist[u][v])

    # Reconstruct the Hamiltonian Path
    min_cost = sys.maxsize
    last_city = -1
    for i in range(1, n):
        if dp[(1 << n) - 1][i] + dist[i][0] < min_cost:
            min_cost = dp[(1 << n) - 1][i] + dist[i][0]
            last_city = i

    # Reconstruct the path
    path = []
    mask = (1 << n) - 1
    while last_city != -1:
        path.append(last_city)
        next_city = -1
        for i in range(n):
            if (mask & (1 << i)) != 0 and dp[mask][last_city] == dp[mask ^ (1 << last_city)][i] + dist[i][last_city]:
                next_city = i
                break
        mask ^= (1 << last_city)
        last_city = next_city

    path.reverse()
    return path, min_cost

# Example usage
cities = ["dog", "cat", "bat", "rat", "mat"]

# Step 1: Create the distance matrix based on Levenshtein distance
dist_matrix = create_distance_matrix(cities)

# Step 2: Solve TSP and get the Hamiltonian Path
n = len(cities)
path, min_cost = tsp(dist_matrix, n)

# Output the result
print("Hamiltonian Path:", [cities[i] for i in path])
print("Minimum cost (Levenshtein distance):", min_cost)
