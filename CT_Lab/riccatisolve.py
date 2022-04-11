import scipy.linalg 
import numpy as np

data = []
with open("ABQR.txt") as f:
    for line in f:
        data.append([float(x) for x in line.split()])

A = np.array(data[0:4])
B = np.array(data[5:9])
Q = np.array(data[10:14])
R = np.array(data[15:17])

'''
print(A)
print(B)
print(Q)
print(R)
'''

X = scipy.linalg.solve_continuous_are(A, B, Q, R)
np.savetxt("X.txt", X)
