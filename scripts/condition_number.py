
import numpy as np

def conditionnement(a, b, c):

    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)
    c = np.array(c, dtype=float)

    N = len(b)
    A = np.zeros((N,N))

    for i in range(N):
        A[i,i] = b[i]
        if i > 0:
            A[i, i-1] = a[i-1]
        if i < N-1:
            A[i, i+1] = c[i]

    return np.linalg.cond(A)

a = [-1]*19
b = [2]*20
c = [-1]*19

cond = conditionnement(a, b, c)
print(cond)



