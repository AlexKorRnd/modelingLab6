#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import math

Nx = 100
Ny = 100

N = Nx * Ny

Hx = 1
Hy = 1
Ht = 0.1
Sigma = 0.5

count = 10000

Lt = 10

Cn = []


def lab5():
    O = [1 for i in range(0, count)]
    B1 = [0 for i in range(0, count)]
    B2 = [0 for i in range(0, count)]
    B3 = [0 for i in range(0, count)]
    B4 = [0 for i in range(0, count)]
    B5 = [0 for i in range(0, count)]
    A = [0 for i in range(0, count)]
    Ux = [3 for i in range(0, count)]
    Uy = [4 for i in range(0, count)]
    mu = [2 for i in range(0, count)]

    F = []

    C = []

    O = []
    u = []
    v = []
    mu = []

    F = [0 for i in range(0, 10000)]

    for k in range(0, N):
        O.insert(k, 1)
        u.insert(k, 3)
        v.insert(k, 4)
        mu.insert(k, 2)

    for i in range(0, Nx):
        line = []
        for j in range(0, Ny):
            line.append(0)
        Cn.append(line)

    for x in range(10, 20):
        for y in range(10, 20):
            Cn[x][y] = math.sin(math.pi * (x - 10) / 10) * math.sin(math.pi * (y - 10) / 10)

    l = 0
    for i in range(0, Nx):
        for j in range(0, Ny):
            C.insert(l, Cn[i][j])
            l += 1

    t = 0
    while True:
        for i in range(1, Nx - 2):
            for j in range(1, Ny - 2):
                m0 = i + j * Nx
                m1 = m0 + 1
                m2 = m0 - 1
                m3 = m0 + Nx
                m4 = m0 - Nx
                m24 = m0 - 1 - Nx

                q1 = (O[m0] + O[m4]) / 2
                q2 = (O[m2] + O[m24]) / 2
                q3 = (O[m0] + O[m2]) / 2
                q4 = (O[m4] + O[m24]) / 2
                q0 = (q1 + q2) / 2

                B1[m0] = q1 * (-(u[m1] + u[m0]) / (4 * Hx) + (mu[m1] + mu[m0]) / (2 * Hx * Hx))
                B2[m0] = q2 * ((u[m2] + u[m0]) / (4 * Hx) + (mu[m2] + mu[m0]) / (2 * Hx * Hx))
                B3[m0] = q3 * (-(v[m3] + v[m0]) / (4 * Hy) + (mu[m3] + mu[m0]) / (2 * Hy * Hy))
                B4[m0] = q4 * ((v[m4] + v[m0]) / (4 * Hy) + (mu[m4] + mu[m0]) / (2 * Hy * Hy))

                B6 = (1 - Sigma) * B1[m0]
                B7 = (1 - Sigma) * B2[m0]
                B8 = (1 - Sigma) * B3[m0]
                B9 = (1 - Sigma) * B4[m0]

                B1[m0] *= Sigma
                B2[m0] *= Sigma
                B3[m0] *= Sigma
                B4[m0] *= Sigma

                A[m0] = q0 / Ht + B1[m0] + B2[m0] + B3[m0] + B4[m0]
                B5[m0] = q0 / Ht - B6 - B7 - B8 - B9

                F[m0] = B5[m0] * C[m0] + B6 * C[m1] + B7 * C[m2] + B8 * C[m3] + B9 * C[m4]

        C = seidel(B1, B2, B3, B4, A, F, C)
        t += Ht
        if t >= Lt:
            break
    return C


pogr = 1 / math.pow(10, 8)


def seidel(B1, B2, B3, B4, A, F, C):
    maxPogr = 2 * pogr
    # Cres = copyArray(C)
    while True:
        maxPogr = 0
        for i in range(1, Nx - 2):
            for j in range(1, Ny - 2):
                m0 = i + j * Nx
                m1 = m0 + 1
                m2 = m0 - 1
                m3 = m0 + Nx
                m4 = m0 - Nx
                w = C[m0]
                C[m0] = (F[m0] + B1[m0] * C[m1] + B2[m0] * C[m2] + B3[m0] * C[m3] + B4[m0] * C[m4]) / A[m0]
                w = math.fabs(w - C[m0])
                if w > maxPogr:
                    maxPogr = w
        if maxPogr < pogr:
            break

    return C


def copyArray(array):
    copy = []
    for i in range(0, len(array)):
        copy.append(array[i])
    return copy


def convertArrayToMatrix(c):
    Matrix = []
    index = 0
    for i in range(0, Nx):
        line = []
        for j in range(0, Ny):
            line.append(c[index])
            index += 1
        Matrix.append(line)
    return Matrix


C = convertArrayToMatrix(lab5())

data = np.loadtxt('data.txt')

f1 = open('./testfile', 'w+')
for i in range(0, len(C)):
    for j in range(0, len(C[i])):
        f1.write(str(C[i][j]))
        f1.write("  ")
    f1.write("\n")

plt.contour(C)
plt.show()
