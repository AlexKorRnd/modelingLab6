#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from numpy import array

import math
from operator import lt

Nx = 20
Ny = 20

N = Nx * Ny

Hx = 1
Hy = 1
Ht = 1
Sigma = 0.5
ro = 0.9982
pogr = 1 / math.pow(10, 8)

count = N

Lt = 10

r = 1000

Cn = []

u = []
v = []
O = []
U = []
V = []
mu = []
B1 = []
B2 = []
B3 = []
B4 = []
B5 = []
B6 = []
B7 = []
B7 = []
B8 = []
B9 = []
A = []
F = []
Fu = []
Fv = []
P = []


def initArrays():
    for i in range(0, N):
        O.append(0.0)
        u.append(0.0)
        v.append(0.0)
        U.append(0.0)
        V.append(0.0)
        mu.append(0.0)
        B1.append(0.0)
        B2.append(0.0)
        B3.append(0.0)
        B4.append(0.0)
        B5.append(0.0)
        B6.append(0.0)
        B7.append(0.0)
        B8.append(0.0)
        B9.append(0.0)
        A.append(0.0)
        F.append(0.0)
        Fu.append(0.0)
        Fv.append(0.0)
        P.append(0.0)


def replaceArray(dest, source):
    if len(dest) < len(source):
        for i in range(len(dest), len(source)):
            dest.append(0)

    for i in range(0, len(source)):
        dest[i] = source[i]


def seidel(np, vp, A, B1, B2, B3, B4, UV, F, pogr):
    maxPogr = 2 * pogr
    while True:
        maxPogr = 0
        for i in range(np, vp):
            for j in range(np, vp):
                m0 = i + j * Nx
                m1 = m0 + 1
                m2 = m0 - 1
                m3 = m0 + Nx
                m4 = m0 - Nx
                w = UV[m0]
                UV[m0] = (F[m0] + B1[m0] * UV[m1] + B2[m0] * UV[m2] + B3[m0] * UV[m3] + B4[m0] * UV[m4]) / float(A[m0])
                w = math.fabs(w - UV[m0])
                if w > maxPogr:
                    maxPogr = w
        if maxPogr < pogr:
            break

    return UV


def func1():
    #print "func1"
    for i in range(2, Nx - 2):
        for j in range(2, Ny - 2):
            m0 = i + j * Nx
            m1 = m0 + 1
            m2 = m0 - 1
            m3 = m0 + Nx
            m4 = m0 - Nx
            m24 = m0 - 1 - Nx

            q1 = float(O[m0] + O[m4]) / 2
            q2 = float(O[m2] + O[m24]) / 2
            q3 = float(O[m0] + O[m2]) / 2
            q4 = float(O[m4] + O[m24]) / 2
            q0 = float(q1 + q2) / 2

            B1[m0] = q1 * (-(U[m1] + U[m0]) / float(4 * Hx) + (mu[m1] + mu[m0]) / (2 * Hx * Hx))
            B2[m0] = q2 * ((U[m2] + U[m0]) / float(4 * Hx) + (mu[m2] + mu[m0]) / (2 * Hx * Hx))
            B3[m0] = q3 * (-(V[m3] + V[m0]) / float(4 * Hy) + (mu[m3] + mu[m0]) / (2 * Hy * Hy))
            B4[m0] = q4 * ((V[m4] + V[m0]) / float(4 * Hy) + (mu[m4] + mu[m0]) / (2 * Hy * Hy))

            B6[m0] = (1 - Sigma) * B1[m0]
            B7[m0] = (1 - Sigma) * B2[m0]
            B8[m0] = (1 - Sigma) * B3[m0]
            B9[m0] = (1 - Sigma) * B4[m0]

            B1[m0] *= Sigma
            B2[m0] *= Sigma
            B3[m0] *= Sigma
            B4[m0] *= Sigma

            A[m0] = float(q0) / Ht + B1[m0] + B2[m0] + B3[m0] + B4[m0]
            B5[m0] = float(q0) / Ht - B6[m0] - B7[m0] - B8[m0] - B9[m0]

            Fu[m0] = B5[m0] * U[m0] + B6[m0] * U[m1] + B7[m0] * U[m2] + B8[m0] * U[m3] + B9[m0] * U[m4]
            Fv[m0] = B5[m0] * V[m0] + B6[m0] * V[m1] + B7[m0] * V[m2] + B8[m0] * V[m3] + B9[m0] * V[m4]

    replaceArray(U, seidel(2, Nx - 2, A, B1, B2, B3, B4, U, Fu, pogr))
    replaceArray(V, seidel(2, Nx - 2, A, B1, B2, B3, B4, V, Fv, pogr))


def func2():
    #print "func2"
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            m0 = i + j * Nx
            m1 = m0 + 1
            m2 = m0 - 1
            m3 = m0 + Nx
            m4 = m0 - Nx
            m24 = m0 - 1 - Nx

            q1 = float(O[m0] + O[m4]) / 2
            q2 = float(O[m2] + O[m24]) / 2
            q3 = float(O[m0] + O[m2]) / 2
            q4 = float(O[m4] + O[m24]) / 2
            q0 = float(q1 + q2) / 2

            B1[m0] = q1 / float(Hx * Hx)
            B2[m0] = q2 / float(Hx * Hx)
            B3[m0] = q3 / float(Hy * Hy)
            B4[m0] = q4 / float(Hy * Hy)

            A[m0] = B1[m0] + B2[m0] + B3[m0] + B4[m0]
            F[m0] = (-ro / float(Ht)) * ((q1 * (U[m1] + U[m0]) -
                 q2 * (U[m2] + U[m0])) / float(2 * Hx) + (q3 * (V[m3] + V[m0])
                 - q4 * (V[m4] + V[m0])) / float(2 * Hy))
    replaceArray(P, seidel(1, Nx - 1, A, B1, B2, B3, B4, P, F, pogr))


def func3():
    #print "func3"
    for i in range(2, Nx - 2):
        for j in range(2, Ny - 2):
            m0 = i + j * Nx
            m1 = m0 + 1
            m2 = m0 - 1
            m3 = m0 + Nx
            m4 = m0 - Nx
            m24 = m0 - 1 - Nx

            q1 = float(O[m0] + O[m4]) / 2
            q2 = float(O[m2] + O[m24]) / 2
            q3 = float(O[m0] + O[m2]) / 2
            q4 = float(O[m4] + O[m24]) / 2
            q0 = float(q1 + q2) / 2

            U[m0] -= (Ht / (ro * q0)) * (q1 * (P[m1] - P[m0]) / float(2 * Hx) + q2 * (P[m0] - P[m2]) / float(2 * Hx))
            V[m0] -= (Ht / (ro * q0)) * (q3 * (P[m3] - P[m0]) / float(2 * Hy) + q4 * (P[m0] - P[m4]) / float(2 * Hy))


def lab6():
    initArrays()
    for k in range(0, N):
        u[k] = 3
        v[k] = 4
        mu[k] = 2
        O[k] = 0
        U[k] = 0

    j = 1
    for i in range(1, Nx - 2):
        U[i + j * Nx] = 1

    for i in range(1, Nx - 2):
        for j in range(1, Ny - 2):
            O[i + j * Nx] = 1

    for k in range(0, N):
        u[k] = 3
        v[k] = 4
        mu[k] = 2

    t = 0
    while True:
        func1()
        func2()
        func3()
        t += Ht
        print "t = ", t
        if t > Lt:
            break


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
            line.append(c[index] * 100000)
            index += 1
        Matrix.append(line)
    return Matrix


def printMatrix(matrix):
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix[i])):
            print "matrix[", i, "][", j, "]", matrix[i][j]



lab6()
C_U = convertArrayToMatrix(U)
C_V = convertArrayToMatrix(V)
C_P = convertArrayToMatrix(P)

Y, X = np.mgrid[1:Nx:20j, 1:Ny:20j]
U = -1 - X**2 + Y

U1 = array(C_U)
V1 = array(C_V)

# fig0, ax0 = plt.subplots()
# strm = ax0.streamplot(X, Y, U1, V1, color=U, linewidth=2, cmap=plt.cm.autumn)
#
#
# plt.show()



plt.figure()
Q = plt.quiver(U1, V1)
qk = plt.quiverkey(Q, 0, 0, 0, '')
l, r, b, t = plt.axis()
dx, dy = r - l, t - b
plt.axis([l, r, b, t])

plt.show()
