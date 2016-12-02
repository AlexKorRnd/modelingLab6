#!/usr/bin/python
import matplotlib.pyplot as plt
import math

k = 1
Hx = 1
Ht = 0.1
N = 101
Lx = 100
Lt = 100
Nt = int(Lt / Ht)
sigma = 0.5


def init():
    T = []
    for i in range(0, N):
        if 40 <= i <= 60:
            T.append(1)
        else:
            T.append(0)
    return T


def func3(count):
    T = init()
    for j in range(0, count):
        a = []
        b = []
        c = []
        f = []
        for i in range(0, N):
            a.insert(i, k / pow(Hx, 2))
            b.insert(i, k / pow(Hx, 2))
            c.insert(i, (1 / Ht) + (2 * k / pow(Hx, 2)))
            f.insert(i, T[i] / Ht)
        alpha = [0]
        beta = [0]
        for i in range(0, N - 1):
            alpha.insert(i+i, b[i] / (c[i] - a[i] * alpha[i]))
            beta.insert(i+i, (f[i] + a[i] * beta[i]) / (c[i] - a[i] * alpha[i]))
        i = N - 1
        T[N - 1] = (f[i] + a[i] * beta[i]) / (c[i] - a[i] * alpha[i])
        for i in range(N - 2, 0, -1):
            T[i] = beta[i + 1] + alpha[i + 1] * T[i + 1]
    return T


def func4(count):
    T = init()
    for j in range(0, count):
        a = [0]
        b = [0]
        c = [1]
        f = [math.sin(math.pi * j * Ht / 100)]
        for i in range(1, N - 1):
            a.insert(i, k / pow(Hx, 2) * sigma)
            b.insert(i, k / pow(Hx, 2) * sigma)
            c.insert(i, 1 / Ht + 2 * k / pow(Hx, 2) * sigma)
            f.insert(i, T[i] / Ht + (T[i + 1] - 2 * T[i] + T[i - 1]) / pow(Hx, 2) * (1 - sigma))
        i = N - 1
        alpha = [0]
        beta = [0]
        c.insert(i, 1 / Ht + 0.1 / Hx + k / pow(Hx, 2))
        b.insert(i, k / pow(Hx, 2))
        a.insert(i, k / pow(Hx, 2))
        f.insert(i, T[i] / Ht + k * (T[i] - T[i - 1]) / pow(Hx, 2) + k * (0.2 - 0.1 * T[i]) / pow(Hx, 2))
        for i in range(0, N - 1):
            alpha.insert(i + 1, b[i] / (c[i] - a[i] * alpha[i]))
            beta.insert(i + 1, (f[i] + a[i] * beta[i]) / (c[i] - a[i] * alpha[i]))
        i = N - 1
        T[i] = (f[i] + a[i] * beta[i]) / (c[i] - a[i] * alpha[i])
        #T[N-1] =
        for i in range(N - 2, 0, -1):
            T[i] = beta[i + 1] + alpha[i + 1] * T[i + 1]
    return T


x = range(0, N)
y = func4(500)
plt.plot(x, y, color=(0, 0, 0))
plt.show()
