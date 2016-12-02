import matplotlib.pyplot as plt

k = 1
Hx = 1
Ht = 0.1
N = 100
Lx = 100
Lt = 100
Nt = int(Lt / Ht)


def init():
    T = []
    for i in range(0, 100):
        if 40 <= i <= 60:
            T.append(1)
        else:
            T.append(0)
    return T


def func2(Nt):
    T = init()
    T1 = []
    for j in range(1, Nt):
        for i in range(1, Lx - 1):
            T1.insert(i, Ht * k * ((T[i + 1] - 2 * T[i] + T[i - 1]) / pow(Hx, 2)) + T[i])
        for i in range(1, Lx - 2):
            T[i] = T1[i]
    return T


x = range(0, 100)

plt.plot(x, func2(0), color=(1, 0, 0))
plt.plot(x, func2(400), color=(0, 1, 0))
plt.plot(x, func2(1000), color=(0, 0, 1))
plt.plot(x, func2(50), color=(0, 0, 0))
plt.show()
