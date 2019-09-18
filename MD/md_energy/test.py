import scipy.optimize as optimize
import numpy as np
import membrane
import matplotlib.pyplot as plt

plt.axis([0, 10, 0, 1])
plt.ion()

y = np.zeros(10)
x = np.linspace(1, 10, 10)

for i in range(100):
    y = np.append(np.delete(y ,0), np.random.random())
    plt.clf()
    plt.plot(x, y)
    plt.pause(0.05)
    
while True:
    plt.pause(0.05)


def fun(x, a, b, m):
    print(m)
    return (x[0]*a[0] - 1)**2 + (x[1] - 2.5)**2*a[1]

aa = np.ones(2)
b = 1
m = membrane.Membrane()

res = optimize.fmin(fun, aa, args=(aa, b, m))
