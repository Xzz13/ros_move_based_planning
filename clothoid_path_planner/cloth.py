# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot as plt
from sympy import *
import math

fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax1.set_title('Fig')
plt.xlabel('x(m)')
plt.ylabel('y(m)')

if __name__ == '__main__':
    pi = 3.14159265
    a = 0.5
    x_list = [0]
    y_list = [0]
    t = symbols('t')
    for i in range(500):
        theta = (i + 1) * 2 * 3.14159265 / 360
        x = a * integrate(cos(t * t), (t, 0, theta))
        y = a * integrate(sin(t * t), (t, 0, theta))
        x_list.append(x)
        y_list.append(y)
        print(i, x, y)
    plt.plot(x_list, y_list, '.')
    plt.show()
