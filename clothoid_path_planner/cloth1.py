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
    dl = 0.01
    a = 2
    k0 = 0
    x = 0
    y = 0
    theta= 0
    x_list = [0]
    y_list = [0]
    for i in range(500):
        r = 2 / (2 * k0 + a * dl)
        k0 += a * dl
        theta += dl / r
        x += dl * cos(theta)
        y += dl * sin(theta)
        x_list.append(x)
        y_list.append(y)
    dl = 0.01
    a = 1
    k0 = 0
    x = 0
    y = 0
    theta= 0
    x_list2 = [0]
    y_list2 = [0]
    for i in range(500):
        r = 2 / (2 * k0 + a * dl)
        k0 += a * dl
        theta += dl / r
        x += dl * cos(theta)
        y += dl * sin(theta)
        x_list2.append(x)
        y_list2.append(y)
    plt.plot(x_list, y_list, '-')
    plt.plot(x_list2, y_list2, '.')
    plt.show()
