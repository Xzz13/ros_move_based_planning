"""
Cubic spline planner

Author: Atsushi Sakai(@Atsushi_twi)

"""
import math
import numpy as np
import bisect


class Spline:
    """
    Cubic Spline class
    """

    def __init__(self, x, y):
        self.b, self.c, self.d, self.w = [], [], [], []

        self.x = x
        self.y = y

        self.nx = len(x)  # dimension of x
        h = np.diff(x)

        # calc coefficient c
        self.a = [iy for iy in y]

        # calc coefficient c
        A = self.__calc_A(h)
        B = self.__calc_B(h)
        self.c = np.linalg.solve(A, B)
        #  print(self.c1)

        # calc spline coefficient b and d
        for i in range(self.nx - 1):
            self.d.append((self.c[i + 1] - self.c[i]) / (3.0 * h[i]))
            tb = (self.a[i + 1] - self.a[i]) / h[i] - h[i] * \
                (self.c[i + 1] + 2.0 * self.c[i]) / 3.0
            self.b.append(tb)

    def calc(self, t):
        """
        Calc position

        if t is outside of the input x, return None

        """

        if t < self.x[0]:
            return None
        elif t > self.x[-1]:
            return None

        i = self.__search_index(t)
        dx = t - self.x[i]
        result = self.a[i] + self.b[i] * dx + \
            self.c[i] * dx ** 2.0 + self.d[i] * dx ** 3.0

        return result

    def calcd(self, t):
        """
        Calc first derivative

        if t is outside of the input x, return None
        """

        if t < self.x[0]:
            return None
        elif t > self.x[-1]:
            return None

        i = self.__search_index(t)
        dx = t - self.x[i]
        result = self.b[i] + 2.0 * self.c[i] * dx + 3.0 * self.d[i] * dx**2.0
        return result

    def calcdd(self, t):
        """
        Calc second derivative
        """

        if t < self.x[0]:
            return None
        elif t > self.x[-1]:
            return None

        i = self.__search_index(t)
        dx = t - self.x[i]
        result = 2.0 * self.c[i] + 6.0 * self.d[i] * dx
        return result

    def __search_index(self, x):
        """
        search data segment index
        """
        return bisect.bisect(self.x, x) - 1

    def __calc_A(self, h):
        """
        calc matrix A for spline coefficient c
        """
        A = np.zeros((self.nx, self.nx))
        A[0, 0] = 1.0
        for i in range(self.nx - 1):
            if i != (self.nx - 2):
                A[i + 1, i + 1] = 2.0 * (h[i] + h[i + 1])
            A[i + 1, i] = h[i]
            A[i, i + 1] = h[i]

        A[0, 1] = 0.0
        A[self.nx - 1, self.nx - 2] = 0.0
        A[self.nx - 1, self.nx - 1] = 1.0
        #  print(A)
        return A

    def __calc_B(self, h):
        """
        calc matrix B for spline coefficient c
        """
        B = np.zeros(self.nx)
        for i in range(self.nx - 2):
            B[i + 1] = 3.0 * (self.a[i + 2] - self.a[i + 1]) / \
                h[i + 1] - 3.0 * (self.a[i + 1] - self.a[i]) / h[i]
        return B


class Spline2D:
    """
    2D Cubic Spline class

    """

    def __init__(self, x, y):
        self.s = self.__calc_s(x, y)
        self.sx = Spline(self.s, x)
        self.sy = Spline(self.s, y)

    def __calc_s(self, x, y):
        dx = np.diff(x)
        dy = np.diff(y)
        self.ds = np.hypot(dx, dy)
        s = [0]
        s.extend(np.cumsum(self.ds))
        return s

    def calc_position(self, s):
        """
        calc position
        """
        x = self.sx.calc(s)
        y = self.sy.calc(s)

        return x, y

    def calc_curvature(self, s):
        """
        calc curvature
        """
        dx = self.sx.calcd(s)
        ddx = self.sx.calcdd(s)
        dy = self.sy.calcd(s)
        ddy = self.sy.calcdd(s)
        k = (ddy * dx - ddx * dy) / ((dx**2 + dy**2)**(3 / 2))
        return k

    def calc_yaw(self, s):
        """
        calc yaw
        """
        dx = self.sx.calcd(s)
        dy = self.sy.calcd(s)
        yaw = math.atan2(dy, dx)
        return yaw


def calc_spline_course(x, y, ds=0.1):
    sp = Spline2D(x, y)
    s = list(np.arange(0, sp.s[-1], ds))

    rx, ry, ryaw, rk = [], [], [], []
    for i_s in s:
        ix, iy = sp.calc_position(i_s)
        rx.append(ix)
        ry.append(iy)
        ryaw.append(sp.calc_yaw(i_s))
        rk.append(sp.calc_curvature(i_s))

    return rx, ry, ryaw, rk, s


def main():  # pragma: no cover
    print("Spline 2D test")
    import matplotlib.pyplot as plt

    x = [
        -1 + 1.6, -1 + 1.6, -1 + 1.6, -1 + 1.6, 0.25 + 1.6, 1.5 + 1.6,
        1.307 + 1.6, 0.7 + 1.6, 0.4 + 1.6, 0.3 + 1.6, 0.25 + 1.6, 1.15, 0.62,
        0.6
    ]
    y = [
        2.5 + 1.6, 2 + 1.6, 1 + 1.6, 0 + 1.6, -1 + 1.6, 0 + 1.6, 0.707 + 1.6,
        1.1 + 1.6, 1.5 + 1.6, 2 + 1.6, 2.5 + 1.6, 4.7, 4.2, 4.1
    ]
    # x = [
    #     0.291666667, 0.291666667, 0.32, 0.4, 0.48, 0.56, 0.64, 0.72, 0.8, 0.88,
    #     0.96, 1.04, 1.12, 1.2, 1.28, 1.36, 1.44, 1.489583333, 1.489583333,
    #     1.44, 1.36, 1.28, 1.2, 1.12, 1.04, 0.96, 0.88, 0.8, 0.72, 0.64, 0.56,
    #     0.48, 0.4, 0.32, 0.24, 0.16, 0.08, 0, -0.08, -0.16, -0.24, -0.32, -0.4,
    #     -0.48, -0.56, -0.64, -0.72, -0.8, -0.88, -0.96, -1.04, -1.0625,
    #     -1.0625, -1.04, -1.03125, -1.03125
    # ]
    # y = [
    #     1.956360947, 1.806213018, 1.609301775, 1.436458807, 1.326345016,
    #     1.227242604, 1.146059172, 1.084177515, 1.032485207, 0.984147929,
    #     0.940739645, 0.899071006, 0.85135503, 0.796183432, 0.720301775,
    #     0.61691716, 0.463579882, 0.24260355, 0.035502959, -0.141775148,
    #     -0.310023669, -0.43991716, -0.538579882, -0.621088757, -0.688686391,
    #     -0.747213018, -0.799236686, -0.844467456, -0.876278107, -0.908088757,
    #     -0.935508876, -0.956384615, -0.968313609, -0.976473373, -0.979289941,
    #     -0.975189349, -0.963260355, -0.951331361, -0.927639053, -0.898065089,
    #     -0.865177515, -0.825414201, -0.779023669, -0.724970414, -0.662591716,
    #     -0.591017751, -0.508384615, -0.409970414, -0.291053254, -0.146745562,
    #     0.101360947, 0.294378698, 0.522189349, 0.71852071, 0.770710059,
    #     0.905325444
    # ]

    ds = 0.01  # [m] distance of each intepolated points

    sp = Spline2D(x, y)
    s = np.arange(0, sp.s[-1], ds)

    rx, ry, ryaw, rk = [], [], [], []
    for i_s in s:
        ix, iy = sp.calc_position(i_s)
        rx.append(ix)
        ry.append(iy)
        ryaw.append(sp.calc_yaw(i_s))
        rk.append(sp.calc_curvature(i_s))

    plt.subplots(1)
    plt.plot(x, y, "xb", label="input")
    plt.plot(rx, ry, "-r", label="spline")
    plt.grid(True)
    plt.axis("equal")
    plt.xlabel("x[m]")
    plt.ylabel("y[m]")
    plt.legend()

    plt.subplots(1)
    plt.plot(s, [np.rad2deg(iyaw) for iyaw in ryaw], "-g", label="yaw")
    plt.grid(True)
    plt.legend()
    plt.xlabel("line length[m]")
    plt.ylabel("yaw angle[deg]")

    plt.subplots(1)
    plt.plot(s, rk, "-b", label="curvature")
    plt.grid(True)
    plt.legend()
    plt.xlabel("line length[m]")
    plt.ylabel("curvature [1/m]")

    plt.show()


if __name__ == '__main__':
    main()
