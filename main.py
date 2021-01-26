from geostatspy.GSLIB import varmapv_2d
from pyrcz_kriging import *


def main():
    nug = 0.01
    it1 = "Spherical"
    azi = 1
    hmaj1 = 5
    hmin1 = 5

    x1 = 1
    y1 = 2
    v1 = 2

    x2 = 1
    y2 = 1
    v2 = 2

    x3 = 2
    y3 = 1
    v3 = -1

    xq = 1.99
    yq = 1

    estim, var, weights = f_make_krige(
        nug, it1, azi, hmaj1, hmin1, x1, y1, v1, x2, y2, v2, x3, y3, v3, xq, yq
    )

    # breakpoint()

    print("Estimation:")
    print(estim)
    print()
    print("Variance:")
    print(var)
    print()
    print("Weights:")
    print(weights)


if __name__ == "__main__":
    main()
