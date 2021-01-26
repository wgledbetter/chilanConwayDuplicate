import geostatspy.GSLIB as GSLIB  # GSLIB utilies, visualization and wrapper
import geostatspy.geostats as geostats  # GSLIB methods convert to Python

import os  # to set current working directory
import sys  # supress output to screen for interactive variogram modeling
import io
import numpy as np  # arrays and matrix math
import pandas as pd  # DataFrames
import matplotlib.pyplot as plt  # plotting
from matplotlib.pyplot import cm  # color maps
from matplotlib.patches import Ellipse  # plot an ellipse
import math  # sqrt operator
from scipy.stats import norm


def simple_simple_krige(df, xcol, ycol, vcol, dfl, xlcol, ylcol, vario, skmean):
    # load the variogram
    nst = vario["nst"]
    pmx = 9999.9
    cc = np.zeros(nst)
    aa = np.zeros(nst)
    it = np.zeros(nst)
    ang = np.zeros(nst)
    anis = np.zeros(nst)
    nug = vario["nug"]
    sill = nug
    cc[0] = vario["cc1"]
    sill = sill + cc[0]
    it[0] = vario["it1"]
    ang[0] = vario["azi1"]
    aa[0] = vario["hmaj1"]
    anis[0] = vario["hmin1"] / vario["hmaj1"]
    if nst == 2:
        cc[1] = vario["cc2"]
        sill = sill + cc[1]
        it[1] = vario["it2"]
        ang[1] = vario["azi2"]
        aa[1] = vario["hmaj2"]
        anis[1] = vario["hmin2"] / vario["hmaj2"]

    # set up the required matrices
    rotmat, maxcov = geostats.setup_rotmat(nug, nst, it, cc, ang, pmx)
    ndata = len(df)
    a = np.zeros([ndata, ndata])
    r = np.zeros(ndata)
    s = np.zeros(ndata)
    rr = np.zeros(ndata)
    nest = len(dfl)

    est = np.zeros(nest)
    var = np.full(nest, sill)
    weights = np.zeros([nest, ndata])

    # Make and solve the kriging matrix, calculate the kriging estimate and variance
    for iest in range(0, nest):
        for idata in range(0, ndata):
            for jdata in range(0, ndata):
                a[idata, jdata] = geostats.cova2(
                    df[xcol].values[idata],
                    df[ycol].values[idata],
                    df[xcol].values[jdata],
                    df[ycol].values[jdata],
                    nst,
                    nug,
                    pmx,
                    cc,
                    aa,
                    it,
                    ang,
                    anis,
                    rotmat,
                    maxcov,
                )
            r[idata] = geostats.cova2(
                df[xcol].values[idata],
                df[ycol].values[idata],
                dfl[xlcol].values[iest],
                dfl[ylcol].values[iest],
                nst,
                nug,
                pmx,
                cc,
                aa,
                it,
                ang,
                anis,
                rotmat,
                maxcov,
            )
            rr[idata] = r[idata]

        s = geostats.ksol_numpy(ndata, a, r)
        sumw = 0.0
        for idata in range(0, ndata):
            sumw = sumw + s[idata]
            weights[iest, idata] = s[idata]
            est[iest] = est[iest] + s[idata] * df[vcol].values[idata]
            var[iest] = var[iest] - s[idata] * rr[idata]
        est[iest] = est[iest] + (1.0 - sumw) * skmean
    return est, var, weights


def convert_type(it):
    if it == "Spherical":
        return 1
    elif it == "Exponential":
        return 2
    else:
        return 3


def f_make_krige(
    nug, it1, azi, hmaj1, hmin1, x1, y1, v1, x2, y2, v2, x3, y3, v3, xq, yq
):  # function to take parameters, make sample and plot

    it1 = convert_type(it1)
    nst = 1
    xlag = 10
    nlag = int(hmaj1 / xlag)
    c1 = 1.0 - nug
    vario = GSLIB.make_variogram(
        nug, nst, it1, c1, azi, hmaj1, hmin1
    )  # make model object
    # index_maj, h_maj, gam_maj, cov_maj, ro_maj = geostats.vmodel(
    #     nlag, xlag, azi, vario
    # )  # project the model in the major azimuth
    # index_min, h_min, gam_min, cov_min, ro_min = geostats.vmodel(
    #     nlag, xlag, azi + 90.0, vario
    # )  # project the model in the minor azimuth

    x = [x1, x2, x3]
    y = [y1, y2, y3]
    value = [v1, v2, v3]
    skm = np.mean(value)
    df = pd.DataFrame({"X": x, "Y": y, "Value": value})

    xl = [xq, 0, 1]
    yl = [yq, 0, 1]
    value1 = [0, 0, 0]
    dfl = pd.DataFrame({"X": xl, "Y": yl, "Value": value1})

    sk_est, sk_var, sk_weights = simple_simple_krige(
        df, "X", "Y", "Value", dfl, "X", "Y", vario, skmean=skm
    )
    if sk_var[0] == 0:
        sk_std = 0.0
    else:
        sk_std = math.sqrt(sk_var[0])

    return sk_est, sk_var, sk_weights

    # xlag = 10.0
    # nlag = int(hmaj1 / xlag)

    # plt.subplot(1, 3, 1)
    # plt.plot([0, hmaj1 * 1.5], [1.0, 1.0], color="black")
    # plt.plot(h_maj, gam_maj, color="black", label="Major " + str(azi))
    # plt.plot(h_min, gam_min, color="black", label="Minor " + str(azi + 90.0))
    # deltas = [22.5, 45, 67.5]
    # ndelta = len(deltas)
    # hd = np.zeros(ndelta)
    # gamd = np.zeros(ndelta)
    # color = iter(cm.plasma(np.linspace(0, 1, ndelta)))
    # for delta in deltas:
    #     index, hd, gamd, cov, ro = geostats.vmodel(nlag, xlag, azi + delta, vario)
    #     c = next(color)
    #     plt.plot(hd, gamd, color=c, label="Azimuth " + str(azi + delta))
    # plt.xlabel(r"Lag Distance $\bf(h)$, (m)")
    # plt.ylabel(r"$\gamma \bf(h)$")
    # plt.title("Interpolated NSCORE Porosity Variogram Models")
    # plt.xlim([0, hmaj1 * 1.5])
    # plt.ylim([0, 1.4])
    # plt.legend(loc="upper left")

    # plt.subplot(1, 3, 2)
    # plt.scatter(x1, y1, color="blue", edgecolors="black", s=sk_weights[0, 0] * 1000)
    # plt.scatter(x2, y2, color="red", edgecolors="black", s=sk_weights[0, 1] * 1000)
    # plt.scatter(x3, y3, color="green", edgecolors="black", s=sk_weights[0, 2] * 1000)
    # scatter = plt.scatter(
    #     500, 500, color="gray", edgecolors="black", s=(1 - sk_std) * 1000
    # )
    # ax = plt.gca()
    # plt.xlabel("X(m)")
    # plt.ylabel("Y(m)")
    # plt.title("Simple Kriging - Data and Unknown Locations")
    # plt.xlim([0, 1000])
    # plt.ylim([0, 1000])
    # for i, txt in enumerate(np.round(sk_weights[0], 2)):
    #     plt.annotate(txt, (x[i] + 20, y[i] + 20))
    # for i, txt in enumerate([1, 2, 3]):
    #     plt.annotate(txt, (x[i] - 40, y[i] - 40))
    # plt.annotate(np.round(sk_est[0], 2), (500 - 40, 500 - 40))
    # plt.annotate(
    #     "Mean Weight = " + str(np.round(1.0 - np.sum(sk_weights[0]), 2)), (20, 20)
    # )
    # plt.annotate("Unknown Location", (500 - 150, 500 + 50))

    # ellipse = Ellipse(
    #     (500, 500),
    #     width=hmin1 * 2.0,
    #     height=hmaj1 * 2.0,
    #     angle=360 - azi,
    #     facecolor="gray",
    #     alpha=0.1,
    # )
    # ax = plt.gca()
    # ax.add_patch(ellipse)

    # samples = norm.rvs(sk_est[0], sk_std, 1000, random_state=73073)
    # plt.subplot(1, 3, 3)
    # plt.hist(
    #     samples, bins=np.linspace(0, 4.0, 20), alpha=0.2, color="red", edgecolor="black"
    # )
    # plt.xlim([0.0, 4.0])
    # plt.ylim([0, 300])
    # plt.title("Uncertainty Model at Unknown Location")
    # plt.xlabel("Value")
    # plt.ylabel("Frequency")

    # ax = plt.gca()
    # ax.annotate("Simple Kriging Estimate = " + str(np.round(sk_est[0], 2)), (0.2, 17.5))
    # ax.annotate("Simple Kriging Variance = " + str(np.round(sk_var[0], 2)), (0.2, 6.5))
    # plt.subplots_adjust(
    #     left=0.0, bottom=0.0, right=2.2, top=0.9, wspace=0.3, hspace=0.3
    # )
    # plt.show()
