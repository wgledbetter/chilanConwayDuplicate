from geostatspy.GSLIB import varmapv_2d
from pyrcz_kriging import *

import numpy as np

import plotly.offline as py
import plotly.graph_objs as go


def main():
    ## Base Parameters
    nug = 0.01
    it1 = "Spherical"
    azi = 1
    hmaj1 = 5
    hmin1 = 5

    xs = [-1, -1, 1, 1, 0, 0.25, -0.5, -0.45, -0.55]
    ys = [-1, 1, -1, 1, 0, -0.25, 0.5, -0.45, -0.55]
    vs = [1, 3, 3, -1, 4, -3, -1, -2, 2]

    ############################################################################
    ## Query

    nQ = 300

    X = np.linspace(-0.99, 0.99, nQ)
    Y = np.linspace(-0.99, 0.99, nQ)

    xqs = []
    yqs = []
    for x in X:
        for y in Y:
            xqs.append(x)
            yqs.append(y)

    estim, var, weights = f_make_krige(
        nug, it1, azi, hmaj1, hmin1, xs, ys, vs, xqs, yqs
    )

    ############################################################################
    ## Plotting
    krigTrace = go.Contour(
        x=xqs,
        y=yqs,
        z=estim,
        colorscale="Earth",
        showscale=False,
        line_width=3,
        colorbar=dict(tickangle=0, tickfont=dict(size=15)),
        contours=dict(showlabels=True, labelfont=dict(size=20, color="black")),
    )
    sampleTrace = go.Scatter(x=xs, y=ys, mode="markers", marker=dict(size=8), text=vs)
    layout = go.Layout(
        xaxis=dict(
            fixedrange=False,
            tickfont=dict(size=20),
            title_text="X",
            title_font=dict(size=20),
        ),
        yaxis=dict(
            fixedrange=False,
            tickfont=dict(size=20),
            title_text="Y",
            title_font=dict(size=20),
        ),
        showlegend=True,
    )
    fig = go.Figure(data=[krigTrace, sampleTrace], layout=layout)
    py.plot(fig, filename="./krig.html")
    fig.write_image("./krig.pdf", width=1000, height=650)

    # surfTrace = go.Surface(x=xqs, y=yqs, z=estim)
    # sample3d = go.Scatter3d(x=xs, y=ys, z=vs, mode="markers", marker=dict(size=8))
    # fig2 = go.Figure(data=[surfTrace, sample3d])
    # py.plot(fig2, filename="./krig3d.html")


if __name__ == "__main__":
    main()
