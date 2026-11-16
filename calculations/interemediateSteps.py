from pathlib import Path, PurePath
import re
import matplotlib.pyplot as plt

import matplotlib_format

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

from scipy.interpolate import CubicSpline
from scipy.signal import savgol_filter

def _zero_boundary_linear(y):
    """
    Impose y[0] = y[-1] = 0. The leading NaN run is filled by a straight line
    from 0 (at index 0) to the first valid value, and the trailing NaN run by a
    straight line from the last valid value to 0 (at the last index).
    """
    y = y.copy()
    good = np.where(~np.isnan(y))[0]
    if len(good) == 0:
        return y
    first, last = good[0], good[-1]
    if first > 0:
        y[:first] = np.linspace(0.0, y[first], first + 1)[:-1]
    if last < len(y) - 1:
        y[last + 1:] = np.linspace(y[last], 0.0, len(y) - last)[1:]
    return y

def _segment_edges(row, jump_tol, jump_factor, kink_tol):
    """Return the jump indices (discontinuity between i and i+1), including optional kinks."""
    diffs = np.abs(np.diff(row))
    tol = jump_tol if jump_tol is not None else \
        jump_factor * max(np.median(diffs), np.finfo(float).eps)
    jumps = set(np.where(diffs > tol)[0])

    if kink_tol is not None:  # kinks: |second difference| large -> isolate that sample
        d2 = np.abs(row[2:] - 2 * row[1:-1] + row[:-2])
        for i in np.where(d2 > kink_tol)[0]:      # kink at sample i+1
            jumps.update((i, i + 1))
    return np.array(sorted(jumps), dtype=int)


def _savgol_derivatives(row, dz, jump_idx, window, polyorder, trim):
    """
    Savitzky-Golay first and second derivatives, computed separately on each
    smooth segment. `trim` samples at each segment end are set to NaN, since
    one-sided fits are the least reliable (they get bridged afterwards).
    """
    m = len(row)
    d1 = np.full(m, np.nan)
    d2 = np.full(m, np.nan)
    edges = np.concatenate(([0], jump_idx + 1, [m]))
    for a, b in zip(edges[:-1], edges[1:]):
        seg = row[a:b]
        w = min(window, len(seg))
        if w % 2 == 0:
            w -= 1
        if w <= polyorder:
            continue                                  # too short: left as NaN
        d1[a:b] = savgol_filter(seg, w, polyorder, deriv=1, delta=dz, mode="interp")
        d2[a:b] = savgol_filter(seg, w, polyorder, deriv=2, delta=dz, mode="interp")
        if trim > 0:
            for d in (d1, d2):
                d[a:a + trim] = np.nan
                d[b - trim:b] = np.nan
    return d1, d2


def _fill_nans_cubic(y, window=4):
    """Fill interior NaN runs with a local cubic spline; hold the ends constant."""
    y = y.copy()
    bad = np.isnan(y)
    if bad.all() or not bad.any():
        return y
    idx = np.arange(len(y))
    good = idx[~bad]
    first, last = good[0], good[-1]
    y[:first] = y[first]
    y[last + 1:] = y[last]
    bad = np.isnan(y)
    i = first
    while i < last:
        if bad[i]:
            j = i
            while bad[j]:
                j += 1
            left = good[good < i][-window:]
            right = good[good >= j][:window]
            nodes = np.concatenate((left, right))
            y[i:j] = CubicSpline(nodes, y[nodes])(np.arange(i, j))
            i = j
        else:
            i += 1
    return y


def plot_iterates(filename, delimiter=",", kappa_start=1, cmap="viridis",
                  jump_tol=None, jump_factor=10.0, kink_tol=None,
                  sg_window=21, sg_polyorder=3, trim=None,
                  fill_gaps=True, fill_window=4,
                  yscale_deriv=None, yscale_deriv2=None):
    """
    Returns (figF, figdFdz, figd2Fdz2).

    sg_window    : Savitzky-Golay window (odd, in samples). Larger = smoother.
                   The 2nd derivative usually needs a much larger window than the 1st.
    sg_polyorder : polynomial order of the local fit (>= 2, typically 3 or 4)
    trim         : samples dropped at each segment end before bridging
                   (default sg_window // 2)
    kink_tol     : if set, also split where |f[i+1]-2f[i]+f[i-1]| > kink_tol
                   (for maps with discontinuous f', e.g. tent maps)
    """
    data = np.loadtxt(filename, delimiter=delimiter, ndmin=2)
    N, m = data.shape
    z = np.linspace(0, 1, m)
    dz = z[1] - z[0]
    kappas = np.arange(kappa_start, kappa_start + N)
    if trim is None:
        trim = sg_window // 2

    dF = np.empty_like(data)
    d2F = np.empty_like(data)
    Fplot = data.copy()

    for i, row in enumerate(data):
        jump_idx = _segment_edges(row, jump_tol, jump_factor, kink_tol)
        dF[i], d2F[i] = _savgol_derivatives(-row, dz, jump_idx, sg_window,
                                            sg_polyorder, trim)

        # f'' = 0 at z = 0 and z = 1: linear ramp from 0 to the first/last valid value
        dF[i] = _zero_boundary_linear(dF[i])
        d2F[i] = _zero_boundary_linear(d2F[i])

        if fill_gaps:
            dF[i] = _fill_nans_cubic(dF[i], fill_window)
            d2F[i] = _fill_nans_cubic(d2F[i], fill_window)   # only interior gaps remain
        Fplot[i, jump_idx] = np.nan

    norm = Normalize(vmin=kappas.min(), vmax=kappas.max())
    colormap = plt.get_cmap(cmap)
    sm = cm.ScalarMappable(norm=norm, cmap=colormap)

    def make_fig(Y, ylabel, yscale=None):
        fig, ax = plt.subplots()
        for k, y in zip(kappas, Y):
            ax.plot(z, y, color=colormap(norm(k)), lw=1)
        fig.colorbar(sm, ax=ax).set_label(r"$\kappa$")
        ax.set_xlabel(r"$z$")
        ax.set_ylabel(ylabel)
        if yscale:
            ax.set_yscale(yscale)
        return fig

    figF = make_fig(Fplot, r"$A_{0, \lambda}^{\text{br}, \kappa}(z)")
    figdFdz = make_fig(dF, r"$E_{\lambda}^\kappa$", yscale_deriv)
    figd2Fdz2 = make_fig(d2F, r"$\mathrm{d}^2 f^{\kappa}/\mathrm{d}z^2$", yscale_deriv2)
    return figF, figdFdz, figd2Fdz2



thingsRe = re.compile(r"([a-zA-Z0]*)_lambda_([\d.]*)_([a-z]+)")

dirs = Path("convergenceLog/Self-consistent solutions")

figOmega, axOmega = plt.subplots()
figA0, axA0 = plt.subplots()


figDict = {
            "eigenvalues":figOmega,
            "A0Induced":figA0
          }

axDict = {
            "eigenvalues":axOmega,
            "A0Induced":axA0
          }



axTitle = {
    "eigenvalues": r"$\omega^\kappa_\lambda$",
    "A0Induced":   r"$A_{0, \lambda}^{\text{br},\kappa}(1)$",
    # "eigenvalues": r"omega",
    # "A0Induced":   r"A",
}

for filename in list(dirs.glob("*"))[3:]:
    quantity, lambdaValue, bcs = re.findall(thingsRe, str(filename))[0]

    print(filename)
 
    lambdaValue = round(float(lambdaValue), 3)

    kappa = 1
    ax = axDict[quantity]
    fig = figDict[quantity]
    fig.suptitle(f"$\lambda = {lambdaValue}$")
    ax.cla()
    ax.set_xlabel("$\kappa$")
    ax.set_ylabel(axTitle[quantity])

    if quantity == "A0Induced":
        inv = -1
        figA0, figE, figrho = plot_iterates(filename)
        figA0.tight_layout()
        figE.tight_layout()
        figrho.tight_layout()
        figA0.savefig(f'convergenceLogFigs/A0(z)/{filename.stem}_function.pdf')
        figE.savefig(f'convergenceLogFigs/E(z)/{filename.stem}_function.pdf')
        figrho.savefig(f'convergenceLogFigs/rho(z)/{filename.stem}_rho.pdf')
    else:
        inv = 1
    with open(filename, "r") as myFile:
        for line in myFile:
            ax.scatter(kappa, inv*float(line.split(", ")[0]), color="#0E103D")
            kappa += 1

    fig.tight_layout()

    fig.savefig(f'convergenceLogFigs/{quantity}/{filename.stem}.pdf')

