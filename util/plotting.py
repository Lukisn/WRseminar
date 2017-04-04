#!/usr/bin/env python3

MARKER_SIZE = 5

NO_LINE = 0
SMALL_LINE = 0.75
MEDIUM_LINE = 1
THICK_LINE = 2

ANALYTIC_COLOR = "orange"
INITIAL_COLOR = "green"
FIXED_PIVOT_COLOR = "blue"
CELL_AVERAGE_COLOR = "red"

legend_style = {
    "loc": "best", "fontsize": "small", "numpoints": 1,
    "framealpha": 0.75, "fancybox": True
}

ana_style = {
    "linewidth": THICK_LINE, "linestyle": "solid", "color": ANALYTIC_COLOR,
    "marker": None, "markersize": MARKER_SIZE
}
ana_style0 = {
    "linewidth": THICK_LINE, "linestyle": "solid", "color": ANALYTIC_COLOR,
    "marker": "o", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
ana_style1 = {
    "linewidth": THICK_LINE, "linestyle": "solid", "color": ANALYTIC_COLOR,
    "marker": "s", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}

initial_style = {
    "linewidth": MEDIUM_LINE, "linestyle": "dotted", "color": INITIAL_COLOR,
    "marker": ".", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}

ca_style = {
    "linewidth": NO_LINE, "linestyle": "dashed", "color": CELL_AVERAGE_COLOR,
    "marker": "o", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
ca_style0 = {
    "linewidth": SMALL_LINE, "linestyle": "dotted", "color": CELL_AVERAGE_COLOR,
    "marker": "^", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
ca_style1 = {
    "linewidth": SMALL_LINE, "linestyle": "dotted", "color": CELL_AVERAGE_COLOR,
    "marker": "v", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}

fp_style = {
    "linewidth": NO_LINE, "linestyle": "dashed", "color": FIXED_PIVOT_COLOR,
    "marker": "x", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
fp_style0 = {
    "linewidth": SMALL_LINE, "linestyle": "dotted",  "color": FIXED_PIVOT_COLOR,
    "marker": "v", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
fp_style1 = {
    "linewidth": SMALL_LINE, "linestyle": "dotted", "color": FIXED_PIVOT_COLOR,
    "marker": "^", "markersize": MARKER_SIZE, "markerfacecolor": "None"
}
