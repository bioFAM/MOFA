import numpy as np

from bokeh.plotting import figure, show, output_file

from bokeh.charts import Scatter, output_file, show, HeatMap, BoxPlot
from bokeh.models import ColumnDataSource
from bokeh.models import HoverTool

import seaborn as sns
import matplotlib.pyplot as plt

# from bokeh.sampledata.autompg import autompg as df

# TODO look into crossfilter
# TODO look into the widgets to change colour interactivly etc
# TODO look at other plots to do (heatmap etc)
# TODO for continuous colors use matplotlib with bokeh
# TODO get density plots on the side of teh scatter plot
def scatter_plot(data, _x, _y, identifier=None, size=None, invert_size=False):

    _data = data.copy()
    # invert size scale
    if size is not None and invert_size:
        _data[size] = 1./_data[size]
    # rescale size
    if size is not None:
        _data[size] /= _data[size].max()

    # configure hoover tool
    _data =_data.reset_index(level=['protein'])
    _data =_data.reset_index(level=['image'])
    _data =_data.reset_index(level=['normalisation'])
    _data =_data.reset_index(level=['focal_type'])
    _data =_data.reset_index(level=['env_type'])

    tooltips=[
            ("index", "$index"),
            ("(mean,LLR)", "($x, $y)"),
            ('image', "@image"),
            ('protein', "@protein"),
            ('normalisation', "@normalisation"),
            ('focal type', "@focal_type"),
            ('env type', "@env_type")
        ]

    # TODO figure out how to add alpha value
    p = Scatter(_data, x=_x, y=_y, color=identifier, legend="top_right",
                tooltips=tooltips, toolbar_location='above')

    # output_file("scatter.html")

    show(p)
