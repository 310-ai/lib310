import matplotlib.pyplot as plt
import seaborn as sns
import math
from .charts.bubble_chart import BubbleChart
import plotly.express as px
import numpy as np
import pandas as pd



def visualize_all_datasets(df):
    bubble_chart = BubbleChart(area=[math.log(x) for x in df['num_rows']],
                               bubble_spacing=0.1)
    bubble_chart.collapse()

    fig, ax = plt.subplots(subplot_kw=dict(aspect="equal"))
    bubble_chart.plot(ax, list(df.index), sns.color_palette("Spectral", len(list(df.index))).as_hex())
    ax.axis("off")
    ax.relim()
    ax.autoscale_view()
    ax.set_title('Dataset (All)')

    plt.show()


def visualize_dataset(df, name):
    fig, ax = plt.subplots()
    ax.bar(df['table'], df['num_rows'], color=sns.color_palette("Spectral", len(df['table'])).as_hex())
    ax.axis("on")
    ax.relim()
    ax.autoscale_view()
    ax.set_title('Dataset (' + name + ')')
    plt.yscale('log')
    plt.show()


def visualize_all(df):
    df = pd.DataFrame({'310 db': ['310 DB'] * len(df['dataset']),
                       'dataset': df['dataset'],
                       'table': df['table'],
                       'size': df['size'],
                       'weight': np.log10(df['num_rows'])
                       })
    leaves = df.select_dtypes("object").apply("/".join, axis=1).values

    fig = px.treemap(df, path=['310 db', 'dataset', 'table'], values='weight', hover_data=["size"])
    # fig.update_traces(hovertemplate='%{customdata[0]}')
    # fig.data[0].customdata = [
    #     v if fig.data[0].ids[i]
    #          in leaves else [""]
    #     for i, v in enumerate(fig.data[0].customdata)
    # ]
    fig.show()
