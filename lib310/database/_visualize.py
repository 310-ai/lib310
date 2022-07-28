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
    bg_colors = []
    font_sizes = []
    for sector in fig.data[0]['ids']:
        if len(sector.split('/')) == 3:
            bg_colors.append(color_dataset(sector.split('/')[1]))
            font_sizes.append(18)
            continue
        bg_colors.append('transparent')
        font_sizes.append(14)

    fig.data[0]['marker']['colors'] = bg_colors
    fig.data[0].textposition = 'middle center'
    fig.data[0].texttemplate = "<span style='font-size: 18px;'><b> %{label} </b></span>"
    fig.data[0]['textfont']['color'] = ['#FFFFFF' for sector in fig.data[0]['ids'] if len(sector.split("/")) == 3]
    fig.update_layout(margin=dict(t=50, l=25, r=25, b=25))
    # fig.update_traces(hovertemplate='%{customdata[0]}')
    # fig.data[0].customdata = [
    #     v if fig.data[0].ids[i]
    #          in leaves else [""]
    #     for i, v in enumerate(fig.data[0].customdata)
    # ]
    fig.show()


def color_dataset(name):
    pallete = {
        'UNIPROT': '#D0534E',
        'GO': '#F29D38',
        'INTERPRO': '#AFD562',
        'STRINGS': '#51A2D1',
        'MID': '#693BD7',
        'METACLUST': '#A74CD8',
        'SANDBOX': 'CD507B'
    }
    if name not in pallete:
        return '#D0534E'
    return pallete[name]