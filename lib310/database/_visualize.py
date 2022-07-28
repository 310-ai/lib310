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
                       'weight': df['num_rows'].apply(lambda x: math.log10(x) - 4 if math.log10(x) - 4 > 1 else 1),
                       'rows': [number_to_abbrevation(n) for n in df['num_rows']]
                       })
    leaves = df.select_dtypes("object").apply("/".join, axis=1).values

    fig = px.treemap(df, path=['310 db', 'dataset', 'table'], values='weight', hover_data=["size", "rows"])
    bg_colors = []
    for sector in fig.data[0]['ids']:
        if len(sector.split('/')) == 3:
            bg_colors.append(color_dataset(sector.split('/')[1]))
            continue
        bg_colors.append('transparent')

    fig.data[0]['marker']['colors'] = bg_colors
    fig.data[0].textposition = 'middle center'
    fig.data[0]['marker']['line']['color'] = '#D3D3D3'
    fig.data[0].texttemplate = "<span style='font-size: 16px;'> %{label} </span> <br> <span style='font-size: 12px;'>%{customdata[0]} <br> %{customdata[1]}<span>"
    fig.data[0]['textfont']['color'] = ['#FFFFFF' for sector in fig.data[0]['ids'] if len(sector.split("/")) == 3]
    fig.update_layout(margin=dict(t=50, l=25, r=25, b=25))
    # fig.update_traces(hovertemplate='%{customdata[0]}')
    # fig.data[0].customdata = [
    #     v if fig.data[0].ids[i]
    #          in leaves else [""]
    #     for i, v in enumerate(fig.data[0].customdata)
    # ]
    fig.show()


def number_to_abbrevation(num):
    label = ['', 'K', 'M', 'B', 'T']
    i = 0
    while num > 1000:
        num /= 1000
        i += 1
    return f'{num:.1f} {label[i]}'

def color_dataset(name):
    pallete = {
        'UNIPROT': '#E04848',
        'INTERPRO': '#FF9900',
        'STRINGS': '#A7DD4F',
        'GO': '#2FA7DB',
        'METACLUST': '#7330DF',
        'MID': '#693BD7',
        'SANDBOX': 'CD507B'
    }
    if name not in pallete:
        return '#D0534E'
    return pallete[name]