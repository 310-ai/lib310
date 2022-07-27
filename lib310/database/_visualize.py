import matplotlib.pyplot as plt
import seaborn as sns
import math
from .charts.bubble_chart import BubbleChart


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
    fig = plt.figure()
    ax = fig.add_axes([0.06, 0.06, 0.94, 0.94])
    ax.bar(df['table'], df['num_rows'], color=sns.color_palette("Spectral", len(df['table'])).as_hex())
    ax.axis("on")
    ax.relim()
    ax.autoscale_view()
    ax.set_title('Dataset (' + name + ')')
    ax.set_yscale('log')
    plt.show()
