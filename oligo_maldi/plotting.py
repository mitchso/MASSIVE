import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


# colour_dict_old = {'dark grey': '#3a3a3a',
#                'grey': '#747474',
#                'light grey': '#d9d9d9',
#                'lighter grey': '#f2f2f2',
#
#                'dark blue': '#4472c4',
#                'blue': '#5089ee',
#                'light blue': '#a6caec',
#
#                'dark red': '#af4c3a',
#                'red': '#da5e48',
#                'light red': '#ff6f55',
#
#                'dark green': '#006837',
#                'green': '#009245',
#                'light green': '#8cc63f'}


colour_dict = {'dark grey': '#333333',
               'grey': '#737373',
               'light grey': '#b3b3b3',
               'lighter grey': '#e6e6e6',

               'dark blue': '#007aed',
               'blue': '#5fa8ed',
               'light blue': '#bed6ed',

               'dark red': '#ff2700',
               'red': '#ff7d66',
               'light red': '#ffd4cc',

               'dark green': '#44cc00',
               'green': '#7acc52',
               'light green': '#b1cca3'}


def group_constructor(rows, cols, groups, group_by: str='column'):
    """Creates a dictionary.
        Keys = groups
        Values = lists of wells associated with group.
         Wells are grouped by column"""
    d = {}

    for i, g in enumerate(groups):
        if group_by == 'column':
            wells = [row + str(cols[i]) for row in rows]
        elif group_by == 'row':
            wells = [rows[i] + str(col) for col in cols]
        else:
            raise ValueError("group_by must be row or column")

        d[g] = wells


    return d


def stacked_bar(df,
                x_order: list, x_labels: list,
                y_order: list, y_colours: list, y_labels: list | None=None,
                norm_to: str | None = None, mode: str = 'relative',
                sort_by: str | None = None, sort_order: str = 'ascending',
                errorbars: bool = False,
                title: str | None = None) -> plt.axes:
    """
    :param df: dataframe, must contain all wells (rows) from x_order and all columns from y_order
    :param x_order: nested list of wells to include, in the desired order on the x-axis
    :param x_labels: categorical x-axis labels
    :param y_order: list of columns to include, in the desired order on the y-axis (bottom to top)
    :param y_colours: bar colours (bottom to top)
    :param y_labels: option to change the labels for the bars
    :param norm_to: column that all other columns will be normalized to
    :param mode: 'relative' will present data of all columns relative to 'norm_to'.
                 'fraction' will present data of all columns as a fraction of 100%.
    :param title: option to add a title.
    :return: plt.fig, plt.ax
    """

    # make sure that the correct params are given.
    if mode == 'relative':
        if norm_to is None:
            raise ValueError("parameter \'norm_to\' must be provided if mode is \'relative\'")

    # pull rows that correspond to wells
    df_slice = pd.DataFrame()   # creates an empty dataframe to populate
    wells = sum(x_order, [])    # creates a list of all wells needed
    for well in wells:
        row = df.loc[df.well == well]
        df_slice = pd.concat([df_slice, row])

    # create a nested dictionary with arrangement: x_category -> moi -> values and avg
    ddict = {}
    i = 0
    for group in x_order:
        if mode == 'fraction':  # All values will be expressed as a percentage of the total signal
            cols = y_order
            rows = df_slice.loc[df_slice['well'].isin(group)][cols]
            rows['total'] = rows.sum(axis=1)
            rows = rows.div(rows['total'], axis=0)

        elif mode == 'relative':    # All values will be expressed relative to 'norm_to'
            cols = y_order + [norm_to]
            rows = df_slice.loc[df_slice['well'].isin(group)][cols]
            rows = rows.div(rows[norm_to], axis=0)

        else:
            raise ValueError("mode must be \'fraction\' or \'relative\'")

        x = x_labels[i]
        ddict[x] = {}
        for moi in y_order:
            ddict[x][moi] = {}
            ddict[x][moi]['vals'] = list(rows[moi])
            ddict[x][moi]['avg'] = np.average(ddict[x][moi]['vals'])
            ddict[x][moi]['std'] = np.std(ddict[x][moi]['vals'])
        i += 1

    if sort_by:
        if sort_order == 'ascending' or sort_order is None:
            reverse = False
        elif sort_order == 'descending':
            reverse = True
        else:
            raise ValueError("sort_order must be \'ascending\', \'descending\', or \'None\'")
        # Sort dictionary based on the avg value of MOI designated as 'sort_by'
        ddict = {k: v for k, v in sorted(ddict.items(), key=lambda item: item[1][sort_by]['avg'], reverse=reverse)}
        # Update labels to reflect the sorting
        x_labels = [k for k in ddict.keys()]

    # Generate plot
    plt.rcParams.update({'font.family': 'Arial',
                         'font.size': 13,
                         'text.color': 'black',
                         'axes.labelcolor': 'black'})

    sns.set_style('ticks')  # Necessary to see minor ticks

    fig, ax = plt.subplots()

    # populate bar chart by calculating means
    bottom = [0] * len(x_labels)    # starting point for the bar
    for i, moi in enumerate(y_order):
        if y_labels:
            label = y_labels[i]
        else:
            label = moi
        avg = [ddict[x][moi]['avg'] for x in x_labels]
        std = [ddict[x][moi]['std'] for x in x_labels]

        bar = {'x': x_labels,
               'height': avg,
               'bottom': bottom,
               'label': label,
               'color': y_colours[i],
               'edgecolor': 'black',
               'linewidth': 0.5}

        if errorbars:
            bar['yerr'] = std   # default error bars are 1 standard deviation
            bar['capsize'] = 3

        plt.bar(**bar)

        bottom = [v1 + v2 for v1, v2 in zip(bottom, avg)]   # update bottom for next bar

    # add individual data points as a scatter
    scatter_df = pd.DataFrame(columns=['x', 'moi', 'val'])  # df with 1 row per point to plot
    for x in x_labels:
        bottom = 0
        for moi in y_order:
            vals = ddict[x][moi]['vals']
            vals = [v + bottom for v in vals]
            bottom = bottom + ddict[x][moi]['avg']
            entry = pd.DataFrame({'x': x, 'moi': moi, 'val': vals})
            scatter_df = pd.concat([scatter_df, entry], ignore_index=True)

    sns.swarmplot(data=scatter_df, x='x', y='val',
                  hue='moi', palette=y_colours,
                  linewidth=0.5, size=3,
                  edgecolor='black', legend=False, clip_on=False)

    # ticks
    if mode == 'fraction':
        plt.ylim(0, 1)
        plt.yticks(ticks=[0, 0.25, 0.5, 0.75, 1.0], labels=['0%', '25%', '50%', '75%', '100%'])
        plt.yticks(ticks=np.linspace(0, 1, num=20, endpoint=False), minor=True)
        ax.yaxis.get_ticklocs(minor=True)
        ax.minorticks_on()
        ax.xaxis.set_tick_params(which='minor', bottom=False)   # turn off x-axis minor ticks

    plt.xticks(rotation=45, ha='right')
    ax.tick_params(color='black', labelcolor='black')
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
    plt.grid(color='black', axis='y', linewidth=0.5)

    #labels, legend, title
    plt.xlabel('Condition')
    plt.ylabel('Signal intensity (au)')
    if title:
        plt.title(title, pad=20)

    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', reverse=True)
    # plt.tight_layout()

    return fig, ax

