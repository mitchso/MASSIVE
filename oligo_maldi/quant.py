import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import warnings
import sys
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action="ignore", category=UserWarning)

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


def get_aa_position(label):
    """
    returns integer value of amino acid position from single character labels.
    EXAMPLE:    "D83A" return int(83)
                "S340P" returns int(340)
                "blah" returns 0
    """
    try:
        pos = int(label[1:-1])
        return pos
    except ValueError:
        return 0


class DataProcessor:
    """
    Class built for processing the output of Experiment.write_to_excel().
    Multiple datasets can be processed together by passing a dictionary of dataframes.
    """
    def __init__(self, dataframes: dict,    #
                 attributes: list,
                 collapse_mois: dict | None = None):

        """

        :param dataframes: keys are whatever you would like to name each dataset, each value is a pandas dataframe
                                i.e. frames = {'plate1': plate1_dataframe,
                                                'plate2': plate2_dataframe}

        :param attributes: a list of ALL experimental conditions which are shared across technical replicates.
                            these should correspond to column headers in the dataframes
                            it is CRITICAL to include ALL otherwise samples will be lumped together that should not be
                                i.e. attributes = ['ntp', 'enz_name', 'enz_id']

        :param collapse_mois: keys are the new name for an MOI, each value is a list of MOIs that should be renamed
                                i.e. collapse = {'+N, PPP': ['+A, PPP', '+G, PPP', '+C, PPP', '+T, PPP']}

        """

        self.dataframes = dataframes
        self.data_ungrouped = self.collapse_input_dataframes()
        self.data_mois_collapsed = self.collapse_mois(df=self.data_ungrouped, mapping_dict=collapse_mois)
        self.data_replicates_averaged = self.collapse_replicates(self.data_mois_collapsed, attributes)
        self.data_activity_scores = None    # TODO, just score in the dataframe nothing else

    def collapse_input_dataframes(self):
        """All dataframes from self.dataframes will now be in self.data_ungrouped,
        with an additional column "source" referencing the name of the original dataframe"""

        df_list = []
        for k, v in self.dataframes.items():
            new_v = v.copy()
            new_v.insert(0, 'source', k)
            df_list.append(new_v)

        collapsed = pd.concat(df_list, ignore_index=True)
        return collapsed

    def collapse_mois(self, df: pd.DataFrame, mapping_dict: dict) -> pd.DataFrame:
        """
        :param mapping_dict: keys are the new name for an MOI, each value is a list of MOIs that should be renamed
                                i.e. collapse = {'+N, PPP': ['+A, PPP', '+G, PPP', '+C, PPP', '+T, PPP']}
        :return: dataframe where all data for MOIs to be removed (values in mapping_dict) are placed
                    under the new column headers (keys in mapping_dict)
        """

        if mapping_dict is None:
            return df

        new_df = df.copy()  # make a new dataframe

        for i, row in df.iterrows():    # for each row in dataframe
            for new_key, old_keys in mapping_dict.items():
                for old_key in old_keys:
                    # fina actual values
                    if not np.isnan(row[old_key]):  # evaluates to true if there is a real number in the column associated with old_key
                        # place the number with the new key in the new df
                        new_df.at[i, new_key] = row[old_key]    # i = row index

        for old_keys in mapping_dict.values():  # remove all the old_keys from the new dict
            new_df.drop(columns=old_keys, inplace=True)

        return new_df

    def collapse_replicates(self, df: pd.DataFrame, attributes: list):
        """Averages technical replicates into average values.
        Requires a list of attributes that should be shared across all replicates."""
        return df.groupby(attributes).mean(numeric_only=True)

    def activity_score(self):
        """
        Score must be >0 if the enzyme does something
        Going slower but being highly specific is better
        All nonspecific products can be scored the same

        + for desired product
        reacted substrate = 1 - unreacted substrate


        :return:
        """

        return 0


    def full_figure(self):

        # fig, axs = subplots() then call other functions for each subplot
        pass

    def get_stacked_bar_ddict(self,
                              df: pd.DataFrame,
                              x_wells: dict,
                              y_labels: list,
                              x_overwrite: list | None,
                              y_overwrite: list | None) -> dict:
        """
        :param df: Dataframe to pull from
        :param x_wells: Nested dictionary
                            # keys = x_labels
                            # values = {source: which dataframe the data is located in
                            #           wells: which wells the data is from
        :param y_labels: Columns to collect data from
        :param x_overwrite: Optional new names for x_labels
        :param y_overwrite: Optional new names for column data
        :return: nested dictionary of the following style:
                    ddict = {
                    'x_label_1': {
                    'MOI_1': {'vals': [data], 'avg': avg of data},
                    'MOI_2': {'vals': [data], 'avg': avg of data}
                    },
                    'x_label_2': {
                    'MOI_1': {'vals': [data], 'avg': avg of data} ,
                    'MOI_2': {'vals': [data], 'avg': avg of data}
                    },
                    etc
                    }
        """
        ddict = {}
        for i, (x_label, label_info) in enumerate(x_wells.items()):
            source = label_info['source']
            wells = label_info['wells']
            cols = y_labels

            rows = df.loc[df['source'] == source] # slice only the correct source plate
            rows = rows.loc[df['well'].isin(wells)] #slice only the correct wells for a particular x_label
            rows = rows[cols]   # slice only the columns corresponding to MOIs you want to work with
            rows['total'] = rows.sum(axis=1)    # total ion intensity across all MOIs for a row
            rows = rows.div(rows['total'], axis=0)  # converts every column to a percentage of total signal

            if x_overwrite:
                x_label = x_overwrite[i]    # change name to the new label from y_overwrite

            ddict[x_label] = {}

            for j, moi in enumerate(y_labels):
                if y_overwrite:
                    moi = y_overwrite[j]    # change name to the new label from y_overwrite

                ddict[x_label][moi] = {}
                ddict[x_label][moi]['vals'] = list(rows[y_labels[j]])
                ddict[x_label][moi]['avg'] = np.average(ddict[x_label][moi]['vals'])

        return ddict

    def sort_ddict(self, ddict: dict, sort_by: list, sort_order: str) -> dict:
        """
        :param ddict: see self.get_stacked_bar_ddict() for more info
        :param sort_by: list of MOIs to sort by, in order of priority
        :param sort_order: whether to sort by having the most or the least of each MOI
        :return:
        """

        if sort_order == 'ascending' or sort_order is None:
            reverse = False
        elif sort_order == 'descending':
            reverse = True
        else:
            raise ValueError("sort_order must be \'ascending\', \'descending\', or \'None\'")

        # First, sort dictionary alphabetically
        ddict = {k: v for k, v in sorted(ddict.items(), key=lambda item: item[0])}  # item[0] corresponds to the x_label

        # Then, sort dictionary by amino acid position
        ddict = {k: v for k, v in sorted(ddict.items(), key=lambda item: get_aa_position(item[0]))} # get_aa_pos defaults to 0 if the x_label is not in the form 'A28P'

        # check if it's sorting on a single MOI, if it is, put it in a list to standardize next step
        if isinstance(sort_by, str):
            sort_by = [sort_by]

        # if it's multiple MOIs, reverse them so you sort by the lowest priority first
        # this results in a final order that prioritizes the first item in sort_by
        elif isinstance(sort_by, list):
            sort_by.reverse()

        # Finally, sort dictionary based on the avg value of MOI designated as 'sort_by'
        for moi in sort_by:
            ddict = {k: v for k, v in
                     sorted(ddict.items(), key=lambda item: item[1][moi]['avg'], reverse=reverse)}

        return ddict

    def stacked_bar(self,
                    global_var: tuple,  # (column, column_value) tuple. Only plots data for rows that have column_value
                    x_category: str,    # which categorical variable to plot as individual samples on the x axis
                    y_species: list,    # which MOIs to show on the y-axis
                    y_colours: list,    # colours for each MOI in the order they appear in y_species
                    x_labels: list | None = None,   # can select a specific subset of x_category
                    x_labels_overwrite: list | None = None, # allows for changing the x_labels
                    y_legend_overwrite: list | None = None, # allows for changing the y_labels
                    sort_by: list | None = None,    # list of MOIs to sort by, in order of priority
                    sort_order: str | None = None, # whether to sort by having the most or the least of each MOI
                    figsize=None,
                    title=None,
                    show_legend=True,
                    ax=None):

        df = self.data_mois_collapsed   # start with a dataframe where all MOIs are named the same thing
        df = df.loc[df[global_var[0]] == global_var[1]] # take the slice that corresponds to global_var

        if x_labels is None:    # if no specific labels have been provided, get all labels from x_category
            x_labels = list(set(df[x_category].tolist()))

        # x_wells is a nested dict.
        # keys = x_labels
        # values = {source: which dataframe the data is located in
        #           wells: which wells the data is from
        x_wells = {}
        for x_label in x_labels:
            try:
                source = df.loc[df[x_category] == x_label]['source'].unique()[0]
            except IndexError:
                print(f"IndexError: '{x_label}' is not a column in the DataFrame.")
                sys.exit(0)

            wells = df.loc[df[x_category] == x_label]['well'].tolist()
            x_wells[x_label] = {'source': source,
                                'wells': wells}

        ddict = self.get_stacked_bar_ddict(df=df,
                                           x_wells=x_wells,
                                           x_overwrite=x_labels_overwrite,
                                           y_labels=y_species,
                                           y_overwrite=y_legend_overwrite)

        if sort_by:
            ddict = self.sort_ddict(ddict=ddict, sort_by=sort_by, sort_order=sort_order)

        # if an ax is specified when this is called, that implies that a Figure already exists, so only make one if ax is not called
        if not ax:
            fig, ax = plt.subplots()

        ax = self.generate_stacked_bar_ax(ddict=ddict,
                                          y_order=y_species if not y_legend_overwrite else y_legend_overwrite,
                                          y_colours=y_colours,
                                          figsize=figsize,
                                          title=title,
                                          show_legend=show_legend,
                                          ax=ax)
        return ax


    def generate_stacked_bar_ax(self,
                                # x_labels: list,
                                y_order: list,
                                y_colours: list,
                                ddict: dict,
                                figsize: tuple,
                                title: str,
                                show_legend: bool,
                                ax=None) -> plt.axes:

        # Generate plot
        plt.rcParams.update({'font.family': 'Arial',
                             'font.size': 13,
                             'text.color': 'black',
                             'axes.labelcolor': 'black'})

        sns.set_style('ticks')  # Necessary to see minor ticks

        if not ax:
            fig, ax = plt.subplots(figsize=figsize)

        # populate bar chart by calculating means
        x_labels = ddict.keys()
        bottom = [0] * len(x_labels)  # starting point for the bar
        for i, moi in enumerate(y_order):
            avg = [ddict[x][moi]['avg'] for x in x_labels]

            bar = {'x': x_labels,
                   'height': avg,
                   'bottom': bottom,
                   'label': moi,
                   'color': y_colours[i],
                   'edgecolor': 'black',
                   'linewidth': 0.5}

            ax.bar(**bar)

            bottom = [v1 + v2 for v1, v2 in zip(bottom, avg)]  # update bottom for next bar

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
                      edgecolor='black', legend=False, clip_on=False, ax=ax)

        # ticks
        ax.set_ylim(0, 1)
        ax.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1.0], labels=['0%', '25%', '50%', '75%', '100%'])
        ax.set_yticks(ticks=np.linspace(0, 1, num=20, endpoint=False), minor=True)
        ax.yaxis.get_ticklocs(minor=True)
        ax.minorticks_on()
        ax.xaxis.set_tick_params(which='minor', bottom=False)  # turn off x-axis minor ticks

        ax.set_xticklabels(labels=x_labels, ha='right', rotation=45)
        # ax.tick_params(axis='x', labelrotation=45) #, ha='right')
        ax.tick_params(color='black', labelcolor='black')
        for spine in ax.spines.values():
            spine.set_edgecolor('black')
        ax.grid(color='black', axis='y', linewidth=0.5)

        # labels, legend, title
        ax.set_xlabel('')
        ax.set_ylabel('Signal intensity (au)')
        if title:
            ax.set_title(title, pad=20)
        if show_legend:
            ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', reverse=True)
        # plt.tight_layout()

        return ax




#     def stacked_bar(self,
#                     x_order: list, x_labels: list,
#                     y_order: list, y_colours: list, y_labels: list | None = None,
#                     norm_to: str | None = None, mode: str = 'relative',
#                     sort_by: str | None | list = None, sort_order: str = 'ascending',
#                     errorbars: bool = False,
#                     title: str | None = None,
#                     figsize=None) -> plt.axes:
#         """
#         :param df: dataframe, must contain all wells (rows) from x_order and all columns from y_order
#         :param x_order: nested list of wells to include, in the desired order on the x-axis
#         :param x_labels: categorical x-axis labels
#         :param y_order: list of columns to include, in the desired order on the y-axis (bottom to top)
#         :param y_colours: bar colours (bottom to top)
#         :param y_labels: option to change the labels for the bars
#         :param norm_to: column that all other columns will be normalized to
#         :param mode: 'relative' will present data of all columns relative to 'norm_to'.
#                      'fraction' will present data of all columns as a fraction of 100%.
#         :param sort_by: designates which MOI to sort the columns by. You can designate 1 as a string,
#                         or multiple as a list. If multiple, the function will prioritize starting from
#                         the first, ending at the last.
#         :param title: option to add a title.
#         :return: plt.fig, plt.ax
#         """
#
#         # make sure that the correct params are given.
#         if mode == 'relative':
#             if norm_to is None:
#                 raise ValueError("parameter \'norm_to\' must be provided if mode is \'relative\'")
#
#         # pull rows that correspond to wells
#         df_slice = pd.DataFrame()  # creates an empty dataframe to populate
#         wells = sum(x_order, [])  # creates a list of all wells needed
#         for well in wells:
#             row = self.data_ungrouped.loc[self.data_ungrouped.well == well]
#             df_slice = pd.concat([df_slice, row])
#
#         # create a nested dictionary with arrangement: x_category -> moi -> values and avg
#         ddict = {}
#         i = 0
#         for group in x_order:
#             if mode == 'fraction':  # All values will be expressed as a percentage of the total signal
#                 cols = y_order
#                 rows = df_slice.loc[df_slice['well'].isin(group)][cols]
#                 rows['total'] = rows.sum(axis=1)
#                 rows = rows.div(rows['total'], axis=0)
#
#             elif mode == 'relative':  # All values will be expressed relative to 'norm_to'
#                 cols = y_order + [norm_to]
#                 rows = df_slice.loc[df_slice['well'].isin(group)][cols]
#                 rows = rows.div(rows[norm_to], axis=0)
#
#             else:
#                 raise ValueError("mode must be \'fraction\' or \'relative\'")
#
#             x = x_labels[i]
#             ddict[x] = {}
#             for moi in y_order:
#                 ddict[x][moi] = {}
#                 ddict[x][moi]['vals'] = list(rows[moi])
#                 ddict[x][moi]['avg'] = np.average(ddict[x][moi]['vals'])
#                 ddict[x][moi]['std'] = np.std(ddict[x][moi]['vals'])
#             i += 1
#
#         if sort_by:
#             if sort_order == 'ascending' or sort_order is None:
#                 reverse = False
#             elif sort_order == 'descending':
#                 reverse = True
#             else:
#                 raise ValueError("sort_order must be \'ascending\', \'descending\', or \'None\'")
#
#             # First, sort dictionary alphabetically
#             ddict = {k: v for k, v in sorted(ddict.items(), key=lambda item: item[0])}
#
#             # Then, sort dictionary by amino acid position
#             ddict = {k: v for k, v in sorted(ddict.items(), key=lambda item: get_aa_position(item[0]))}
#
#             # check if it's sorting on a single MOI, if it is, put it in a list to standardize next step
#             if isinstance(sort_by, str):
#                 sort_by = [sort_by]
#
#             # if it's multiple MOIs, reverse them so you sort by the lowest priority first
#             # this results in a final order that prioritizes the first item in sort_by
#             elif isinstance(sort_by, list):
#                 sort_by.reverse()
#
#             # Finally, sort dictionary based on the avg value of MOI designated as 'sort_by'
#             for moi in sort_by:
#                 ddict = {k: v for k, v in
#                          sorted(ddict.items(), key=lambda item: item[1][moi]['avg'], reverse=reverse)}
#
#             # Update labels to reflect the sorting
#             x_labels = [k for k in ddict.keys()]
#
#         # Generate plot
#         plt.rcParams.update({'font.family': 'Arial',
#                              'font.size': 13,
#                              'text.color': 'black',
#                              'axes.labelcolor': 'black'})
#
#         sns.set_style('ticks')  # Necessary to see minor ticks
#
#         fig, ax = plt.subplots(figsize=figsize)
#
#         # populate bar chart by calculating means
#         bottom = [0] * len(x_labels)  # starting point for the bar
#         for i, moi in enumerate(y_order):
#             if y_labels:
#                 label = y_labels[i]
#             else:
#                 label = moi
#             avg = [ddict[x][moi]['avg'] for x in x_labels]
#             std = [ddict[x][moi]['std'] for x in x_labels]
#
#             bar = {'x': x_labels,
#                    'height': avg,
#                    'bottom': bottom,
#                    'label': label,
#                    'color': y_colours[i],
#                    'edgecolor': 'black',
#                    'linewidth': 0.5}
#
#             if errorbars:
#                 bar['yerr'] = std  # default error bars are 1 standard deviation
#                 bar['capsize'] = 3
#
#             plt.bar(**bar)
#
#             bottom = [v1 + v2 for v1, v2 in zip(bottom, avg)]  # update bottom for next bar
#
#         # add individual data points as a scatter
#         scatter_df = pd.DataFrame(columns=['x', 'moi', 'val'])  # df with 1 row per point to plot
#         for x in x_labels:
#             bottom = 0
#             for moi in y_order:
#                 vals = ddict[x][moi]['vals']
#                 vals = [v + bottom for v in vals]
#                 bottom = bottom + ddict[x][moi]['avg']
#                 entry = pd.DataFrame({'x': x, 'moi': moi, 'val': vals})
#                 scatter_df = pd.concat([scatter_df, entry], ignore_index=True)
#
#         sns.swarmplot(data=scatter_df, x='x', y='val',
#                       hue='moi', palette=y_colours,
#                       linewidth=0.5, size=3,
#                       edgecolor='black', legend=False, clip_on=False)
#
#         # ticks
#         if mode == 'fraction':
#             plt.ylim(0, 1)
#             plt.yticks(ticks=[0, 0.25, 0.5, 0.75, 1.0], labels=['0%', '25%', '50%', '75%', '100%'])
#             plt.yticks(ticks=np.linspace(0, 1, num=20, endpoint=False), minor=True)
#             ax.yaxis.get_ticklocs(minor=True)
#             ax.minorticks_on()
#             ax.xaxis.set_tick_params(which='minor', bottom=False)  # turn off x-axis minor ticks
#
#         plt.xticks(rotation=45, ha='right')
#         ax.tick_params(color='black', labelcolor='black')
#         for spine in ax.spines.values():
#             spine.set_edgecolor('black')
#         plt.grid(color='black', axis='y', linewidth=0.5)
#
#         # labels, legend, title
#         plt.xlabel('')
#         plt.ylabel('Signal intensity (au)')
#         if title:
#             plt.title(title, pad=20)
#
#         plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', reverse=True)
#         # plt.tight_layout()
#
#         return fig, ax
#
#
#
# def group_constructor(rows, cols, groups, group_by: str='column'):
#     """Creates a dictionary.
#         Keys = groups
#         Values = lists of wells associated with group.
#          Wells are grouped by column"""
#     d = {}
#
#     for i, g in enumerate(groups):
#         if group_by == 'column':
#             wells = [row + str(cols[i]) for row in rows]
#         elif group_by == 'row':
#             wells = [rows[i] + str(col) for col in cols]
#         else:
#             raise ValueError("group_by must be row or column")
#
#         d[g] = wells
#
#
#     return d
#

#
#
# def stacked_bar(df,
#                 x_order: list, x_labels: list,
#                 y_order: list, y_colours: list, y_labels: list | None=None,
#                 norm_to: str | None = None, mode: str = 'relative',
#                 sort_by: str | None | list = None, sort_order: str = 'ascending',
#                 errorbars: bool = False,
#                 title: str | None = None,
#                 figsize = None) -> plt.axes:
#     """
#     :param df: dataframe, must contain all wells (rows) from x_order and all columns from y_order
#     :param x_order: nested list of wells to include, in the desired order on the x-axis
#     :param x_labels: categorical x-axis labels
#     :param y_order: list of columns to include, in the desired order on the y-axis (bottom to top)
#     :param y_colours: bar colours (bottom to top)
#     :param y_labels: option to change the labels for the bars
#     :param norm_to: column that all other columns will be normalized to
#     :param mode: 'relative' will present data of all columns relative to 'norm_to'.
#                  'fraction' will present data of all columns as a fraction of 100%.
#     :param sort_by: designates which MOI to sort the columns by. You can designate 1 as a string,
#                     or multiple as a list. If multiple, the function will prioritize starting from
#                     the first, ending at the last.
#     :param title: option to add a title.
#     :return: plt.fig, plt.ax
#     """
#
#     # make sure that the correct params are given.
#     if mode == 'relative':
#         if norm_to is None:
#             raise ValueError("parameter \'norm_to\' must be provided if mode is \'relative\'")
#
#     # pull rows that correspond to wells
#     df_slice = pd.DataFrame()   # creates an empty dataframe to populate
#     wells = sum(x_order, [])    # creates a list of all wells needed
#     for well in wells:
#         row = df.loc[df.well == well]
#         df_slice = pd.concat([df_slice, row])
#
#     # create a nested dictionary with arrangement: x_category -> moi -> values and avg
#     ddict = {}
#     i = 0
#     for group in x_order:
#         if mode == 'fraction':  # All values will be expressed as a percentage of the total signal
#             cols = y_order
#             rows = df_slice.loc[df_slice['well'].isin(group)][cols]
#             rows['total'] = rows.sum(axis=1)
#             rows = rows.div(rows['total'], axis=0)
#
#         elif mode == 'relative':    # All values will be expressed relative to 'norm_to'
#             cols = y_order + [norm_to]
#             rows = df_slice.loc[df_slice['well'].isin(group)][cols]
#             rows = rows.div(rows[norm_to], axis=0)
#
#         else:
#             raise ValueError("mode must be \'fraction\' or \'relative\'")
#
#         x = x_labels[i]
#         ddict[x] = {}
#         for moi in y_order:
#             ddict[x][moi] = {}
#             ddict[x][moi]['vals'] = list(rows[moi])
#             ddict[x][moi]['avg'] = np.average(ddict[x][moi]['vals'])
#             ddict[x][moi]['std'] = np.std(ddict[x][moi]['vals'])
#         i += 1
#
#     if sort_by:
#         if sort_order == 'ascending' or sort_order is None:
#             reverse = False
#         elif sort_order == 'descending':
#             reverse = True
#         else:
#             raise ValueError("sort_order must be \'ascending\', \'descending\', or \'None\'")
#
#         # First, sort dictionary alphabetically
#         ddict = {k: v for k, v in sorted(ddict.items(), key=lambda item: item[0])}
#
#         # Then, sort dictionary by amino acid position
#         ddict = {k: v for k, v in sorted(ddict.items(), key=lambda item: get_aa_position(item[0]))}
#
#         # check if it's sorting on a single MOI, if it is, put it in a list to standardize next step
#         if isinstance(sort_by, str):
#             sort_by = [sort_by]
#
#         # if it's multiple MOIs, reverse them so you sort by the lowest priority first
#         # this results in a final order that prioritizes the first item in sort_by
#         elif isinstance(sort_by, list):
#             sort_by.reverse()
#
#         # Finally, sort dictionary based on the avg value of MOI designated as 'sort_by'
#         for moi in sort_by:
#             ddict = {k: v for k, v in sorted(ddict.items(), key=lambda item: item[1][moi]['avg'], reverse=reverse)}
#
#         # Update labels to reflect the sorting
#         x_labels = [k for k in ddict.keys()]
#
#     # Generate plot
#     plt.rcParams.update({'font.family': 'Arial',
#                          'font.size': 13,
#                          'text.color': 'black',
#                          'axes.labelcolor': 'black'})
#
#     sns.set_style('ticks')  # Necessary to see minor ticks
#
#     fig, ax = plt.subplots(figsize=figsize)
#
#     # populate bar chart by calculating means
#     bottom = [0] * len(x_labels)    # starting point for the bar
#     for i, moi in enumerate(y_order):
#         if y_labels:
#             label = y_labels[i]
#         else:
#             label = moi
#         avg = [ddict[x][moi]['avg'] for x in x_labels]
#         std = [ddict[x][moi]['std'] for x in x_labels]
#
#         bar = {'x': x_labels,
#                'height': avg,
#                'bottom': bottom,
#                'label': label,
#                'color': y_colours[i],
#                'edgecolor': 'black',
#                'linewidth': 0.5}
#
#         if errorbars:
#             bar['yerr'] = std   # default error bars are 1 standard deviation
#             bar['capsize'] = 3
#
#         plt.bar(**bar)
#
#         bottom = [v1 + v2 for v1, v2 in zip(bottom, avg)]   # update bottom for next bar
#
#     # add individual data points as a scatter
#     scatter_df = pd.DataFrame(columns=['x', 'moi', 'val'])  # df with 1 row per point to plot
#     for x in x_labels:
#         bottom = 0
#         for moi in y_order:
#             vals = ddict[x][moi]['vals']
#             vals = [v + bottom for v in vals]
#             bottom = bottom + ddict[x][moi]['avg']
#             entry = pd.DataFrame({'x': x, 'moi': moi, 'val': vals})
#             scatter_df = pd.concat([scatter_df, entry], ignore_index=True)
#
#     sns.swarmplot(data=scatter_df, x='x', y='val',
#                   hue='moi', palette=y_colours,
#                   linewidth=0.5, size=3,
#                   edgecolor='black', legend=False, clip_on=False)
#
#     # ticks
#     if mode == 'fraction':
#         plt.ylim(0, 1)
#         plt.yticks(ticks=[0, 0.25, 0.5, 0.75, 1.0], labels=['0%', '25%', '50%', '75%', '100%'])
#         plt.yticks(ticks=np.linspace(0, 1, num=20, endpoint=False), minor=True)
#         ax.yaxis.get_ticklocs(minor=True)
#         ax.minorticks_on()
#         ax.xaxis.set_tick_params(which='minor', bottom=False)   # turn off x-axis minor ticks
#
#     plt.xticks(rotation=45, ha='right')
#     ax.tick_params(color='black', labelcolor='black')
#     for spine in ax.spines.values():
#         spine.set_edgecolor('black')
#     plt.grid(color='black', axis='y', linewidth=0.5)
#
#     #labels, legend, title
#     plt.xlabel('')
#     plt.ylabel('Signal intensity (au)')
#     if title:
#         plt.title(title, pad=20)
#
#     plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', reverse=True)
#     # plt.tight_layout()
#
#     return fig, ax

