import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import warnings
import sys
import os
from collections.abc import Callable
import xlsxwriter

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
                 exclude_mois: list = None,
                 collapse_mois: dict | None = None,
                 activity_scoring_function: Callable = None,
                 activity_scoring_df_variables: dict = None,
                 activity_scoring_other_variables: dict = None):

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
        self.attributes = attributes
        self.excluded_mois = exclude_mois

        self.activity_scoring_function = activity_scoring_function
        self.activity_scoring_df_variables = activity_scoring_df_variables
        self.activity_scoring_other_variables = activity_scoring_other_variables

        # process data
        self.data_ungrouped = self.collapse_input_dataframes()
        self.data_mois_collapsed = self.collapse_mois(df=self.data_ungrouped, mapping_dict=collapse_mois)
        self.data_replicates_averaged = self.collapse_replicates(self.data_mois_collapsed, attributes)
        self.data_as_percentage = self.convert_to_percentage(self.data_replicates_averaged)

        # add activity scores
        if activity_scoring_function:
            self.data_as_percentage = self.add_activity_scores()
            self.final_activity_sums = self.sum_activities()
            self.evolvepro_formatted = self.format_for_evolvepro()

        self.dfs = {'ungrouped': self.data_ungrouped,
                    'mois_collapsed': self.data_mois_collapsed,
                    'replicates_averaged': self.data_replicates_averaged,
                    'percentage': self.data_as_percentage,
                    'activity_sums': self.final_activity_sums,
                    'evolvepro': self.evolvepro_formatted}

    def collapse_input_dataframes(self) -> pd.DataFrame:
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

        new_df.drop(columns=self.excluded_mois, inplace=True, errors='raise')

        return new_df

    def collapse_replicates(self, df: pd.DataFrame, attributes: list) -> pd.DataFrame:
        """Averages technical replicates into average values.
        Requires a list of attributes that should be shared across all replicates."""
        try:
            new_df = df.groupby(attributes).mean(numeric_only=True)
            new_df.reset_index(inplace=True)
            return new_df
        except KeyError:
            print(f"Failed attempting to group dataframe by attributes. Double check that each attribute matches a column header.")

    def convert_to_percentage(self, df: pd.DataFrame) -> pd.DataFrame:
        """Converts MOIs in each row to a percentage of total rather than absolute values."""

        labels_to_remove = self.attributes + ['chip', 'background', 'mz_offset', 'noise', 'noise_cutoff']

        #this df contains only numerical values
        df_for_summing = df.drop(labels=labels_to_remove, axis='columns', inplace=False, errors='raise')
        row_totals = df_for_summing.sum(numeric_only=True, axis='columns')
        df_percentage = df_for_summing.div(row_totals, axis='rows')

        # combine the categorical labels of the original dataframe with the percentage data of the new dataframe
        df_final = df[self.attributes].join(df_percentage)

        return df_final

    def add_activity_scores(self) -> pd.DataFrame:
        """
        Adds a column to the dataframe self.data_as_percentage which is the score for that row.
        Uses    self.activity_scoring_function,
                self.activity_scoring_df_variables,
                self.activity_scoring_other_variables
        """

        new_df = self.data_as_percentage.copy()
        scores = []
        for index, row in new_df.iterrows():
            # collect values from df
            df_dict = {}
            for k, v in self.activity_scoring_df_variables.items():
                df_dict[k] = row[v]
            score = self.activity_scoring_function(**df_dict, **self.activity_scoring_other_variables)
            scores.append(score)
        new_df['activity_score'] = scores

        return new_df

    def sum_activities(self) -> pd.DataFrame:
        new_df = self.data_as_percentage.copy()[['enz_id','enz_mutation','activity_score']]
        df_summed = new_df.groupby(['enz_id', 'enz_mutation'])['activity_score'].sum().reset_index()
        return df_summed

    def format_for_evolvepro(self) -> pd.DataFrame:
        old_df = self.final_activity_sums.copy()
        old_df.loc[old_df['enz_mutation'] == 'WT', 'enz_mutation'] = 'A82S'
        old_df.loc[old_df['enz_mutation'] == 'Base (S82A)', 'enz_mutation'] = 'WT'
        rows_to_drop = old_df[old_df['enz_mutation'] == 'No enzyme'].index
        old_df = old_df.drop(rows_to_drop)
        old_df.reset_index(inplace=True)

        new_df = pd.DataFrame()
        new_df['Variant'] = old_df['enz_mutation'].str[1:]
        new_df.loc[new_df['Variant'] == 'T', 'Variant'] = 'WT'  # fixes the WT from getting messed up from the above operation
        new_df['activity'] = old_df['activity_score']
        new_df.sort_values(by=['activity'], inplace=True, ascending=False)

        return new_df

    def sort_sample_labels(self,
                           df_name: str,
                           label_category: str,
                           sort_by: list | str,
                           global_var: tuple | None = None,
                           ascending: bool=False) -> list:
        """General purpose function for returning sample labels in a specific order."""

        df = self.dfs[df_name].copy()

        if global_var:
            df = df.loc[df[global_var[0]] == global_var[1]]  # take the slice that corresponds to global_var

        df = df.sort_values(by=sort_by, ascending=ascending, inplace=False)
        return list(df[label_category])

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

            rows = df.loc[df['source'].isin(source)]  # slice only the correct source plate(s)
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
                    hide_points: bool = False,
                    figsize=None,
                    title=None,
                    show_legend=True,
                    ax=None):   # can call this command to generate a standalone figure (ax=None) or to build an ax object inside an existing figure

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
                source = df.loc[df[x_category] == x_label]['source'].unique()   # this can be 1 or more plates
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
            fig, ax = plt.subplots(figsize=figsize)

        ax = self.generate_stacked_bar_ax(ddict=ddict,
                                          y_order=y_species if not y_legend_overwrite else y_legend_overwrite,
                                          y_colours=y_colours,
                                          figsize=figsize,
                                          title=title,
                                          show_legend=show_legend,
                                          hide_points=hide_points,
                                          ax=ax)
        return ax

    def generate_bar(self, ax: plt.axes, y_order: list, y_colours: list, ddict: dict):

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

    def generate_scatter(self, ax: plt.axes, y_order: list, y_colours: list, ddict: dict):
        # add individual data points as a scatter
        x_labels = ddict.keys()
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
                      linewidth=0.5, size=2.5,
                      edgecolor='black', legend=False, clip_on=False, ax=ax)

    def generate_stacked_bar_ax(self,
                                # x_labels: list,
                                y_order: list,
                                y_colours: list,
                                ddict: dict,
                                figsize: tuple,
                                title: str,
                                show_legend: bool,
                                hide_points: bool = False,
                                ax=None) -> plt.axes:

        """
        Helper function for self.stacked_bar(). This function skips over the data prep, and generates the actual plot.
        For param descriptions see self.stacked_bar()
        """
        self.generate_bar(ax=ax, y_order=y_order, y_colours=y_colours, ddict=ddict)

        if not hide_points:
            self.generate_scatter(ax=ax, y_order=y_order, y_colours=y_colours, ddict=ddict)
        else:   # necessary because seaborn and vanilla matplotlib have different spacing defaults
            ax.set_xlim(-0.5, len(ddict) - 0.5)

        # Formatting
        plt.rcParams.update({'font.family': 'Arial',
                             'font.size': 13,
                             'text.color': 'black',
                             'axes.labelcolor': 'black'})
        sns.set_style('ticks')  # Necessary to see minor ticks

        # ticks
        x_labels = ddict.keys()
        ax.set_ylim(0, 1)
        ax.set_yticks(ticks=[0, 0.25, 0.5, 0.75, 1.0], labels=['0%', '25%', '50%', '75%', '100%'])
        ax.set_yticks(ticks=np.linspace(0, 1, num=20, endpoint=False), minor=True)
        ax.yaxis.get_ticklocs(minor=True)
        ax.minorticks_on()
        ax.xaxis.set_tick_params(which='minor', bottom=False)  # turn off x-axis minor ticks
        ax.set_xticklabels(labels=x_labels, ha='right', rotation=45)
        ax.tick_params(color='black', labelcolor='black')

        # exterior and grid
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

        return ax

    def score_bar(self,
                  x_category,
                  x_labels: list | None = None, # determines both the species to include and their order
                  y_species: list = ('ntp', ['G','C','A','T']),
                  y_colours: list = ["#82cce7", "#e5b5bf", "#aadab4", "#ddcb9c"],
                  title='Activity score',
                  ax=None,
                  show_legend: bool = True):

        """Function only holds for my specific experimental setup because it assumes the score is a certain composite of
        different nucleotide additions"""
        plt.rcParams.update({'font.family': 'Arial',
                             'font.size': 13,
                             'text.color': 'black',
                             'axes.labelcolor': 'black'})

        df = self.data_as_percentage

        if x_labels is None:
            x_labels = list(df[x_category].unique())

        if ax is None:
            fig, ax = plt.subplots()

        ddict = {}
        for label in x_labels:
            ddict[label] = {}
            df_slice = df[df[x_category] == label]  # slice to correspond to x_label

            for y in y_species[1]:
                df_slice_slice = df_slice[df[y_species[0]] == y]    # slice to correspond to correct nucleotide
                ddict[label][y] = float(df_slice_slice['activity_score'])

        bottom = [0] * len(x_labels)  # starting point for the bar
        for i, y in enumerate(y_species[1]):
            height = [ddict[x][y] for x in x_labels]

            bar = {'x': x_labels,
                   'height': height,
                   'bottom': bottom,
                   'label': y,
                   'color': y_colours[i],
                   'edgecolor': 'black',
                   'linewidth': 0.5}

            ax.bar(**bar)
            for i, h in enumerate(height):
                if h < 0.3:
                    continue
                else:
                    ax.text(x=i, y=h/2 + bottom[i], s=y, ha='center', va='center', fontsize=10, color='black')
            bottom = [v1 + v2 for v1, v2 in zip(bottom, height)]  # update bottom for next bar

        # ticks
        ax.set_xlim(-0.5, len(x_labels)-0.5)    # for some reason this is necessary to get the same alignment as self.stacked_bar
        # ax.set_ylim(0, 4)
        # ax.set_yticks(ticks=[0, 1, 2, 3, 4])
        # ax.set_yticks(ticks=np.linspace(0, 4, num=20, endpoint=False), minor=True)
        ax.yaxis.get_ticklocs(minor=True)
        ax.minorticks_on()
        ax.xaxis.set_tick_params(which='minor', bottom=False)  # turn off x-axis minor ticks
        ax.set_xticklabels(labels=x_labels, ha='right', rotation=45)
        ax.tick_params(color='black', labelcolor='black')

        # exterior and grid
        for spine in ax.spines.values():
            spine.set_edgecolor('black')
        ax.grid(color='black', axis='y', linewidth=0.5)

        # labels, legend, title
        ax.set_xlabel('')
        ax.set_ylabel('Activity score')
        if title:
            ax.set_title(title, pad=20)
        if show_legend:
            ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')

        return ax

    def write_to_excel(self, filename: str, overwrite=False) -> None:
        """
        Writes all dataframes to an excel file.
        """

        # create file
        outfile = filename if filename.endswith(".xlsx") else filename + ".xlsx"

        if not overwrite:
            if os.path.exists(outfile):
                print(f"\'write_to_excel()\' did not execute because the file already exists.\n"
                      f"To proceed anyway, use parameter \'overwrite=True\'.")
                return

        writer = pd.ExcelWriter(path=outfile, engine='xlsxwriter')

        for k, v in self.dfs.items():
            v.to_excel(writer, sheet_name=k, index=False)

        writer.close()
