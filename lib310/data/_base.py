from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd
from ._utils import color


class ProteinDataTable(pd.DataFrame):
    @classmethod
    def from_pandas(cls, df):
        return cls(data=df.copy())

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def summary(self, verbosity: int = 1):
        if verbosity == 1:
            print(f'* Shape: ({color.BOLD}{self.shape[0]}{color.END}, {color.BOLD}{self.shape[1]}{color.END})\n')
            print(f'* Columns information:')
            for i, column in enumerate(self.columns):
                print(
                    f'\t{i + 1}. {color.BOLD}{column}{color.END}: data type {color.BLUE}{self.dtypes[column]}{color.END} with {color.BOLD}{self[column].isna().sum()}{color.END} Missing values')

            print(f'\n* Total number of missing values in the data: {color.BOLD}{self.isna().sum().sum()}{color.END}')
            print(f'* Number of rows with missing values: {color.BOLD}{self.shape[0] - self.dropna().shape[0]}{color.END}')
        else:
            raise ValueError('Wrong verbosity level has been set! You have to set it only to `1` for now!')

    def plot(self, column: str, plot_type: Optional[str] = None):
        assert column in self.columns

        import seaborn as sns
        sns.set(rc={"figure.dpi": 150, 'savefig.dpi': 150})
        sns.set_context('notebook')
        sns.set_style("ticks")

        dtype = str(self.dtypes[column]).lower()

        if dtype in ['int64', 'int32', 'int']:
            if plot_type is None or plot_type == 'hist':
                sns.displot(
                    self, x=column,
                    binwidth=100, height=3, aspect=2, facet_kws=dict(margin_titles=True), kind='hist',
                )
                plt.show()
            elif plot_type in ['dist', 'kde']:
                sns.displot(
                    self, x=column,
                    binwidth=100, height=3, aspect=2, facet_kws=dict(margin_titles=True), kind='kde',
                )
                plt.show()
        elif dtype in ['object', 'str'] and plot_type == 'population':
            self['temp'] = self[column].apply(lambda cell: ''.join(c for c in cell if c not in "\'\"[]").split(', '))

            all_terms = set()
            for term_group in self['temp'].to_list():
                all_terms |= set(term_group)
            all_terms = sorted(list(all_terms))

            print(f'Found {len(all_terms)} different values')

            term_population = {column: [], 'population': []}
            for term in all_terms:
                term_df = self.loc[self['temp'].apply(lambda x: term in x)]

                term_population[column].append(term)
                term_population['population'].append(len(term_df))

            plt.figure(figsize=(15, 6))
            sns.barplot(x=column, y='population', palette="deep", data=pd.DataFrame.from_dict(term_population))
            plt.xticks(rotation=90)
            plt.grid(axis='y')
            plt.show()
        else:
            print('Invalid function usage')
