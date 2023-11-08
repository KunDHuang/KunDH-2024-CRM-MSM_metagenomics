#!/usr/bin/env python


"""
NAME: step_curve_drawer.py
DESCRIPTION: This script is to analyze the co-prsense of multiple species in different categories,
             by drawing step curves.
AUTHOR: Kun D. Huang
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import argparse
import textwrap

def read_args(args):
    # This function is to parse arguments

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description = textwrap.dedent('''\
                                     This program is to do draw step curves to analyze co-presense of multiple species in different groups.
                                     '''),
                                    epilog = textwrap.dedent('''\
                                    examples:step_curve_drawer.py --abundance_table <abundance_table_w_md.tsv> --variable <variable_name> --species_number <nr_sps> --output <output.svg>
                                    '''))
    parser.add_argument('--abundance_table',
                        nargs = '?',
                        help = 'Input the MetaPhlAn4 abundance table which contains only a group of species one wants to analyze their co-presense state, with metadata being wedged.',
                        type = str,
                        default = None)

    parser.add_argument('--variable',
                        nargs = '?',
                        help = 'Specify the header of the variable in the metadata table you want to assess. For example, \
                        [Diet] variable columns has three categries - [vegan]/[Flexitarian]/[Omnivore].',
                        type = str,
                        default = None)

    parser.add_argument('--minimum_abundance',
                        nargs = '?',
                        help = 'Specify the minimum abundance used for determining presense. note: [0, 100] and [0.0] by default',
                        type = float,
                        default = 0.0)

    parser.add_argument('--species_number',
                        nargs = '?',
                        help = 'Specify the total number of multiple species in the analysis.',
                        type = int)


    parser.add_argument('--output',
                        nargs = '?',
                        help = 'Specify the output figure name.',
                        type = str,
                        default = None)
    parser.add_argument('--palette',
                        nargs = '?',
                        help = 'Input a tab-delimited mapping file where values are group names and keys are color codes.',
                        type = str,
                        default = None)

    return vars(parser.parse_args())

class PandasDealer:
    """
    This is an object for dealing pandas dataframe.
    """

    def __init__(self, df_):

        self.df_ = df_

    def read_csv(self):
        # Ths fucntion will read tab-delimitted file into a pandas dataframe.

        return pd.read_csv(self.df_, sep = '\t', index_col = False, low_memory=False)

    def rotate_df(self):
        # this function is to rotate the metaphlan-style table into tidy dataframe to ease searching work,

        df_ = self.read_csv()
        df_rows_lists = df_.values.tolist()
        rotated_df_dict = {df_.columns[0]: df_.columns[1:]}
        for i in df_rows_lists:
            rotated_df_dict[i[0]] = i[1:]

        rotated_df = pd.DataFrame.from_dict(rotated_df_dict)
        
        return rotated_df

class CopEstimator:

    def __init__(self, sub_df_md):
        self.sub_df_md = sub_df_md # sub_df_md: a subset of dataframe which contains only a group of species one wants to do co-presence analysis.

    def make_copresense_df(self, factor, total_species_nr, threshold = 0.0):
        # factor: the factor you want to assess the category percentage.
        # total_species_nr: specify the total number of species you want to do co-presense analysis.


        rotated_df = PandasDealer(self.sub_df_md)
        rotated_df = rotated_df.rotate_df()
        cols = rotated_df.columns[-total_species_nr: ].to_list() 
        categories = list(set(rotated_df[factor].to_list()))
        

        copresense = []
        cate_name = []
        ratios = []
        for c in categories:
            sub_df = rotated_df[rotated_df[factor] == c]
            species_group_df = sub_df[cols]
            species_group_df = species_group_df.apply(pd.to_numeric)
            species_group_df['total'] = species_group_df[cols].gt(threshold).sum(axis=1)
            for i in range(1, total_species_nr + 1):
                ratio = count_non_zero_rows(species_group_df, i)
                copresense.append(i)
                cate_name.append(c)
                ratios.append(ratio)

        return pd.DataFrame.from_dict({"copresense": copresense,
                                        factor: cate_name,
                                        "percentage": ratios})

def count_non_zero_rows(df_, nr):
    total_rows = len(df_.index)
    
    sub_df = df_[df_['total'] >= nr]
    ratio = len(sub_df.index)/total_rows

    return ratio
    

class VisualTools:
    def __init__(self, processed_df, factor):
        self.processed_df = processed_df
        self.factor = factor

    def step_curves(self, opt_name, palette = None):
        categories = list(set(self.processed_df[self.factor].to_list()))
        if palette:
            palette_dict = {i.rstrip().split('\t')[0]: i.rstrip().split('\t')[1] for i in open(palette).readlines()}
            for c in categories:
                sub_df = self.processed_df[self.processed_df[self.factor] == c]
                plt.step(sub_df["percentage"]*100, sub_df["copresense"], label = c, color = palette_dict[c])
        else:
            for c in categories:
                sub_df = self.processed_df[self.processed_df[self.factor] == c]
                plt.step(sub_df["percentage"]*100, sub_df["copresense"], label = c)

        plt.title("Number of species in an individual if present")
        plt.xlabel("Percentage")
        plt.ylabel("Co-presense")
        plt.legend(title = self.factor)
        plt.savefig(opt_name, bbox_inches = "tight")


if __name__ == "__main__":

    pars = read_args(sys.argv)
    cop_obj = CopEstimator(pars['abundance_table'])
    p_df = cop_obj.make_copresense_df(pars['variable'], pars['species_number'], pars['minimum_abundance'])
    vis_obj = VisualTools(p_df, pars['variable'])
    vis_obj.step_curves(pars['output'], palette = pars['palette'])
