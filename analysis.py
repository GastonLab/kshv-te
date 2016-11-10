import csv
import sys
import argparse
import scipy.stats
import plotly.plotly as py
import plotly.graph_objs as go

from collections import defaultdict
from plotly.tools import FigureFactory as FF


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input file name with CPM/RPKM Values as a Matrix')
    parser.add_argument('-o', '--outroot', help="Root for output files")
    args = parser.parse_args()

    expression = defaultdict(lambda: defaultdict(dict))
    expression_vectors = defaultdict(list)

    with open(args.input, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            expression[row['Transcript']]['PD']['11'] = row['11_PD']
            expression[row['Transcript']]['PD']['20'] = row['20_PD']
            expression[row['Transcript']]['PT']['11'] = row['11_PT']
            expression[row['Transcript']]['PT']['20'] = row['20_PT']
            expression[row['Transcript']]['TD']['11'] = row['11_TD']
            expression[row['Transcript']]['TD']['20'] = row['20_TD']
            expression[row['Transcript']]['TT']['11'] = row['11_TT']
            expression[row['Transcript']]['TT']['20'] = row['20_TT']

            expression_vectors['11_PD'].append(row['11_PD'])
            expression_vectors['11_PT'].append(row['11_PT'])
            expression_vectors['11_TD'].append(row['11_TD'])
            expression_vectors['11_TT'].append(row['11_TT'])
            expression_vectors['20_PD'].append(row['20_PD'])
            expression_vectors['20_PT'].append(row['20_PT'])
            expression_vectors['20_TD'].append(row['20_TD'])
            expression_vectors['20_TT'].append(row['20_TT'])

    # Calculate the Spearman's Correlation Coefficient between samples
    pd_rho, pd_pvalue = scipy.stats.spearmanr(expression_vectors['11_PD'], expression_vectors['20_PD'])
    pt_rho, pt_pvalue = scipy.stats.spearmanr(expression_vectors['11_PT'], expression_vectors['20_PT'])
    td_rho, td_pvalue = scipy.stats.spearmanr(expression_vectors['11_TD'], expression_vectors['20_TD'])
    tt_rho, tt_pvalue = scipy.stats.spearmanr(expression_vectors['11_TT'], expression_vectors['20_TT'])

    # Calculate the Pearson's Correlation Coefficient between samples
    pd_r, pd_ppvalue = scipy.stats.pearsonr(expression_vectors['11_PD'], expression_vectors['20_PD'])
    pt_r, pt_ppvalue = scipy.stats.pearsonr(expression_vectors['11_PT'], expression_vectors['20_PT'])
    td_r, td_ppvalue = scipy.stats.pearsonr(expression_vectors['11_TD'], expression_vectors['20_TD'])
    tt_r, tt_ppvalue = scipy.stats.pearsonr(expression_vectors['11_TT'], expression_vectors['20_TT'])

    pd_trace = go.Scatter(
        x=expression_vectors['11_PD'],
        y=expression_vectors['20_PD'],
        mode='markers'
    )

    pt_trace = go.Scatter(
        x=expression_vectors['11_PT'],
        y=expression_vectors['20_PT'],
        mode='markers'
    )

    td_trace = go.Scatter(
        x=expression_vectors['11_TD'],
        y=expression_vectors['20_TD'],
        mode='markers'
    )

    tt_trace = go.Scatter(
        x=expression_vectors['11_TT'],
        y=expression_vectors['20_TT'],
        mode='markers'
    )

    pd_data = [pd_trace]
    pt_data = [pt_trace]
    td_data = [td_trace]
    tt_data = [tt_trace]

    pd_layout = dict(title='11_PD vs 20_PD',
                     yaxis=dict(zeroline=False),
                     xaxis=dict(zeroline=False)
                     )

    pt_layout = dict(title='11_PT vs 20_PT',
                     yaxis=dict(zeroline=False),
                     xaxis=dict(zeroline=False)
                     )
    td_layout = dict(title='11_TD vs 20_TD',
                     yaxis=dict(zeroline=False),
                     xaxis=dict(zeroline=False)
                     )

    tt_layout = dict(title='11_TT vs 20_TT',
                     yaxis=dict(zeroline=False),
                     xaxis=dict(zeroline=False)
                     )

    pd_fig = dict(data=pd_data, layout=pd_layout)
    py.iplot(pd_fig, filename="{}_PD_correlation_plots".format(args.outroot))

    pt_fig = dict(data=pt_data, layout=pt_layout)
    py.iplot(pd_fig, filename="{}_PD_correlation_plots".format(args.outroot))

    td_fig = dict(data=td_data, layout=td_layout)
    py.iplot(pd_fig, filename="{}_PD_correlation_plots".format(args.outroot))

    tt_fig = dict(data=tt_data, layout=tt_layout)
    py.iplot(pd_fig, filename="{}_PD_correlation_plots".format(args.outroot))
