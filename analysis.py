#!/usr/bin/env python

import csv
import sys
import math
import plotly
import argparse
import scipy.stats
import plotly.plotly as py
import plotly.graph_objs as go

from collections import defaultdict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input file name with CPM/RPKM Values as a Matrix')
    parser.add_argument('-o', '--outroot', help="Root for output files")
    args = parser.parse_args()

    expression = defaultdict(lambda: defaultdict(dict))
    expression_vectors = defaultdict(list)

    sys.stdout.write("Reading input data from {} and calculating TE Values\n".format(args.input))
    with open(args.input, 'r') as infile:
        reader = csv.DictReader(infile, dialect='excel-tab')
        for row in reader:

            # Filter out any entry with 0 in a sample
            flag = 0
            if float(row['11_PD']) == 0.0:
                flag = 1
            elif float(row['20_PD']) == 0.0:
                flag = 1
            elif float(row['11_PT']) == 0.0:
                flag = 1
            elif float(row['20_PT']) == 0.0:
                flag = 1
            elif float(row['11_TD']) == 0.0:
                flag = 1
            elif float(row['20_TD']) == 0.0:
                flag = 1
            elif float(row['11_TT']) == 0.0:
                flag = 1
            elif float(row['20_TT']) == 0.0:
                flag = 1

            if flag == 1:
                sys.stderr.write("Found zero values for transcript {}. Skipping...\n".format(row['Transcript']))
                continue

            expression[row['Transcript']]['PD']['11'] = row['11_PD']
            expression[row['Transcript']]['PD']['20'] = row['20_PD']
            expression[row['Transcript']]['PT']['11'] = row['11_PT']
            expression[row['Transcript']]['PT']['20'] = row['20_PT']
            expression[row['Transcript']]['TD']['11'] = row['11_TD']
            expression[row['Transcript']]['TD']['20'] = row['20_TD']
            expression[row['Transcript']]['TT']['11'] = row['11_TT']
            expression[row['Transcript']]['TT']['20'] = row['20_TT']

            expression[row['Transcript']]['Control_TE']['11'] = math.log((float(row['11_PD']) / float(row['11_TD'])), 2)
            expression[row['Transcript']]['Control_TE']['20'] = math.log((float(row['20_PD']) / float(row['20_TD'])), 2)
            expression[row['Transcript']]['Treated_TE']['11'] = math.log((float(row['11_PT']) / float(row['11_TT'])), 2)
            expression[row['Transcript']]['Treated_TE']['20'] = math.log((float(row['20_PT']) / float(row['20_TT'])), 2)

            expression[row['Transcript']]['Raw_Control_TE']['11'] = float(row['11_PD']) / float(row['11_TD'])
            expression[row['Transcript']]['Raw_Control_TE']['20'] = float(row['20_PD']) / float(row['20_TD'])
            expression[row['Transcript']]['Raw_Treated_TE']['11'] = float(row['11_PT']) / float(row['11_TT'])
            expression[row['Transcript']]['Raw_Treated_TE']['20'] = float(row['20_PT']) / float(row['20_TT'])

            expression[row['Transcript']]['Delta_TE']['11'] = expression[row['Transcript']]['Treated_TE']['11'] - \
                                                              expression[row['Transcript']]['Control_TE']['11']
            expression[row['Transcript']]['Delta_TE']['20'] = expression[row['Transcript']]['Treated_TE']['20'] - \
                                                              expression[row['Transcript']]['Control_TE']['20']

            expression[row['Transcript']]['TE_FC']['11'] = \
                math.log(expression[row['Transcript']]['Raw_Treated_TE']['11'] /
                         expression[row['Transcript']]['Raw_Control_TE']['11'], 2)
            expression[row['Transcript']]['TE_FC']['20'] = \
                math.log(expression[row['Transcript']]['Raw_Treated_TE']['20'] /
                         expression[row['Transcript']]['Raw_Control_TE']['20'], 2)

            expression_vectors['11_PD'].append(row['11_PD'])
            expression_vectors['11_PT'].append(row['11_PT'])
            expression_vectors['11_TD'].append(row['11_TD'])
            expression_vectors['11_TT'].append(row['11_TT'])
            expression_vectors['20_PD'].append(row['20_PD'])
            expression_vectors['20_PT'].append(row['20_PT'])
            expression_vectors['20_TD'].append(row['20_TD'])
            expression_vectors['20_TT'].append(row['20_TT'])

            expression_vectors['11_control_te'].append(expression[row['Transcript']]['Control_TE']['11'])
            expression_vectors['20_control_te'].append(expression[row['Transcript']]['Control_TE']['20'])
            expression_vectors['11_treated_te'].append(expression[row['Transcript']]['Treated_TE']['11'])
            expression_vectors['20_treated_te'].append(expression[row['Transcript']]['Treated_TE']['20'])

    # Output TE Data
    sys.stdout.write("Outputting TE Data\n")
    with open("{}_TE_Data.txt".format(args.outroot), 'w') as te_outfile:
        te_outfile.write("Transcript\t11_Control_Raw_TE\t20_Control_Raw_TE\t11_Treated_Raw_TE\t20_Treated_Raw_TE\t"
                         "11_Control_Log2TE\t20_Control_Log2TE\t11_Treated_Log2TE\t20_Treated_Log2TE\t"
                         "11_Delta_Log2TE\t20_Delta_Log2TE\t11_FC_Log2TE\t20_FC_Log2TE\n")
        for transcript in expression.keys():
            te_outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
                             "\n".format(transcript,
                                         expression[transcript]['Raw_Control_TE']['11'],
                                         expression[transcript]['Raw_Control_TE']['20'],
                                         expression[transcript]['Raw_Treated_TE']['11'],
                                         expression[transcript]['Raw_Treated_TE']['20'],
                                         expression[transcript]['Control_TE']['11'],
                                         expression[transcript]['Control_TE']['20'],
                                         expression[transcript]['Treated_TE']['11'],
                                         expression[transcript]['Treated_TE']['20'],
                                         expression[transcript]['Delta_TE']['11'],
                                         expression[transcript]['Delta_TE']['20'],
                                         expression[transcript]['TE_FC']['11'],
                                         expression[transcript]['TE_FC']['20']
                                         ))

    sys.stdout.write("Calculate correlation coeffients\n")
    # Calculate the Spearman's Correlation Coefficient between samples
    pd_rho, pd_pvalue = scipy.stats.spearmanr(expression_vectors['11_PD'], expression_vectors['20_PD'])
    pt_rho, pt_pvalue = scipy.stats.spearmanr(expression_vectors['11_PT'], expression_vectors['20_PT'])
    td_rho, td_pvalue = scipy.stats.spearmanr(expression_vectors['11_TD'], expression_vectors['20_TD'])
    tt_rho, tt_pvalue = scipy.stats.spearmanr(expression_vectors['11_TT'], expression_vectors['20_TT'])

    control_te_rho, control_te_pvalue = scipy.stats.spearmanr(expression_vectors['11_control_te'],
                                                              expression_vectors['20_control_te'])
    treated_te_rho, treated_te_pvalue = scipy.stats.spearmanr(expression_vectors['11_treated_te'],
                                                              expression_vectors['20_treated_te'])

    # # Calculate the Pearson's Correlation Coefficient between samples
    # pd_r, pd_ppvalue = scipy.stats.pearsonr(expression_vectors['11_PD'], expression_vectors['20_PD'])
    # pt_r, pt_ppvalue = scipy.stats.pearsonr(expression_vectors['11_PT'], expression_vectors['20_PT'])
    # td_r, td_ppvalue = scipy.stats.pearsonr(expression_vectors['11_TD'], expression_vectors['20_TD'])
    # tt_r, tt_ppvalue = scipy.stats.pearsonr(expression_vectors['11_TT'], expression_vectors['20_TT'])

    sys.stdout.write("Plotting\n")
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

    control_te_trace = go.Scatter(
        x=expression_vectors['11_control_te'],
        y=expression_vectors['20_control_te'],
        mode='markers'
    )

    treated_te_trace = go.Scatter(
        x=expression_vectors['11_treated_te'],
        y=expression_vectors['20_treated_te'],
        mode='markers'
    )

    pd_data = [pd_trace]
    pt_data = [pt_trace]
    td_data = [td_trace]
    tt_data = [tt_trace]

    control_te_data = [control_te_trace]
    treated_te_data = [treated_te_trace]

    pd_layout = dict(title="11_PD vs 20_PD Spearman's Correlation: {}, P-value: {}".format(pd_rho, pd_pvalue),
                     yaxis=dict(zeroline=False),
                     xaxis=dict(zeroline=False)
                     )

    pt_layout = dict(title="11_PT vs 20_PT Spearman's Correlation: {}, P-value: {}".format(pt_rho, pt_pvalue),
                     yaxis=dict(zeroline=False),
                     xaxis=dict(zeroline=False)
                     )
    td_layout = dict(title="11_TD vs 20_TD Spearman's Correlation: {}, P-value: {}".format(td_rho, td_pvalue),
                     yaxis=dict(zeroline=False),
                     xaxis=dict(zeroline=False)
                     )

    tt_layout = dict(title="11_TT vs 20_TT Spearman's Correlation: {}, P-value: {}".format(tt_rho, tt_pvalue),
                     yaxis=dict(zeroline=False),
                     xaxis=dict(zeroline=False)
                     )

    control_te_layout = dict(title="11_ControlTE vs 20_ControlTE Spearman's "
                                   "Correlation: {}, P-value: {}".format(control_te_rho, control_te_pvalue),
                             yaxis=dict(zeroline=False),
                             xaxis=dict(zeroline=False)
                             )

    treated_te_layout = dict(title="11_TreatedTE vs 20_TreatedTE Spearman's "
                                   "Correlation: {}, P-value: {}".format(treated_te_rho, treated_te_pvalue),
                             yaxis=dict(zeroline=False),
                             xaxis=dict(zeroline=False)
                             )

    pd_fig = dict(data=pd_data, layout=pd_layout)
    plotly.offline.plot(pd_fig, filename="{}_PD_correlation_plots".format(args.outroot))

    pt_fig = dict(data=pt_data, layout=pt_layout)
    plotly.offline.plot(pt_fig, filename="{}_PT_correlation_plots".format(args.outroot))

    td_fig = dict(data=td_data, layout=td_layout)
    plotly.offline.plot(td_fig, filename="{}_TD_correlation_plots".format(args.outroot))

    tt_fig = dict(data=tt_data, layout=tt_layout)
    plotly.offline.plot(tt_fig, filename="{}_TT_correlation_plots".format(args.outroot))

    control_te_fig = dict(data=control_te_data, layout=control_te_layout)
    plotly.offline.plot(control_te_fig, filename="{}_Control_TE_correlation_plots".format(args.outroot))

    treated_te_fig = dict(data=treated_te_data, layout=treated_te_layout)
    plotly.offline.plot(treated_te_fig, filename="{}_Treated_TE_correlation_plots".format(args.outroot))
