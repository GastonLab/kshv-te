#!/usr/bin/env python

import csv
import sys
import math
import plotly
import argparse
import scipy.stats
import plotly.plotly as py
import plotly.graph_objs as go

from math import erf
from math import sqrt
from collections import defaultdict


def z_to_p(z):
    left_p = 0.5 * (1 + erf(z / sqrt(2)))
    p = 2 * (1 - left_p)

    return p


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--threshold', help="Threshold CPM count for removal of low count transcripts",
                        type=float, default=1.0)
    parser.add_argument('-i', '--input', help='Input file name with CPM/RPKM Values as a Matrix')
    parser.add_argument('-o', '--outroot', help="Root for output files")
    args = parser.parse_args()

    expression = defaultdict(lambda: defaultdict(dict))
    transcript2gene = dict()
    expression_vectors = defaultdict(list)

    sys.stdout.write("Reading input data from {} and calculating TE Values\n".format(args.input))
    with open(args.input, 'r') as infile:
        reader = csv.DictReader(infile, dialect='excel-tab')
        for row in reader:

            # Filter out any entry with less than 10 in any sample
            flag = 0
            if float(row['11_PD']) < args.threshold:
                flag = 1
            elif float(row['20_PD']) < args.threshold:
                flag = 1
            elif float(row['11_PT']) < args.threshold:
                flag = 1
            elif float(row['20_PT']) < args.threshold:
                flag = 1
            elif float(row['11_TD']) < args.threshold:
                flag = 1
            elif float(row['20_TD']) < args.threshold:
                flag = 1
            elif float(row['11_TT']) < args.threshold:
                flag = 1
            elif float(row['20_TT']) < args.threshold:
                flag = 1
            else:
                flag = 0

            if flag == 1:
                # sys.stderr.write("Found count less than 10 or transcript {}. Skipping...\n".format(row['Transcript']))
                continue

            if row['Transcript'].startswith("ENST"):
                temp = row['Transcript'].split('|')
                transcript = temp[0]
                transcript2gene[temp[0]] = temp[5]
            else:
                transcript = row['Transcript']
                transcript2gene[transcript] = transcript

            expression[transcript]['PD']['11'] = row['11_PD']
            expression[transcript]['PD']['20'] = row['20_PD']
            expression[transcript]['PT']['11'] = row['11_PT']
            expression[transcript]['PT']['20'] = row['20_PT']
            expression[transcript]['TD']['11'] = row['11_TD']
            expression[transcript]['TD']['20'] = row['20_TD']
            expression[transcript]['TT']['11'] = row['11_TT']
            expression[transcript]['TT']['20'] = row['20_TT']

            expression[transcript]['Control_TE']['11'] = math.log((float(row['11_PD']) / float(row['11_TD'])), 2)
            expression[transcript]['Control_TE']['20'] = math.log((float(row['20_PD']) / float(row['20_TD'])), 2)
            expression[transcript]['Treated_TE']['11'] = math.log((float(row['11_PT']) / float(row['11_TT'])), 2)
            expression[transcript]['Treated_TE']['20'] = math.log((float(row['20_PT']) / float(row['20_TT'])), 2)

            expression[transcript]['Raw_Control_TE']['11'] = float(row['11_PD']) / float(row['11_TD'])
            expression[transcript]['Raw_Control_TE']['20'] = float(row['20_PD']) / float(row['20_TD'])
            expression[transcript]['Raw_Treated_TE']['11'] = float(row['11_PT']) / float(row['11_TT'])
            expression[transcript]['Raw_Treated_TE']['20'] = float(row['20_PT']) / float(row['20_TT'])

            expression[transcript]['Control_TE']['Avg'] = (expression[transcript]['Control_TE']['11'] +
                                                           expression[transcript]['Control_TE']['20']) / 2
            expression[transcript]['Treated_TE']['Avg'] = (expression[transcript]['Treated_TE']['11'] +
                                                           expression[transcript]['Treated_TE']['20']) / 2

            expression[transcript]['Delta_TE']['11'] = expression[transcript]['Treated_TE']['11'] - \
                                                       expression[transcript]['Control_TE']['11']
            expression[transcript]['Delta_TE']['20'] = expression[transcript]['Treated_TE']['20'] - \
                                                       expression[transcript]['Control_TE']['20']

            expression[transcript]['TE_FC']['11'] = \
                math.log(expression[transcript]['Raw_Treated_TE']['11'] /
                         expression[transcript]['Raw_Control_TE']['11'], 2)
            expression[transcript]['TE_FC']['20'] = \
                math.log(expression[transcript]['Raw_Treated_TE']['20'] /
                         expression[transcript]['Raw_Control_TE']['20'], 2)

            expression[transcript]['TE_FC']['Avg'] = (expression[transcript]['TE_FC']['11'] +
                                                      expression[transcript]['TE_FC']['20']) / 2

            expression_vectors['11_PD'].append(row['11_PD'])
            expression_vectors['11_PT'].append(row['11_PT'])
            expression_vectors['11_TD'].append(row['11_TD'])
            expression_vectors['11_TT'].append(row['11_TT'])
            expression_vectors['20_PD'].append(row['20_PD'])
            expression_vectors['20_PT'].append(row['20_PT'])
            expression_vectors['20_TD'].append(row['20_TD'])
            expression_vectors['20_TT'].append(row['20_TT'])

            expression_vectors['11_control_te'].append(expression[transcript]['Control_TE']['11'])
            expression_vectors['20_control_te'].append(expression[transcript]['Control_TE']['20'])
            expression_vectors['11_treated_te'].append(expression[transcript]['Treated_TE']['11'])
            expression_vectors['20_treated_te'].append(expression[transcript]['Treated_TE']['20'])

            expression_vectors['11_TE_FC'].append(expression[transcript]['TE_FC']['11'])
            expression_vectors['20_TE_FC'].append(expression[transcript]['TE_FC']['20'])
            expression_vectors['Avg_TE_FC'].append(expression[transcript]['TE_FC']['Avg'])

    # Calculating Z-Scores
    zscore_11 = scipy.stats.zscore(expression_vectors['11_TE_FC'])
    zscore_20 = scipy.stats.zscore(expression_vectors['20_TE_FC'])
    zscore_avg = scipy.stats.zscore(expression_vectors['Avg_TE_FC'])

    num_tests = len(expression_vectors['11_TE_FC'])

    # Output TE Data
    sys.stdout.write("Outputting TE Data for {} transcripts\n".format(num_tests))
    with open("{}_TE_Data.txt".format(args.outroot), 'w') as te_outfile:
        te_outfile.write("Transcript\tGene\t"
                         "11_Control_Raw_TE\t20_Control_Raw_TE\t11_Treated_Raw_TE\t20_Treated_Raw_TE\t"
                         "11_Control_Log2TE\t20_Control_Log2TE\t11_Treated_Log2TE\t20_Treated_Log2TE\t"
                         "Avg_Control_TE\tAvg_Treated_TE\t"
                         "11_Delta_Log2TE\t20_Delta_Log2TE\t11_FC_Log2TE\t20_FC_Log2TE\tAvg_FC_Log2TE\t"
                         "11_FC_ZScore\t11_PValue\t11_Adj_PValue\t20_FC_ZScore\t20_PValue\t20_Adj_PValue\t"
                         "Avg_FC_ZScore\tAvg_PValue\tAvg_Adj_PValue\n")
        i = 0
        for transcript in expression.keys():
            p_11 = scipy.stats.norm.sf(abs(zscore_11[i])) * 2
            p_11_adj = p_11 * num_tests

            p_20 = scipy.stats.norm.sf(abs(zscore_20[i])) * 2
            p_20_adj = p_20 * num_tests

            p_avg = scipy.stats.norm.sf(abs(zscore_avg[i])) * 2
            p_avg_adj = p_avg * num_tests

            te_outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
                             "\t{}\t{}\t{}\t{}\t{}\t{}"
                             "\n".format(transcript, transcript2gene[transcript],
                                         expression[transcript]['Raw_Control_TE']['11'],
                                         expression[transcript]['Raw_Control_TE']['20'],
                                         expression[transcript]['Raw_Treated_TE']['11'],
                                         expression[transcript]['Raw_Treated_TE']['20'],
                                         expression[transcript]['Control_TE']['11'],
                                         expression[transcript]['Control_TE']['20'],
                                         expression[transcript]['Treated_TE']['11'],
                                         expression[transcript]['Treated_TE']['20'],
                                         expression[transcript]['Control_TE']['Avg'],
                                         expression[transcript]['Treated_TE']['Avg'],
                                         expression[transcript]['Delta_TE']['11'],
                                         expression[transcript]['Delta_TE']['20'],
                                         expression[transcript]['TE_FC']['11'],
                                         expression[transcript]['TE_FC']['20'],
                                         expression[transcript]['TE_FC']['Avg'],
                                         zscore_11[i],
                                         p_11,
                                         p_11_adj,
                                         zscore_20[i],
                                         p_20,
                                         p_20_adj,
                                         zscore_avg[i],
                                         p_avg,
                                         p_avg_adj
                                         ))
            i += 1

    # sys.stdout.write("Calculate correlation coeffients\n")
    #
    # # Calculate the Spearman's Correlation Coefficient between samples
    # pd_rho, pd_pvalue = scipy.stats.spearmanr(expression_vectors['11_PD'], expression_vectors['20_PD'])
    # pt_rho, pt_pvalue = scipy.stats.spearmanr(expression_vectors['11_PT'], expression_vectors['20_PT'])
    # td_rho, td_pvalue = scipy.stats.spearmanr(expression_vectors['11_TD'], expression_vectors['20_TD'])
    # tt_rho, tt_pvalue = scipy.stats.spearmanr(expression_vectors['11_TT'], expression_vectors['20_TT'])
    #
    # control_te_rho, control_te_pvalue = scipy.stats.spearmanr(expression_vectors['11_control_te'],
    #                                                           expression_vectors['20_control_te'])
    # treated_te_rho, treated_te_pvalue = scipy.stats.spearmanr(expression_vectors['11_treated_te'],
    #                                                           expression_vectors['20_treated_te'])
    #
    # rho_11, pvalue_11 = scipy.stats.spearmanr(expression_vectors['11_PD'], expression_vectors['11_TE_FC'])
    # rho_20, pvalue_20 = scipy.stats.spearmanr(expression_vectors['20_PD'], expression_vectors['20_TE_FC'])
    #
    # sys.stdout.write("Plotting\n")
    #
    # fc_abundance_11_trace = go.Scatter(
    #     x=expression_vectors['11_PD'],
    #     y=expression_vectors['11_TE_FC'],
    #     mode='markers'
    # )
    #
    # fc_abundance_20_trace = go.Scatter(
    #     x=expression_vectors['20_PD'],
    #     y=expression_vectors['20_TE_FC'],
    #     mode='markers'
    # )
    #
    # pd_trace = go.Scatter(
    #     x=expression_vectors['11_PD'],
    #     y=expression_vectors['20_PD'],
    #     mode='markers'
    # )
    #
    # pt_trace = go.Scatter(
    #     x=expression_vectors['11_PT'],
    #     y=expression_vectors['20_PT'],
    #     mode='markers'
    # )
    #
    # td_trace = go.Scatter(
    #     x=expression_vectors['11_TD'],
    #     y=expression_vectors['20_TD'],
    #     mode='markers'
    # )
    #
    # tt_trace = go.Scatter(
    #     x=expression_vectors['11_TT'],
    #     y=expression_vectors['20_TT'],
    #     mode='markers'
    # )
    #
    # control_te_trace = go.Scatter(
    #     x=expression_vectors['11_control_te'],
    #     y=expression_vectors['20_control_te'],
    #     mode='markers'
    # )
    #
    # treated_te_trace = go.Scatter(
    #     x=expression_vectors['11_treated_te'],
    #     y=expression_vectors['20_treated_te'],
    #     mode='markers'
    # )
    #
    # pd_data = [pd_trace]
    # pt_data = [pt_trace]
    # td_data = [td_trace]
    # tt_data = [tt_trace]
    #
    # control_te_data = [control_te_trace]
    # treated_te_data = [treated_te_trace]
    #
    # fc_abundance_11_data = [fc_abundance_11_trace]
    # fc_abundance_20_data = [fc_abundance_20_trace]
    #
    # fc_abundance_11_layout = dict(title="11_PD Abundance vs 11_TE_FC Spearman's Correlation: {}, P-value: "
    #                                     "{}".format(rho_11, pvalue_11),
    #                               yaxis=dict(zeroline=False),
    #                               xaxis=dict(zeroline=False)
    #                               )
    #
    # fc_abundance_20_layout = dict(title="20_PD Abundance vs 20_TE_FC Spearman's Correlation: {}, P-value: "
    #                                     "{}".format(rho_20, pvalue_20),
    #                               yaxis=dict(zeroline=False),
    #                               xaxis=dict(zeroline=False)
    #                               )
    #
    # pd_layout = dict(title="11_PD vs 20_PD Spearman's Correlation: {}, P-value: {}".format(pd_rho, pd_pvalue),
    #                  yaxis=dict(zeroline=False),
    #                  xaxis=dict(zeroline=False)
    #                  )
    #
    # pt_layout = dict(title="11_PT vs 20_PT Spearman's Correlation: {}, P-value: {}".format(pt_rho, pt_pvalue),
    #                  yaxis=dict(zeroline=False),
    #                  xaxis=dict(zeroline=False)
    #                  )
    # td_layout = dict(title="11_TD vs 20_TD Spearman's Correlation: {}, P-value: {}".format(td_rho, td_pvalue),
    #                  yaxis=dict(zeroline=False),
    #                  xaxis=dict(zeroline=False)
    #                  )
    #
    # tt_layout = dict(title="11_TT vs 20_TT Spearman's Correlation: {}, P-value: {}".format(tt_rho, tt_pvalue),
    #                  yaxis=dict(zeroline=False),
    #                  xaxis=dict(zeroline=False)
    #                  )
    #
    # control_te_layout = dict(title="11_ControlTE vs 20_ControlTE Spearman's "
    #                                "Correlation: {}, P-value: {}".format(control_te_rho, control_te_pvalue),
    #                          yaxis=dict(zeroline=False),
    #                          xaxis=dict(zeroline=False)
    #                          )
    #
    # treated_te_layout = dict(title="11_TreatedTE vs 20_TreatedTE Spearman's "
    #                                "Correlation: {}, P-value: {}".format(treated_te_rho, treated_te_pvalue),
    #                          yaxis=dict(zeroline=False),
    #                          xaxis=dict(zeroline=False)
    #                          )
    #
    # pd_fig = dict(data=pd_data, layout=pd_layout)
    # plotly.offline.plot(pd_fig, filename="{}_PD_correlation_plots".format(args.outroot))
    #
    # pt_fig = dict(data=pt_data, layout=pt_layout)
    # plotly.offline.plot(pt_fig, filename="{}_PT_correlation_plots".format(args.outroot))
    #
    # td_fig = dict(data=td_data, layout=td_layout)
    # plotly.offline.plot(td_fig, filename="{}_TD_correlation_plots".format(args.outroot))
    #
    # tt_fig = dict(data=tt_data, layout=tt_layout)
    # plotly.offline.plot(tt_fig, filename="{}_TT_correlation_plots".format(args.outroot))
    #
    # control_te_fig = dict(data=control_te_data, layout=control_te_layout)
    # plotly.offline.plot(control_te_fig, filename="{}_Control_TE_correlation_plots".format(args.outroot))
    #
    # treated_te_fig = dict(data=treated_te_data, layout=treated_te_layout)
    # plotly.offline.plot(treated_te_fig, filename="{}_Treated_TE_correlation_plots".format(args.outroot))
    #
    # fc_fig_11 = dict(data=fc_abundance_11_data, layout=fc_abundance_11_layout)
    # plotly.offline.plot(fc_fig_11, filename="{}_11_TE_FC_Abundance_Correlation".format(args.outroot))
    #
    # fc_fig_20 = dict(data=fc_abundance_20_data, layout=fc_abundance_20_layout)
    # plotly.offline.plot(fc_fig_20, filename="{}_20_TE_FC_Abundance_Correlation".format(args.outroot))
