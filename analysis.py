#!/usr/bin/env python

import csv
import sys
import math
import plotly
import argparse
import scipy.stats
import plotly.plotly as py
import plotly.graph_objs as go
import statsmodels.api as sm

from collections import defaultdict


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

    transcript_labels = list()
    human_transcript_labels = list()
    kshv_transcript_labels = list()

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

            transcript_labels.append(transcript)

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

            expression_vectors['log_11_PD'].append(math.log(float(row['11_PD']), 2))
            expression_vectors['log_11_PT'].append(math.log(float(row['11_PT']), 2))
            expression_vectors['log_11_TD'].append(math.log(float(row['11_TD']), 2))
            expression_vectors['log_11_TT'].append(math.log(float(row['11_TT']), 2))
            expression_vectors['log_20_PD'].append(math.log(float(row['20_PD']), 2))
            expression_vectors['log_20_PT'].append(math.log(float(row['20_PT']), 2))
            expression_vectors['log_20_TD'].append(math.log(float(row['20_TD']), 2))
            expression_vectors['log_20_TT'].append(math.log(float(row['20_TT']), 2))

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

            # Set up Human versus KSHV expression vectors
            if transcript.startswith("trans"):
                kshv_transcript_labels.append(transcript)
                expression_vectors['kshv_11_TE_FC'].append(expression[transcript]['TE_FC']['11'])
                expression_vectors['kshv_20_TE_FC'].append(expression[transcript]['TE_FC']['20'])
                expression_vectors['kshv_Avg_TE_FC'].append(expression[transcript]['TE_FC']['Avg'])
            elif transcript.startswith("ENST"):
                human_transcript_labels.append(transcript)
                expression_vectors["human_11_TE_FC"].append(expression[transcript]['TE_FC']['11'])
                expression_vectors["human_20_TE_FC"].append(expression[transcript]['TE_FC']['20'])
                expression_vectors["human_Avg_TE_FC"].append(expression[transcript]['TE_FC']['Avg'])
            else:
                sys.stderr.write("Could not determine organism for transcript {}\n".format(transcript))

    # Calculating Z-Scores
    zscore_11 = scipy.stats.zscore(expression_vectors['11_TE_FC'])
    zscore_20 = scipy.stats.zscore(expression_vectors['20_TE_FC'])
    zscore_avg = scipy.stats.zscore(expression_vectors['Avg_TE_FC'])

    # Calculate P-Values and Adjusted P_Values
    p_11 = list()
    for zscore in zscore_11:
        p_11.append(scipy.stats.norm.sf(abs(zscore)) * 2)
    p_11_adj_results = sm.stats.multipletests(p_11, method='b')
    p_11_adj = p_11_adj_results[1]

    p_20 = list()
    for zscore in zscore_20:
        p_20.append(scipy.stats.norm.sf(abs(zscore)) * 2)
    p_20_adj_results = sm.stats.multipletests(p_20, method='b')
    p_20_adj = p_20_adj_results[1]

    p_avg = list()
    for zscore in zscore_avg:
        p_avg.append(scipy.stats.norm.sf(abs(zscore)) * 2)
    p_avg_adj_results = sm.stats.multipletests(p_avg, method='b')
    p_avg_adj = p_avg_adj_results[1]

    num_tests = len(expression_vectors['Avg_TE_FC'])

    # Calculate Mann-Whitney Rank Tests on values

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
                                         p_11[i],
                                         p_11_adj[i],
                                         zscore_20[i],
                                         p_20[i],
                                         p_20_adj[i],
                                         zscore_avg[i],
                                         p_avg[i],
                                         p_avg_adj[i]
                                         ))
            i += 1

    sys.stdout.write("Calculate correlation coefficients\n")

    # Calculate the Spearman's Correlation Coefficient between samples
    pd_rho, pd_pvalue = scipy.stats.spearmanr(expression_vectors['11_PD'], expression_vectors['20_PD'])
    pt_rho, pt_pvalue = scipy.stats.spearmanr(expression_vectors['11_PT'], expression_vectors['20_PT'])
    td_rho, td_pvalue = scipy.stats.spearmanr(expression_vectors['11_TD'], expression_vectors['20_TD'])
    tt_rho, tt_pvalue = scipy.stats.spearmanr(expression_vectors['11_TT'], expression_vectors['20_TT'])

    control_te_rho, control_te_pvalue = scipy.stats.spearmanr(expression_vectors['11_control_te'],
                                                              expression_vectors['20_control_te'])
    treated_te_rho, treated_te_pvalue = scipy.stats.spearmanr(expression_vectors['11_treated_te'],
                                                              expression_vectors['20_treated_te'])

    rho_11, pvalue_11 = scipy.stats.mannwhitneyu(expression_vectors['11_PD'], expression_vectors['11_TE_FC'],
                                                 alternative='two-sided')
    rho_20, pvalue_20 = scipy.stats.mannwhitneyu(expression_vectors['20_PD'], expression_vectors['20_TE_FC'],
                                                 alternative='two-sided')

    rho_11_control_te_abundance, pvalue_11_control_te_abundance = scipy.stats.mannwhitneyu(
        expression_vectors['11_control_te'], expression_vectors['11_TD'], alternative='two-sided')

    rho_11_treated_te_abundance, pvalue_11_treated_te_abundance = scipy.stats.mannwhitneyu(
        expression_vectors['11_treated_te'], expression_vectors['11_TT'], alternative='two-sided')

    rho_20_control_te_abundance, pvalue_20_control_te_abundance = scipy.stats.mannwhitneyu(
        expression_vectors['20_control_te'], expression_vectors['20_TD'], alternative='two-sided')

    rho_20_treated_te_abundance, pvalue_20_treated_te_abundance = scipy.stats.mannwhitneyu(
        expression_vectors['20_treated_te'], expression_vectors['20_TT'], alternative='two-sided')

    sys.stdout.write("Plotting\n")

    # Plots
    # Scatter
    # Fold Change vs Abundance
    fc_abundance_11_trace = go.Scatter(
        y=expression_vectors['log_11_TD'],
        x=expression_vectors['11_TE_FC'],
        text=transcript_labels,
        mode='markers'
    )

    fc_abundance_20_trace = go.Scatter(
        y=expression_vectors['log_20_TD'],
        x=expression_vectors['20_TE_FC'],
        text=transcript_labels,
        mode='markers'
    )

    # Correlation Plots
    # pd_trace = go.Scatter(
    #     x=expression_vectors['11_PD'],
    #     y=expression_vectors['20_PD'],
    #     text=transcript_labels,
    #     mode='markers'
    # )
    #
    # pt_trace = go.Scatter(
    #     x=expression_vectors['11_PT'],
    #     y=expression_vectors['20_PT'],
    #     text=transcript_labels,
    #     mode='markers'
    # )
    #
    # td_trace = go.Scatter(
    #     x=expression_vectors['11_TD'],
    #     y=expression_vectors['20_TD'],
    #     text=transcript_labels,
    #     mode='markers'
    # )
    #
    # tt_trace = go.Scatter(
    #     x=expression_vectors['11_TT'],
    #     y=expression_vectors['20_TT'],
    #     text=transcript_labels,
    #     mode='markers'
    # )

    te_abundance_control_11_trace = go.Scatter(
        x=expression_vectors['11_control_te'],
        y=expression_vectors['log_11_TD'],
        text=transcript_labels,
        mode='markers'
    )

    te_abundance_treated_11_trace = go.Scatter(
        x=expression_vectors['11_treated_te'],
        y=expression_vectors['log_11_TT'],
        text=transcript_labels,
        mode='markers'
    )

    te_abundance_control_20_trace = go.Scatter(
        x=expression_vectors['20_control_te'],
        y=expression_vectors['log_20_TD'],
        text=transcript_labels,
        mode='markers'
    )

    te_abundance_treated_20_trace = go.Scatter(
        x=expression_vectors['20_treated_te'],
        y=expression_vectors['log_20_TT'],
        text=transcript_labels,
        mode='markers'
    )

    # control_te_trace = go.Scatter(
    #     x=expression_vectors['11_control_te'],
    #     y=expression_vectors['20_control_te'],
    #     text=transcript_labels,
    #     mode='markers'
    # )
    #
    # treated_te_trace = go.Scatter(
    #     x=expression_vectors['11_treated_te'],
    #     y=expression_vectors['20_treated_te'],
    #     text=transcript_labels,
    #     mode='markers'
    # )

    # Histograms
    human_11_TE_hist_trace = go.Histogram(
        x=expression_vectors['human_11_TE_FC'],
        histnorm='probability'
    )

    human_20_TE_hist_trace = go.Histogram(
        x=expression_vectors['human_20_TE_FC'],
        histnorm='probability'
    )

    human_Avg_TE_hist_trace = go.Histogram(
        x=expression_vectors['human_Avg_TE_FC'],
        histnorm='probability'
    )

    kshv_11_TE_hist_trace = go.Histogram(
        x=expression_vectors['kshv_11_TE_FC'],
        histnorm='probability'
    )

    kshv_20_TE_hist_trace = go.Histogram(
        x=expression_vectors['kshv_20_TE_FC'],
        histnorm='probability'
    )

    kshv_Avg_TE_hist_trace = go.Histogram(
        x=expression_vectors['kshv_Avg_TE_FC'],
        histnorm='probability'
    )

    comb_11_TE_hist_trace = go.Histogram(
        x=expression_vectors['11_TE_FC'],
        histnorm='probability'
    )

    comb_20_TE_hist_trace = go.Histogram(
        x=expression_vectors['20_TE_FC'],
        histnorm='probability'
    )

    comb_Avg_TE_hist_trace = go.Histogram(
        x=expression_vectors['Avg_TE_FC'],
        histnorm='probability'
    )

    # pd_data = [pd_trace]
    # pt_data = [pt_trace]
    # td_data = [td_trace]
    # tt_data = [tt_trace]
    #
    # control_te_data = [control_te_trace]
    # treated_te_data = [treated_te_trace]

    fc_abundance_11_data = [fc_abundance_11_trace]
    fc_abundance_20_data = [fc_abundance_20_trace]

    human_11_TE_hist_data = [human_11_TE_hist_trace]
    human_20_TE_hist_data = [human_20_TE_hist_trace]
    human_Avg_TE_hist_data = [human_Avg_TE_hist_trace]

    kshv_11_TE_hist_data = [kshv_11_TE_hist_trace]
    kshv_20_TE_hist_data = [kshv_20_TE_hist_trace]
    kshv_Avg_TE_hist_data = [kshv_Avg_TE_hist_trace]

    comb_11_TE_hist_data = [comb_11_TE_hist_trace]
    comb_20_TE_hist_data = [comb_20_TE_hist_trace]
    comb_Avg_TE_hist_data = [comb_Avg_TE_hist_trace]

    te_abundance_control_11_data = [te_abundance_control_11_trace]
    te_abundance_treated_11_data = [te_abundance_treated_11_trace]
    te_abundance_control_20_data = [te_abundance_control_20_trace]
    te_abundance_treated_20_data = [te_abundance_treated_20_trace]

    te_abundance_control_11_layout = dict(title="11 TE VS Abundance (TD) Control Two-Sided Mann Whitney U: {}, P-value:"
                                                " {}".format(rho_11_control_te_abundance,
                                                             pvalue_11_control_te_abundance),
                                          yaxis=dict(zeroline=False),
                                          xaxis=dict(zeroline=False)
                                          )

    te_abundance_treated_11_layout = dict(title="11 TE VS Abundance (TT) Treated Two-Sided Mann Whitney U: {}, P-value:"
                                                " {}".format(rho_11_treated_te_abundance,
                                                             pvalue_11_treated_te_abundance),
                                          yaxis=dict(zeroline=False),
                                          xaxis=dict(zeroline=False)
                                          )

    te_abundance_control_20_layout = dict(title="20 TE VS Abundance (TD) Control Two-Sided Mann Whitney U: {}, P-value:"
                                                " {}".format(rho_20_control_te_abundance,
                                                             pvalue_20_control_te_abundance),
                                          yaxis=dict(zeroline=False),
                                          xaxis=dict(zeroline=False)
                                          )

    te_abundance_treated_20_layout = dict(title="20 TE VS Abundance (TT) Treated Two-Sided Mann Whitney U: {}, P-value:"
                                                " {}".format(rho_20_treated_te_abundance,
                                                             pvalue_20_treated_te_abundance),
                                          yaxis=dict(zeroline=False),
                                          xaxis=dict(zeroline=False)
                                          )

    fc_abundance_11_layout = dict(title="Sample 11: Translational Efficiency Fold-Change VS Abundance (TD) "
                                        "Two-Sided Mann Whitney U: {}, P-value: {}".format(rho_11, pvalue_11),
                                  yaxis=dict(zeroline=False),
                                  xaxis=dict(zeroline=False)
                                  )

    fc_abundance_20_layout = dict(title="Sample 20: Translational Efficiency Fold-Change VS Abundance (TD) "
                                        "Two-Sided Mann Whitney U: {}, P-value: {}".format(rho_20, pvalue_20),
                                  yaxis=dict(zeroline=False),
                                  xaxis=dict(zeroline=False)
                                  )

    # Histograms
    human_11_TE_hist_layout = dict(title="Frequency of TE Fold Change Values in Sample 11 from {} human transcripts"
                                         "".format(len(expression_vectors['human_11_TE_FC'])),
                                   yaxis=dict(zeroline=False),
                                   xaxis=dict(zeroline=False)
                                   )

    human_20_TE_hist_layout = dict(title="Frequency of TE Fold Change Values in Sample 20 from {} human transcripts"
                                         "".format(len(expression_vectors['human_20_TE_FC'])),
                                   yaxis=dict(zeroline=False),
                                   xaxis=dict(zeroline=False)
                                   )

    human_Avg_TE_hist_layout = dict(title="Frequency of TE Fold Change Values in Sample Average from {} human "
                                          "transcripts".format(len(expression_vectors['human_Avg_TE_FC'])),
                                    yaxis=dict(zeroline=False),
                                    xaxis=dict(zeroline=False)
                                    )

    kshv_11_TE_hist_layout = dict(title="Frequency of TE Fold Change Values in Sample 11 from {} KSHV transcripts"
                                        "".format(len(expression_vectors['kshv_11_TE_FC'])),
                                  yaxis=dict(zeroline=False),
                                  xaxis=dict(zeroline=False)
                                  )

    kshv_20_TE_hist_layout = dict(title="Frequency of TE Fold Change Values in Sample 20 from {} KSHV transcripts"
                                        "".format(len(expression_vectors['kshv_20_TE_FC'])),
                                  yaxis=dict(zeroline=False),
                                  xaxis=dict(zeroline=False)
                                  )

    kshv_Avg_TE_hist_layout = dict(title="Frequency of TE Fold Change Values in Sample Average from {} KSHV transcripts"
                                         "".format(len(expression_vectors['kshv_Avg_TE_FC'])),
                                   yaxis=dict(zeroline=False),
                                   xaxis=dict(zeroline=False)
                                   )

    comb_11_TE_hist_layout = dict(title="Frequency of TE Fold Change Values in Sample 11 from {} KSHV and Human "
                                        "transcripts".format(len(expression_vectors['11_TE_FC'])),
                                  yaxis=dict(zeroline=False),
                                  xaxis=dict(zeroline=False)
                                  )

    comb_20_TE_hist_layout = dict(title="Frequency of TE Fold Change Values in Sample 20 from {} KSHV and Human "
                                        "transcripts".format(len(expression_vectors['20_TE_FC'])),
                                  yaxis=dict(zeroline=False),
                                  xaxis=dict(zeroline=False)
                                  )

    comb_Avg_TE_hist_layout = dict(title="Frequency of TE Fold Change Values in Sample Average from {} KSHV and Human "
                                         "transcripts".format(len(expression_vectors['Avg_TE_FC'])),
                                   yaxis=dict(zeroline=False),
                                   xaxis=dict(zeroline=False)
                                   )

    # Correlation Plots
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

    te_abundance_control_11_fig = dict(data=te_abundance_control_11_data, layout=te_abundance_control_11_layout)
    plotly.offline.plot(te_abundance_control_11_fig,
                        filename="{}_11_TE_Abundance_Control_plots".format(args.outroot))

    te_abundance_treated_11_fig = dict(data=te_abundance_treated_11_data, layout=te_abundance_treated_11_layout)
    plotly.offline.plot(te_abundance_treated_11_fig,
                        filename="{}_11_TE_Abundance_Treated_plots".format(args.outroot))

    te_abundance_control_20_fig = dict(data=te_abundance_control_20_data, layout=te_abundance_control_20_layout)
    plotly.offline.plot(te_abundance_control_20_fig,
                        filename="{}_20_TE_Abundance_Control_plots".format(args.outroot))

    te_abundance_treated_20_fig = dict(data=te_abundance_treated_20_data, layout=te_abundance_treated_20_layout)
    plotly.offline.plot(te_abundance_treated_20_fig,
                        filename="{}_20_TE_Abundance_Treated_plots".format(args.outroot))

    # control_te_fig = dict(data=control_te_data, layout=control_te_layout)
    # plotly.offline.plot(control_te_fig, filename="{}_Control_TE_correlation_plots".format(args.outroot))
    #
    # treated_te_fig = dict(data=treated_te_data, layout=treated_te_layout)
    # plotly.offline.plot(treated_te_fig, filename="{}_Treated_TE_correlation_plots".format(args.outroot))

    fc_fig_11 = dict(data=fc_abundance_11_data, layout=fc_abundance_11_layout)
    plotly.offline.plot(fc_fig_11, filename="{}_11_TE_FC_Abundance_Correlation".format(args.outroot))

    fc_fig_20 = dict(data=fc_abundance_20_data, layout=fc_abundance_20_layout)
    plotly.offline.plot(fc_fig_20, filename="{}_20_TE_FC_Abundance_Correlation".format(args.outroot))

    human_te_hist_fig_11 = dict(data=human_11_TE_hist_data, layout=human_11_TE_hist_layout)
    plotly.offline.plot(human_te_hist_fig_11, filename="{}_human_11_TE_hist".format(args.outroot))

    human_te_hist_fig_20 = dict(data=human_20_TE_hist_data, layout=human_20_TE_hist_layout)
    plotly.offline.plot(human_te_hist_fig_20, filename="{}_human_20_TE_hist".format(args.outroot))

    human_te_hist_fig_Avg = dict(data=human_Avg_TE_hist_data, layout=human_Avg_TE_hist_layout)
    plotly.offline.plot(human_te_hist_fig_Avg, filename="{}_human_Avg_TE_hist".format(args.outroot))

    kshv_te_hist_fig_11 = dict(data=kshv_11_TE_hist_data, layout=kshv_11_TE_hist_layout)
    plotly.offline.plot(kshv_te_hist_fig_11, filename="{}_kshv_11_TE_hist".format(args.outroot))

    kshv_te_hist_fig_20 = dict(data=kshv_20_TE_hist_data, layout=kshv_20_TE_hist_layout)
    plotly.offline.plot(kshv_te_hist_fig_20, filename="{}_kshv_20_TE_hist".format(args.outroot))

    kshv_te_hist_fig_Avg = dict(data=kshv_Avg_TE_hist_data, layout=kshv_Avg_TE_hist_layout)
    plotly.offline.plot(kshv_te_hist_fig_Avg, filename="{}_kshv_Avg_TE_hist".format(args.outroot))

    comb_te_hist_fig_11 = dict(data=comb_11_TE_hist_data, layout=comb_11_TE_hist_layout)
    plotly.offline.plot(comb_te_hist_fig_11, filename="{}_comb_11_TE_hist".format(args.outroot))

    comb_te_hist_fig_20 = dict(data=comb_20_TE_hist_data, layout=comb_20_TE_hist_layout)
    plotly.offline.plot(comb_te_hist_fig_20, filename="{}_comb_20_TE_hist".format(args.outroot))

    comb_te_hist_fig_Avg = dict(data=comb_Avg_TE_hist_data, layout=comb_Avg_TE_hist_layout)
    plotly.offline.plot(comb_te_hist_fig_Avg, filename="{}_comb_Avg_TE_hist".format(args.outroot))
