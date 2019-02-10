import pandas as pd
import scipy.stats as ss
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
import math


def add_family(micro_copy):
    microbe_list = list(micro_copy.columns.values)[911:]
    family_list = []
    genus_list = []
    for microbe in microbe_list:
        family = str(microbe).split()[0]
        genus = str(microbe).split()[:2]
        genus = ' '.join(genus)
        family_list.append(family)
        genus_list.append(genus)
    family_list = list(set(family_list))
    genus_list = list(set(genus_list))
    for family in family_list:
        micro_copy[family] = 0
    for genus in genus_list:
        micro_copy[genus] = 0
    microbe_df = micro_copy.loc[:, 'Acidobacteria Granulicella unclassified':'Verrucomicrobia Akkermansia muciniphila']
    for family in family_list:
        fam_df = microbe_df.filter(like=family)
        micro_copy[family] = fam_df.sum(axis=1)
    for genus in genus_list:
        gen_df = microbe_df.filter(like=genus)
        micro_copy[genus] = gen_df.sum(axis=1)
    micro_copy.to_csv('./microcopy_withfamilyandgenus.tsv', sep='\t')

def wilcoxon_test(micro_veo, micro_control, micro_copy):
    metabolite_list = list(micro_copy.columns.values)[16:911]
    rst_list = []
    for metabolite in metabolite_list:
        rst_dict = {}
        rst_dict['metabolite'] = metabolite
        x = micro_veo[metabolite]
        y = micro_control[metabolite]
        wilcox = ss.ranksums(x, y)
        rst_dict['p-value'] = wilcox[1]
        rst_list.append(rst_dict)
    rst = pd.DataFrame(rst_list)
    rst.to_csv('./wilcoxon_text.tsv', sep='\t', index=False)


def split_correlation(micro_veo, micro_control, metabolite_list, microbe_list, outputfile):
    corr_list = []
    for metabolite in metabolite_list:
        for microbe in microbe_list:
            corr_dict = {}
            corr_dict['metabolite'] = metabolite
            x = micro_veo[metabolite]
            x2 = micro_control[metabolite]
            corr_dict['microbe'] = microbe
            y = micro_veo[microbe]
            y2 = micro_control[microbe]
            corr_dict['veo_correlation_coefficient'] = np.corrcoef(x.astype(float), y.astype(float))[1, 0]
            corr_dict['control_correlation_coefficient'] = np.corrcoef(x2.astype(float), y2.astype(float))[1, 0]
            corr_list.append(corr_dict)
    corr_analysis = pd.DataFrame(corr_list)
    corr_analysis = corr_analysis[
        ['microbe', 'metabolite', 'veo_correlation_coefficient', 'control_correlation_coefficient']]
    # corr_analysis = corr_analysis[pd.notnull(corr_analysis['veo_correlation_coefficient'])]
    corr_analysis = corr_analysis.dropna()
    corr_analysis.to_csv(outputfile, sep='\t', index=False)


def combined_correlation(dataframe, metabolite_list, microbe_list, outputfile):
    corr_list = []
    for metabolite in metabolite_list:
        for microbe in microbe_list:
            corr_dict = {}
            corr_dict['metabolite'] = metabolite
            x = dataframe[metabolite]
            corr_dict['microbe'] = microbe
            y = dataframe[microbe]
            corr_dict['correlation_coefficient'] = np.corrcoef(x.astype(float), y.astype(float))[1, 0]
            corr_list.append(corr_dict)
    corr_analysis = pd.DataFrame(corr_list)
    corr_analysis = corr_analysis[['microbe', 'metabolite', 'correlation_coefficient']]
    corr_analysis = corr_analysis[pd.notnull(corr_analysis['correlation_coefficient'])]
    corr_analysis.to_csv(outputfile, sep='\t', index=False)


def get_sig_pairs(filepath):
    sig_pair = pd.read_csv(filepath, sep='\t', header=None)
    sig_pair_list = []
    for index, row in sig_pair.iterrows():
        sig_tuple = (row[0], row[1])
        sig_pair_list.append(sig_tuple)
    return sig_pair_list


def get_seaborn_readable_data(dataframe, sig_pair_list):
    meta_dict = {}

    for index, row in dataframe.iterrows():
        for pair in sig_pair_list:
            metabolite = pair[1]
            microbe = pair[0]
            meta_dict[(str(row['SUBJECT_ID']), row['Cohort'], metabolite, microbe)] = (row[metabolite], row[microbe])

    seaborn_list = []
    for key in meta_dict:
        seaborn_dict = {}
        seaborn_dict['Subject_ID'] = key[0]
        seaborn_dict['Cohort'] = key[1]
        seaborn_dict['Metabolite'] = key[2]
        seaborn_dict['Microbe'] = key[3]
        seaborn_dict['met_abundance'] = meta_dict[key][0]
        seaborn_dict['mic_abundance'] = meta_dict[key][1]
        seaborn_list.append(seaborn_dict)

    seaborn_data = pd.DataFrame(seaborn_list)
    seaborn_path = './famgenus_seaborn_data.tsv'
    seaborn_data.to_csv(seaborn_path, sep='\t', index=False)


def main():
    # micro_copy_path = '/Users/jadhavt/Desktop/Microbiome/micro_meta_data.tsv'
    # micro_copy = pd.read_csv(micro_copy_path, sep='\t')
    # add_family(micro_copy)

    micro_fam = pd.read_csv('./microcopy_withfamilyandgenus.tsv', sep='\t')

    # sig_meta_list = ['4-hydroxychlorothalonil', 'ergothioneine', 'hippurate', 'perfluorooctanoate (PFOA)*',
    #                  'gentisate', 'lactosyl-N-palmitoyl-sphingosine (d18:1/16:0)', 'perfluorooctanesulfonate (PFOS)',
    #                  'hyocholate', '2-hydroxylaurate', 'maltose', 'N-acetylarginine', 'alpha-CEHC sulfate',
    #                  'Fibrinopeptide A*', 'alpha-tocopherol', 'linolenoylcarnitine (C18:3)*',
    #                  'N-stearoyl-sphinganine (d18:0/18:0)*', 'glycolithocholate sulfate*']
    # microbe_list = list(micro_fam.columns.values)[1554:1568]
    # microbe_genus_list = list(micro_fam.columns.values)[1568:]
    #
    # micro_veo = micro_fam.loc[micro_fam['Cohort'] == 1]
    # micro_control = micro_fam.loc[micro_fam['Cohort'] == 0]
    # output_famfile = './fam_split_correlation.tsv'
    # output_genusfile = './genus_split_correlation.tsv'
    # split_correlation(micro_veo, micro_control, sig_meta_list, microbe_list, output_famfile)
    # split_correlation(micro_veo, micro_control, sig_meta_list, microbe_genus_list, output_genusfile)

    # output_famfile = './fam_combined_correlation.tsv'
    # output_genusfile = './genus_combined_correlation.tsv'
    # combined_correlation(micro_fam, sig_meta_list, microbe_list, output_famfile)
    # combined_correlation(micro_fam, sig_meta_list, microbe_genus_list, output_genusfile)


    sig_pair_list = get_sig_pairs('./sig_pairs.tsv')

    # sig_pair_list = get_sig_pairs('significant_pairs_fam.tsv')
    # for pair in sig_pair_list:
    #     microbe = pair[0]
    #     metabolite = pair[1]
    #     try:
    #         p = sns.lmplot(x=metabolite, y=microbe, hue="Cohort", data=micro_fam)
    #         plotname = str(metabolite).replace(' ', '_') + "_" + str(microbe).replace(' ', '_') + "plot.png"
    #         p.savefig(plotname)
    #         print(plotname)
    #     except:
    #         print(plotname, 'skipped')
    #         pass

    get_seaborn_readable_data(micro_fam, sig_pair_list)
    seaborn_data = pd.read_csv('./famgenus_seaborn_data.tsv', sep='\t')
    g = sns.lmplot(x="met_abundance", y="mic_abundance", hue="Cohort", data=seaborn_data, col="Metabolite",
                   row="Microbe")
    g.savefig('gridplot.png')

    # for pair in sig_pair_list:
    #     microbe = pair[0]
    #     metabolite = pair[1]
    #     try:
    # p = sns.lmplot(x='perfluorooctanoate (PFOA)*', y='Firmicutes', hue="Cohort", data=micro_fam)
    # p.savefig('PFOA_Firmicutes.png')
    #         plotname = str(metabolite).replace(' ', '_').replace('_(d18:1/16:0)', '') + "_" + str(microbe).replace(' ', '_') + "plot.png"
    #         p.savefig(plotname)
    #         print(plotname)
    #     except:
    #         print(plotname, 'skipped')
    #         pass

    # g = sns.lmplot(x="met_abundance", y="mic_abundance", hue="Cohort", data=seaborn_data, col="Metabolite", row="Microbe")
    # g.savefig('gridplot.png')

if __name__ == '__main__':
    main()