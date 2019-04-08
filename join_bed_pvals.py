import pandas as pd
from itertools import islice
import seaborn as sns
import matplotlib.pyplot as plt
from decimal import Decimal

def main():
    # bedfilepath = './MultiPlatform.2.merged.bed'
    # pvalpath = './coverage_pvals_new.tsv'
    #
    # bed = pd.read_csv(bedfilepath, sep='\t', header=None)
    # pval = pd.read_csv(pvalpath, sep='\t', header=0)
    # print(bed.head(5))
    # print(pval.head(5))
    # joined = pd.concat([bed, pval], axis=1)
    # joined.to_csv('./combined_pvals_new.bed', sep= '\t', header=False, index=False)
    pvals = pd.read_csv('./coverage_pvals_differences.tsv', sep='\t', header=None)
    print(pvals.head(5))

    # file_path = './intersect_result.bed'
    # intersect = pd.read_csv(file_path, sep='\t', header=None)
    # print(intersect.head(5))
    # df = intersect[[0, 1, 2, 3, 4, 5, 6, 10]]
    # df.to_csv('./intersected_withgenes_means.tsv', sep='\t', header=False, index=False)

    # df = pd.read_csv('./intersected_withgenes_means.tsv', sep='\t', header=None)
    # df.columns = ['chr', 'start', 'stop', 'IBD_average', 'control_average', 'pval', 'adjusted_pval', 'gene']
    # # print(df.head(5))
    # df = df.groupby(['chr', 'start', 'stop'])
    # df = df['gene'].apply(', '.join).reset_index()
    # print(df.head(5))
    # joined = pd.concat([df, pvals], axis=1)
    # print(joined.head(5))
    # joined = joined[['chr', 'start', 'stop', 0, 1, 2, 3, 'gene']]
    # print(joined.head(5))
    # joined.to_csv('./withgenes_test.tsv', sep= '\t', index=False)

    ##for filtering ranges where difference between averages is >0.1

    dif = pd.read_csv('./withgenes_test.tsv', sep= '\t', header=0)
    print(dif.head(5))


    bed_list =[]
    for index, row in dif.iterrows():
        if row[5] < 0.1 and row[6] > 0.01:
            # print(row[5])
            bed_dict = {}
            bed_dict['chr'] = row['chr']
            bed_dict['start'] = row['start']
            bed_dict['stop'] = row['stop']
            bed_list.append(pd.Series(bed_dict))
    # print(bed_list)

    bed = pd.concat(bed_list, axis=1).T
    # print(bed.shape)
    bed.to_csv('./filteredranges1.bed', sep='\t', header=False, index=False)
    # df = pd.read_csv('./withgenes.tsv', sep= '\t')
    # tot_obs = len(df['adjusted_pval'])
    # list_hundthousand = []
    # fdrlist_hundthousand = []
    # list_tenthousand = []
    # fdrlist_tenthousand = []
    # list_thousand = []
    # fdrlist_thousand = []
    # list_hundred = []
    # fdrlist_hundred = []
    # list_low = []
    # fdrlist_low = []
    # for index, row in df.iterrows():
    #     val = row['pval']
    #     fdr_val = row['adjusted_pval']
    #     if  val < 0.00001:
    #         list_hundthousand.append(val)
    #         fdrlist_hundthousand.append(fdr_val)
    #     elif val < 0.0001:
    #         list_tenthousand.append(val)
    #         fdrlist_tenthousand.append(fdr_val)
    #     elif val < 0.001:
    #         list_thousand.append(val)
    #         fdrlist_thousand.append(fdr_val)
    #     elif val < 0.01:
    #         list_hundred.append(val)
    #         fdrlist_hundred.append(fdr_val)
    #     elif val < 0.05:
    #         list_low.append(val)
    #         fdrlist_low.append(fdr_val)
    # list_tenthousand = list_tenthousand + list_hundthousand
    # list_thousand = list_thousand + list_tenthousand
    # list_hundred = list_hundred + list_thousand
    # list_low = list_low + list_hundred
    # barplot_data = [['<0.05', len(list_low)], ['<0.01', len(list_hundred)], ['<0.001', len(list_thousand)],
    #                 ['<0.0001', len(list_tenthousand)], ['<0.00001', len(list_hundthousand)]]
    # bp = pd.DataFrame(barplot_data, columns=['Categories', 'Num_pvals'])
    # print(bp)
    # sns.set(style='whitegrid')
    # ax = sns.barplot(x='Categories', y='Num_pvals', data=bp)
    # max_fdr = [max(fdrlist_low), max(fdrlist_hundred), max(fdrlist_thousand), max(fdrlist_tenthousand), max(fdrlist_hundthousand)]
    # print(max_fdr)
    # count = 0
    # for p in ax.patches:
    #     maximum = max_fdr[count]
    #     percent = (p.get_height()/tot_obs)*100
    #     disp = "{0:.2f}".format(percent) + '% ' + '\n' + 'max fdr ' + "%.2E" %Decimal(maximum)
    #     ax.text(p.get_x() + p.get_width() / 2., p.get_height(), disp,
    #             fontsize=12, color='red', ha='center', va='bottom')
    #     count = count + 1
    # plt.show()



if __name__ == '__main__':
    main()