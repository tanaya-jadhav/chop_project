import pandas as pd
import scipy.stats as ss
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
import math


def find_lines_toskip(vcf_path):
    with open(vcf_path) as fv:
        counter = 0
        for line in fv:
            if line.startswith('#'):
                counter = counter + 1
    toskip = counter - 1
    return toskip


def get_groupdict(dataframe_path):
    groupsdict = {}
    with open(dataframe_path) as gf:
        for line in gf:
            linelist = line.split('\t')
            groupID = linelist[0]
            markerIDs = linelist[1:]
            # print(markerIDs)
            for marker in markerIDs:
                groupsdict[marker] = groupID
    return groupsdict


def get_top_genes_list(df):
    top_genes_list = []
    for i, r in df.iterrows():
        mark_id = r["MARKER_ID"]
        gene = str(mark_id).split('_')[-1]
        top_genes_list.append(gene)
    return top_genes_list


def main():
    top1000path = './test.gene.skat.epacts.top1000'
    top1000 = pd.read_csv(top1000path, sep='\t')
    # print(top1000)
    top_genes_list = get_top_genes_list(top1000)

    filteredvcf = 'ibd.control.chr1.anno.gnomadanno1.filtered.vcf'
    skip_lines = find_lines_toskip(filteredvcf)

    fv = pd.read_csv(filteredvcf, sep='\t', skiprows=skip_lines)

    gf_path = './ibd.control.chr1.anno.groupfile'
    with open(gf_path) as gf:
        groups = []
        for line in gf:
            group = line.split('\t')[0]
            groups.append(group)

    header_list = list(fv.columns.values)
    sample_list = header_list[9:]
    collapsed_df = pd.DataFrame(columns=sample_list, index=groups)

    groupsdict = get_groupdict(gf_path)

    for index, row in fv.iterrows():
        for sample in sample_list:
            geno = str(row[sample]).split(':')[0]
            if geno.split('/')[0] != '.' and geno.split('/')[1] != '.':
                alleles = int(geno.split('/')[0]) + int(geno.split('/')[1])
                chr = row['#CHROM']
                position = row['POS']
                ref = row['REF']
                alt = row['ALT']
                marker_id = str(chr) + ':' + str(position) + '_' + ref + '/' + alt
                if marker_id in groupsdict:
                    group_id = groupsdict[marker_id]
                    allele_count = collapsed_df.loc[group_id][sample]
                    print(allele_count)
                    if math.isnan(allele_count):
                        collapsed_df.loc[group_id][sample] = alleles
                        print(collapsed_df.loc[group_id][sample])
                    else:
                        collapsed_df.loc[group_id][sample] = allele_count + alleles
                        print(print(collapsed_df.loc[group_id][sample]))

    collapsed_df = collapsed_df.dropna(how='all')
    for i, r in collapsed_df.iterrows():
        if i in top_genes_list:
            pass
        else:
            collapsed_df = collapsed_df.drop(i)
    collapsed_df.to_csv('./collapsed_data2.tsv', sep='\t')














if __name__ == '__main__':
    main()