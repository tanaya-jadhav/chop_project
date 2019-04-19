import pandas as pd


def main():
    bedfile_path = '/Users/jadhavt/Desktop/coveragefiles/new_uscsgenes.bed'
    refbed = pd.read_csv(bedfile_path, sep='\t')

    seqbed_path = '/Users/jadhavt/Desktop/coveragefiles/MultiPlatform.2.merged.bed'
    seqbed = pd.read_csv(seqbed_path, sep='\t', header=None)
    # print(refbed.head(5))
    print(seqbed.head(5))

    ##code for removing "chr" from chromosome column
    # refbed['chr'] = ''
    # for index, row in refbed.iterrows():
    #     refbed.iloc[index, 4] = str(row['#chrom']).strip('chr')
    # refbed = refbed[['chr', 'txStart', 'txEnd', 'name2']]
    # refbed.to_csv('./new_uscsgenes.bed', sep='\t', index=False, header=True)

    genedict = {}
    for index, row in refbed.iterrows():
        gene_tup = (row['chr'], row['txStart'], row['txEnd'])
        if gene_tup not in genedict:
            genedict[gene_tup] = row['name2']
    # print(genedict)
    seqbed['gene'] = ''
    for index, row in seqbed.iterrows():
        for key in genedict:
            if key[0] == str(row[0]):
                if key[1] <= row[1] and key[2] > row[1]:
                    seqbed.iloc[index, 3] = genedict[key]
                elif key[1] < row[2] and key[2] >= row[2]:
                    seqbed.iloc[index, 3] = genedict[key]
    print(seqbed)


if __name__ == '__main__':
    main()