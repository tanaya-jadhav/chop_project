import pandas as pd

def main():
    combined_out_file = "./outcombined.tsv"
    second_out_file = "./out_chr22.tsv"

    f = pd.read_csv(combined_out_file, sep='\t', header=0, index_col=0)
    s = pd.read_csv(second_out_file, sep='\t', header=0, index_col=0)
    print(f.shape)
    print(s.shape)
    f_genes = f.index.values.tolist()
    header_list = list(f)
    s_genes = s.index.values.tolist()

    for gene in s_genes:
        if gene in f_genes:
            print(gene)
            for sample in header_list:
                f.loc[gene][sample] = f.loc[gene][sample] + s.loc[gene][sample]

    f = f.combine_first(s)
    print(f.shape)
    f.to_csv("./outcombined.tsv", sep='\t')





if __name__ == '__main__':
    main()