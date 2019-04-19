import pandas as pd
import numpy as np
from glob import iglob, glob
from scipy import stats
# from rpy2.robjects.packages import importr
# from rpy2.robjects.vectors import FloatVector
from statsmodels.stats.multitest import multipletests


def getListOfCoverage(file):
    sample_name = file.split('.sorted')[0].split('/')[-1]
    df = pd.read_csv(file, sep='\t', header=None, index_col=False)
    coverage = df[6]
    cov_list = []
    for line in coverage:
        cov_list.append(line)
    return sample_name, cov_list


def main():
    coverage_dir = '/Users/jadhavt/Desktop/coveragefiles/'
    IBD_dict = {}
    control_dict = {}
    IBD_filetypes = (coverage_dir + "IR*.coverage", coverage_dir + "*IBD*.coverage")
    for pattern in IBD_filetypes:
        for file in iglob(pattern):
            # print(file)
            sample_name, cov_list = getListOfCoverage(file)
            IBD_dict[sample_name] = cov_list
    for file in iglob(coverage_dir + "11*"):
        sample_name, cov_list = getListOfCoverage(file)
        control_dict[sample_name] = cov_list
    # print(IBD_dict)
    # print(control_dict)

    p_list = []
    average_IBD = []
    average_control = []
    average_diff = []
    for row in range(200451):
        IBD_list = []
        control_list = []
        for sample, cov in IBD_dict.items():
            IBD_list.append(cov[row])
        for sample, cov in control_dict.items():
            control_list.append(cov[row])
        mean_ibd = np.mean(IBD_list)
        mean_cont = np.mean(control_list)
        average_IBD.append(mean_ibd)
        average_control.append(mean_cont)
        average_diff.append(abs(mean_ibd - mean_cont))
        t, p = stats.ttest_ind(IBD_list, control_list)
        p_list.append(p)
        # print(p)
    print("plist", len(p_list))
    print("average diff", len(average_diff), average_diff[:10])

    # rstats = importr('stats')

    # p_adjust = rstats.p_adjust(FloatVector(p_list), method = 'BH')

    nullindex_list = []
    count = 0
    for p in p_list:
        if pd.isnull(p):
            nullindex_list.append(count)
        count = count + 1
    print("null index", nullindex_list)

    cleaned_p_list = [x for x in p_list if ~np.isnan(x)]

    for p in cleaned_p_list:
        if pd.isnull(p):
            print(p)

    print("cleaned p list", len(cleaned_p_list))


    p_sm_adjusted = multipletests(cleaned_p_list, method='fdr_bh')
    p_adjusted = p_sm_adjusted[1]

    print("p adjusted", len(p_adjusted))


    p_adjusted_list = list(p_adjusted)
    for index in nullindex_list:
        p_adjusted_list.insert(index, np.nan)

    # print(len(p_adjusted), p_adjusted)
    # print(p_adjusted_list)

    print("new p adjusted", len(p_adjusted_list))

    data = list(zip(average_IBD, average_control, average_diff, p_list))
    p_df = pd.DataFrame(data, columns=['IBD_average', 'control_average', 'average_difference', 'ttest_p-vals'])
    p_df['adjusted_pvals'] = p_adjusted_list
    p_df.to_csv('./coverage_pvals_differences.tsv', sep='\t', index=False, header=False)











if __name__ == '__main__':
    main()