from collections import deque
import multiprocessing
import warnings
import numpy as np
import pandas as pd

import statsmodels.api as sm

from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr

from time import time

def gamma_run_loop(df):
    """ 
    
    Required columns: ["blockSizes_case*", "blockSizes_ctr*"]

    """
    de = deque()
    # case = df.filter(regex="blockSizes_case")/100
    # ctr = df.filter(regex="blockSizes_ctr")/100
    case = df.filter(regex="blockSizes_case")
    ctr = df.filter(regex="blockSizes_ctr")
    constant = [1] * case.shape[1] + [0] * ctr.shape[1]
    df = pd.concat([df[["chrom", "chromStart"]], case,  ctr], axis=1).to_numpy()
    print(df)
    print(constant)
    print(case)
    print(ctr)
    for row in df:
        y = pd.DataFrame([row[2:], constant])
        # print(y)
        y = y.dropna(axis=1)
        X = y.iloc[1].astype(int)
        y = y.iloc[0]
        # print(y)
        # if values are constant
        if len(y) == 0 or sum(y)/len(y) == y.iloc[0]:
            continue
        X = sm.add_constant(X)
        # print(X, y)
        model = sm.GLM(y, X, family=sm.families.Gamma())
        result = model.fit(disp=0)
        # print(result, result.pvalues, result.pvalues[1])
        de.append((row[0], row[1], result.pvalues[1]))
    return de
# helper function
def position_gamma(data, processes):
    """
    
    Required columns ["chrom", "chromStart", "p-val", "blockSizes_case*", "blockSizes_ctr*"]

    """
    dataframes = [data[data["chrom"] == "chr" + x] for x in [str(y) for y in range(1,23)] + ["X", "Y", "M"]]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with multiprocessing.Pool(processes=processes) as pool:
            results = pool.map(gamma_run_loop, dataframes)

    # merge and put into dataframe format
    results = pd.concat([pd.DataFrame(x).rename({0:"chrom",1:"chromStart",2:"p-val"}, axis=1) for x in results], axis=0, ignore_index=True)
    results = results.merge(data, on=["chrom","chromStart"], how="inner")
    print(results)
    results = results[[x for x in results.columns if x != "p-val"] + ["p-val"]]
    print(results)
    return results.sort_values(["chrom", "chromStart"])
# helper function
def position_limma(data, features, test_factor, model):
    """
    
    Required columns: ["chrom", "chromStart", "blockSizes_case*", "blockSizes_ctr*"]

    Columns in the features csv must match the name given in the blockSizes columns (for example, use the name "test" for "blockSizes_ctr_test"). 

    """


    # Enable conversion between pandas DataFrame and R data frame
    pandas2ri.activate()

    limma = importr('limma')
    # prepare dataset
    final = data.set_index(data["chrom"] + "_" + data["chromStart"].astype(str))
    m_df = final.filter(regex=("blockSizes_.*"))

    # prepare features for linear regression
    if features is None:
        design_matrix = pd.DataFrame()
        design_matrix["Sample"] = ["_".join(x.split("_")[2:]) for x in m_df.columns]
        design_matrix = design_matrix.set_index(["Sample"])
    else:
        design_matrix = pd.read_csv(features, index_col=0)
        print("design_matrix", design_matrix)
        print("loc columns", ["_".join(x.split("_")[2:]) for x in m_df.columns])
        design_matrix = design_matrix.loc[["_".join(x.split("_")[2:]) for x in m_df.columns]]
    design_matrix["Group"] = [int("case" in x) for x in m_df.columns]
    design_matrix["Intercept"] = [1] * len(m_df.columns)
    print({x:(x[x.find("_")+1:]) for x in m_df.columns})
    m_df = m_df.rename({x:(x[x.find("_")+1:]) for x in m_df.columns}, axis='columns')
    
    

    df_r = pandas2ri.py2rpy(m_df)
    r_design = pandas2ri.py2rpy(design_matrix)

    if not limma.is_fullrank(r_design)[0]:
        print(f"Column(s) {limma.nonEstimable(r_design)} are linear combinations of each other. This is not allowed. No p-values have been calculated.")
        return

    start_time = time()
    print(f"Checked if design matrix is full rank: {time() - start_time} seconds")

    fit = limma.lmFit(df_r, r_design)
    start_time = time()

    contrast_matrix = limma.makeContrasts(Group=test_factor, levels=r_design)
    start_time = time()

    fit = limma.contrasts_fit(fit, contrast_matrix)
    start_time = time()

    py_fit = pandas2ri.rpy2py(fit)
    stdev_unscaled = np.array(py_fit.rx2("stdev.unscaled"))
    stdev_unscaled[stdev_unscaled < 0.1] = 0.1
    py_fit.rx2["stdev.unscaled"] = stdev_unscaled
    fit = pandas2ri.py2rpy(py_fit)
    start_time = time()

    if model.lower() == "ebayes":
        # results = limma.eBayes(fit, trend=False, robust=False, proportion=1, stdev_coef_lim=r['c'](0.1, np.inf))
        results = limma.eBayes(fit, trend=False, robust=False, proportion=1, stdev_coef_lim=r['c'](0.1, np.inf), winsor_tail_p = 0.03)
    elif model.lower() == "treat":
        results = limma.treat(fit, lfc=0.2, trend=False, robust=False)
    start_time = time()

    results2 = limma.topTable(results, adjust_method='BH', number=data.shape[0], sort_by="p")
    
    final["p-val"] = results.rx2('p.value')

    results_df = pandas2ri.rpy2py(results2)

    final["q-value"] = results_df['adj.P.Val'] #.rx2('adj.P.Val')
    start_time = time()

    final = final.reset_index(drop=True)
    start_time = time()

    final_sorted = final.sort_values(["chrom", "chromStart"])
    start_time = time()

    print(f"Initial data shape: {data.shape}")
    print(f"Final data shape after set_index: {final.shape}")
    print(f"m_df shape after filtering: {m_df.shape}")
    print(f"Design matrix shape: {design_matrix.shape}")
    # print(f"R dataframe shape (m_df): {df_r.shape}")
    # print(f"R dataframe shape (design_matrix): {r_design.shape}")
    print(f"Final shape before p-value assignment: {final.shape}")
    # print(f"Results p-value shape: {results.rx2('p.value').shape}")
    print(f"Final sorted shape: {final_sorted.shape}")

    return final_sorted
