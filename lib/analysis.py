from collections import deque
import multiprocessing
from queue import Queue
from typing import Callable
import numpy as np
import pandas as pd
import math
from collections import defaultdict


from statsmodels.stats import multitest
import polars.selectors as cs
import polars as pl
import re

from scipy.stats import mannwhitneyu, ttest_ind, ranksums, ks_2samp, median_test, brunnermunzel

from lib import InputProcessor

from .position_p_vals import position_gamma, position_limma
from .gene_analysis import window_based_gene, position_based_gene

class Analysis():
    """Example Usage
    
    >>> obj = analysis()
    >>> result = function(parameter_1,parameter_2,...,parameter_n)
    >>> result
    (returned_1,returned_2,...,returned_m)
    """
    def assert_required_columns(self, df, required_columns):
        columns = df.columns
        for col in required_columns:
            col = "^" + col + "$|^chrom" + col + "$"
            if not any(re.match(col, c, re.IGNORECASE) for c in columns):
                raise AssertionError(f"Required column pattern '{col}' not found in filter_samples_ratiocolumns: {columns}")
    
    def assert_one_of_column_pairs(self, df, column_pairs):
        columns = df.columns
        for pair in column_pairs:
            if all(any(re.match(col, c) for c in columns) for col in pair):
                return
        raise AssertionError(f"None of the required column pairs {column_pairs} found in columns: {columns}")

    def merge_tables(self, min_cov_individual = 10, min_cov_group = 15, cov_percentile = 100.0, min_samp_ctr = 2, min_samp_case = 2, filter_samples_ratio=0.6, meth_group_threshold=0.2, case_data: InputProcessor.data_container = None, ctr_data: InputProcessor.data_container = None, small_mean = 0.2):
        """Required columns:

        ["chrom", "chromStart"] + ONE OF [("coverage_*", "blockSizes_*"), ("positive_{name}", "negative_{name}")]
        
        """
        # if input is not a list
        if not isinstance(case_data, list):
            case_data = [case_data]
        if not isinstance(ctr_data, list):
            ctr_data = [ctr_data]
                      
        # assert columns based on above definitions for each data list
        for data in [case_data, ctr_data]:
            for df in data:
                self.assert_required_columns(df, ["chrom", "chromStart"])
                self.assert_one_of_column_pairs(df, [("coverage_.*", "blockSizes_.*"), ("positive_.*", "negative_.*")])

        
        MERGE_LIST = ["chrom", "chromStart"]

        if "chromEnd" in case_data[0].columns:
            MERGE_LIST.append("chromEnd")
        if "strand" in case_data[0].columns:
            MERGE_LIST.append("strand")
        

        d = {}
        d["ctr"] = [pl.from_pandas(x) for x in ctr_data] # Why from pandas
        d["case"] = [pl.from_pandas(x) for x in case_data] # Why from pandas

        # print(d)
        for key in ["ctr", "case"]:
           for i in range(len(d[key])):
                df = d[key][i]
                # print(df)
                for col in df.columns:
                   if "blockSizes" in col or "coverage" in col or "positive" in col or "negative" in col:
                       df = df.with_columns(pl.col(col).cast(pl.Float64, strict=False).alias(col)) # why convert to float64?

                # for positive, negative setup, calcualte coverage and methylation percentage
                positive_col = next((col for col in df.columns if col.startswith("positive_")), None) # why starting with positive
                negative_col = next((col for col in df.columns if col.startswith("negative_")), None) # why starting with negative
                # print(df)
                # print(df.columns)
                # name = "file_"+str(i+1)
                if positive_col and negative_col:
                    name = re.match(r"positive_(.*)", positive_col).group(1)
                    name2 = re.match(r"negative_(.*)", negative_col).group(1)
                    if f"negative_{name}" == negative_col:
                        df = df.with_columns((pl.col(positive_col) + pl.col(negative_col)).alias(f"coverage_{key}_{name}"))
                        df = df.with_columns((pl.col(positive_col) / (pl.col(positive_col) + pl.col(negative_col)) ).alias(f"blockSizes_{key}_{name}")) # why * 100?
                        df = df.drop([positive_col, negative_col])
                        cv_name = f"coverage_{key}_{name}"
                        bs_name = f"blockSizes_{key}_{name}"
                else:
                    bs_name = ''
                    cv_name = ''
                    for _c_col in df.columns:
                        if _c_col[:len(f"blockSizes_{key}")]==f"blockSizes_{key}": bs_name = _c_col
                        if _c_col[:len(f"coverage_{key}")]==f"coverage_{key}": cv_name = _c_col
                    if bs_name=='': print("Error")
                    else:
                        if ((df[bs_name] > 1.5).sum()) > max( (df.shape[0])*0.001, 1): #
                            df = df.with_columns((pl.col(bs_name) / 100).alias(bs_name))
                # if positive_col and negative_col:
                #     print(df)
                #     df = df.with_columns((pl.col(f"positive_{name}") + pl.col(f"negative_{name2}")).alias(f"coverage_{key}_{name}"))
                #     df = df.filter(pl.col(f"coverage_{key}_{name}") >= min_cov_individual) # add cov filter here!!!!!!
                # else:
                df = df.filter(pl.col(cv_name) >= min_cov_individual)
                # select only necessary columns for merge
                d[key][i] = df.select( MERGE_LIST + [col for col in df.columns if "blockSizes" in col or "coverage" in col]) #
        ctr = d["ctr"][0]
        case = d["case"][0]
        case_col = MERGE_LIST.copy()
        ctr_col = MERGE_LIST.copy()
        ctr_col.extend( d["ctr"][0].select(cs.starts_with("coverage_")).columns)
        ctr_col.extend( d["ctr"][0].select(cs.starts_with("blockSizes_")).columns)
        case_col.extend( d["case"][0].select(cs.starts_with("coverage_")).columns)
        case_col.extend( d["case"][0].select(cs.starts_with("blockSizes_")).columns)
        nbr_case = len(d["case"])
        nbr_ctr = len(d["ctr"])
        min_samp_case = max(min_samp_case, math.ceil(nbr_case * filter_samples_ratio))
        min_samp_ctr = max(min_samp_ctr, math.ceil(nbr_ctr * filter_samples_ratio))
        for e in range(1,len(d["ctr"])):
            ctr = ctr.join(d["ctr"][e], on=MERGE_LIST, how='full', coalesce=True)
            ctr_col.extend( d["ctr"][e].select(cs.starts_with("coverage_")).columns)
            ctr_col.extend( d["ctr"][e].select(cs.starts_with("blockSizes_")).columns)
        del(d["ctr"])
        for e in range(1,len(d["case"])):
            case = case.join(d["case"][e], on=MERGE_LIST, how='full', coalesce=True)
            case_col.extend( d["case"][e].select(cs.starts_with("coverage_")).columns)
            case_col.extend( d["case"][e].select(cs.starts_with("blockSizes_")).columns)
        del(d)
        df_all = ctr.join(case, on=MERGE_LIST, how="full")
        case = df_all.select(case_col)
        ctr = df_all.select(ctr_col)
        cov_columns = cs.starts_with("coverage_")
        filtered_meth_df = df_all.with_columns(sum_=pl.sum_horizontal((cov_columns >= min_cov_group))).select(pl.col("sum_") >= (min_samp_case + min_samp_ctr))["sum_"]
        case = case.filter(filtered_meth_df)
        ctr = ctr.filter(filtered_meth_df)
        case1 = case.with_columns(sum_case=pl.sum_horizontal(case.select((cov_columns >= min_cov_group))))["sum_case"]
        ctr1 = ctr.with_columns(sum_ctr=pl.sum_horizontal(ctr.select((cov_columns >= min_cov_group))))["sum_ctr"]
        ctr_cov_cond = ctr.with_columns(sum_ctr=pl.sum_horizontal((cov_columns >= min_cov_group))).select(pl.col("sum_ctr") >= min_samp_ctr)["sum_ctr"]
        case_cov_cond = case.with_columns(sum_case=pl.sum_horizontal((cov_columns >= min_cov_group))).select(pl.col("sum_case") >= min_samp_case)["sum_case"]
        case_mean = case.with_columns(mean_blockSizes=case.select(cs.contains("blockSizes_")).mean_horizontal())['mean_blockSizes']
        ctr_mean = ctr.with_columns(mean_blockSizes=ctr.select(cs.contains("blockSizes_")).mean_horizontal())['mean_blockSizes']
        group_filter = case_cov_cond & ctr_cov_cond
        if min_samp_case + min_samp_ctr < nbr_case + nbr_ctr:
            group_filter = group_filter | ( ( ( case1 >= min_samp_case ) & ( ctr1 < min_samp_ctr) & ( ctr_mean < small_mean ) ) |  ( ( ctr1 >= min_samp_ctr    ) & ( case1 < min_samp_case ) & ( case_mean < small_mean ) ) )
        case = case.filter(group_filter)
        ctr = ctr.filter(group_filter)
        case = case.to_pandas()
        ctr = ctr.to_pandas()
        avg_case = case.filter(like="blockSizes").mean(axis=1, numeric_only=True)
        avg_ctr = ctr.filter(like="blockSizes").mean(axis=1, numeric_only=True)
        case["avg_case"] = avg_case
        ctr["avg_ctr"] = avg_ctr
        f = pd.merge(case,ctr,on=MERGE_LIST, how='inner')
        f["diff"] = f["avg_case"] - f["avg_ctr"]
        self._merged_result = f.sort_values(by=MERGE_LIST).reset_index(drop=True)
        print("final output of filter", self._merged_result)
        return self._merged_result
    def position_based(self, data: InputProcessor.data_container, method="limma", features=None, test_factor="Group", processes=22, model="eBayes", min_std=0.1, fill_na:bool=True):
        """
    
        Required columns: ["chrom", "chromStart", "blockSizes_case.*", "blockSizes_ctr.*"]

        Columns in the features csv must match the name given in the blockSizes columns (for example, use the name "test" for "blockSizes_ctr_test"). 

        """
        assert method in ["gamma", "limma"], "Select a valid method parameter. Options are: [""gamma"", ""limma""]"
        assert isinstance(data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(data, ["chrom", "chromStart", "blockSizes_case.*", "blockSizes_ctr.*"])

        # filter standard deviation
        print(data.shape, min_std)
        # data = data[(data.filter(regex="blockSizes_case.*").std(axis=1) > min_std) & (data.filter(regex="blockSizes_ctr.*").std(axis=1) > min_std)]
        data2 = data.copy()
        data2["blockSizes_std"] = data2.filter(regex="blockSizes_*").std(axis=1)
        # data.to_csv("saved_res.csv")
        # data2.to_csv("std_saved_res.csv")
        del(data2)
        data = data[(data.filter(regex="blockSizes_*").std(axis=1) >= min_std)]
        if fill_na:
            # Define a minimum value for when all values in a row are NaN
            min_value = 10**(-3)  
            # Get columns for each group
            case_columns = data.filter(regex="blockSizes_case.*").columns
            ctr_columns = data.filter(regex="blockSizes_ctr.*").columns
            # Calculate row means for each group with fallback to min_value if all NaN
            row_case_means = data[case_columns].mean(axis=1, skipna=True).fillna(min_value)
            row_ctr_means = data[ctr_columns].mean(axis=1, skipna=True).fillna(min_value)
            # Fill missing values in case columns with their group-specific mean
            for col in case_columns:
                data[col] = data[col].fillna(row_case_means)
            # Fill missing values in control columns with their group-specific mean
            for col in ctr_columns:
                data[col] = data[col].fillna(row_ctr_means)
        # data = data[(data.filter(regex="blockSizes_*").std(axis=1) >= min_std)]
        print("Std filter done!")

        if method == "gamma":
            res = position_gamma(data, processes=processes)
        elif method == "limma":
            res = position_limma(data, features, test_factor, model)

        return res
    
    # Map windows back to positions
    def map_win_2_pos(self, window_data: InputProcessor.data_container, position_data: InputProcessor.data_container, processes=12, sub_window_size = 100, sub_window_step = 100, sub_window_min_diff=0):
        """

        window_data:
            Required Columns: ["chrom", "start", "end"]
        position_data:
            Required Columns: ["chrom", "chromStart", "avg_case", "avg_ctr"]
        
        """

        assert isinstance(window_data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(window_data, ["chrom", "start", "end"])

        assert isinstance(position_data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(position_data, ["chrom", "chromStart", "avg_case", "avg_ctr"])
        
        final_pos = pd.DataFrame()
        wins = window_data.sort_values(["chrom", "start"])
        chroms = position_data["chrom"].unique()
        chroms = [x for x in chroms if "_" not in x]
        
        q_results = deque()

        # TODO MULTI THREAD
        for ch in chroms:
            w = wins[wins["chrom"] == ch].reset_index(drop=True)
            f = position_data[position_data["chrom"] == ch]

            print("Start ", ch)
            for idn in w.index:
                start = w["start"][idn]
                end = w["end"][idn]
                # f["chromEnd"] == f["chromStart"] + 1
                pos_win = f[(f["chromStart"] >= start) & (f["chromStart"]+1 <= end)]
                if pos_win.empty:
                    continue

                valid_window = True
                
                # do sub-window cleaning
                if sub_window_min_diff != 0:
                    sub_start = start
                    sub_end = sub_start + sub_window_size
                    
                    while valid_window and (sub_end <= end):
                        avg_ctr = pos_win[(pos_win["chromStart"] >= sub_start) & (pos_win["chromStart"]+1 <= sub_end)]["avg_ctr"].values
                        avg_case = pos_win[(pos_win["chromStart"] >= sub_start) & (pos_win["chromStart"]+1 <= sub_end)]["avg_case"].values
                        
                        if len(avg_ctr) > 0:
                            avg_ctr = np.mean(avg_ctr)
                            avg_case = np.mean(avg_case)

                            if abs(avg_case - avg_ctr) < sub_window_min_diff:
                                valid_window = False
                        
                        sub_start += sub_window_step
                        sub_end += sub_window_step
                    
                if valid_window:
                    q_results.append(pos_win)
            print("End", ch)
        final_pos = pd.concat(q_results, ignore_index=True)
        final_pos = final_pos.sort_values(["chrom", "chromStart"])
        return final_pos[~final_pos.duplicated()]
        
    def filters(self, data: InputProcessor.data_container, max_q_value=0.05, abs_min_diff=25):
        """Required columns:

        ["q-value", "diff"]
        
        """
        assert isinstance(data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(data, ["q-value", "diff"])

        data = data[data["q-value"] < float(max_q_value)]
        data = data[data["diff"].abs() >= float(abs_min_diff)]

        return data
    
    def generate_DMR(self, significant_position_data: InputProcessor.data_container, position_data: InputProcessor.data_container, min_pos=3, neural_change_limit=7.5, neurl_perc=30, opposite_perc=10):
        """
        
        significant_position_data:
            Required Columns: ["chrom", "chromStart", "diff"]
        position_data:
            Required Columns: ["chrom", "chromStart, "diff"]
        
        """
        neural_change_limit = neural_change_limit/100
        assert isinstance(significant_position_data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(significant_position_data, ["chrom", "chromStart", "diff"])
        assert isinstance(position_data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(position_data, ["chrom", "chromStart", "diff"])
        meth_df = position_data.sort_values(["chrom", "chromStart"]) ## significant_position_data.sort_values(["chrom", "chromStart"]) #
        meth_df.set_index(["chrom", "chromStart"], inplace=True)
        dms_df = significant_position_data.sort_values(["chrom", "chromStart"]) ## position_data.sort_values(["chrom", "chromStart"]) ##
        dms_df.set_index(["chrom", "chromStart"], inplace=True)
        clustered_regions = []
        unclustered_dms = []
        clustered_dms = []
        import tqdm
        for chrom, dms_chr in dms_df.groupby(level=0):  
            data_array = dms_chr.to_numpy()
            chroms = dms_chr.index.get_level_values(0).to_numpy()
            positions = dms_chr.index.get_level_values(1).to_numpy()
            meth_diffs = dms_chr["diff"].to_numpy()
            cluster_start = None
            cluster_end = None
            cluster_sites = []
            cluster_diffs = []
            sign_2 = [0, 0]
            i = -1;
            last_add = i;
            cluster_add_sites = []  
            chrom_mask = (chroms == chrom)
            import time
            start = time.time()
            with tqdm.tqdm(total=len(positions)-1) as pbar:
                # i = 0
                while i<len(positions)-1:
                    i = i + 1;
                    pbar.update(1)
                    if cluster_start is None:
                        cluster_start = positions[i]
                        cluster_end = positions[i]
                        cluster_sites.append(positions[i])
                        if last_add<=i: cluster_add_sites.append(positions[i])
                        cluster_diffs.append(meth_diffs[i])
                        if meth_diffs[i]>0: sign_2[1] = sign_2[1] + 1;
                        else: sign_2[0] = sign_2[0] + 1;
                        continue
                    # Check if the site is within 1kb of the last DMS in cluster
                    if (positions[i] - cluster_end <= 1000) and ((sign_2[1]>0 and sign_2[0]==0 and meth_diffs[i]>0) or (sign_2[0]>0 and sign_2[1]==0 and meth_diffs[i]<0) ):
                        cluster_sites.append(positions[i])
                        if last_add<=i: cluster_add_sites.append(positions[i])
                        cluster_diffs.append(meth_diffs[i])
                        if meth_diffs[i]>0: sign_2[1] = sign_2[1] + 1;
                        else: sign_2[0] = sign_2[0] + 1;
                        #### cluster_end = cluster_sites[0] if len(cluster_sites)<min_pos else cluster_sites[-min_pos] 25/06/2025
                        cluster_end = cluster_sites[-1]
                    else:
                        # Step 2: Process and validate cluster
                        if len(cluster_sites) >= min_pos:
                            trend = np.sign(np.mean(cluster_diffs))
                            if (all(np.sign(cluster_diffs) == trend)):  # Ensure same trend
                                start = time.time()
                                cpg_diff = meth_df.loc[(chrom, slice(cluster_start, cluster_sites[-1])), "diff"]
                                if not cpg_diff.empty:
                                    pos_trend_pct = (cpg_diff>neural_change_limit).sum()*100/cpg_diff.shape[0];
                                    neg_trend_pct = (cpg_diff<-neural_change_limit).sum()*100/cpg_diff.shape[0];
                                    neutral_pct = 100 - pos_trend_pct - neg_trend_pct
                                    start = time.time()
                                    if neutral_pct < neurl_perc and ( (sign_2[1]>0 and neg_trend_pct<opposite_perc) or (sign_2[0]>0 and pos_trend_pct<opposite_perc)):
                                        clustered_regions.append({
                                            "chromosome": chrom,
                                            "start": cluster_start,
                                            "end": cluster_sites[-1],
                                            "num_dms": len(cluster_sites),
                                            "avg_sign_meth": np.mean(cluster_diffs),
                                            "avg_meth": np.mean(cpg_diff),
                                            "std_meth": np.std(cpg_diff),
                                            "pos_trend_pct": pos_trend_pct,
                                            "neg_trend_pct": neg_trend_pct,
                                            "neutral_pct": neutral_pct
                                        })
                                        # cluster_add_sites_set = set(cluster_add_sites) #######################
                                        # mask = chrom_mask & np.array([pos in cluster_add_sites_set for pos in positions]) ################
                                        # clustered_dms.append(mask) ##############################
                                        clustered_dms.append(dms_chr.loc[(chrom, cluster_sites), :])  ##################################
                                        # print("1 ", clustered_dms)
                                        # print("if and big append", time.time() - start) ##############################
                                        # 0, 1, 2, 3, 4; i=5;
                                        # 
                                        if ((sign_2[1]>0 and sign_2[0]==0 and meth_diffs[i]>0) or (sign_2[0]>0 and sign_2[1]==0 and meth_diffs[i]<0) ):
                                            last_add = i
                                            i = i - (min_pos-1)
                                        #if chrom=='chr11':
                                        #   print('ADD', chrom, cluster_sites, cluster_start, cluster_end, len(cluster_sites), neutral_pct, sign_2, neg_trend_pct, pos_trend_pct)
                                    else:
                                        start = time.time()
                                        if len(cluster_add_sites)>0:      
                                            unclustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
                                            # cluster_add_sites_set = set(cluster_add_sites) ###################
                                            # mask = chrom_mask & np.array([pos in cluster_add_sites_set for pos in positions]) ############
                                            # unclustered_dms.append(mask) #################
                                        last_add = i;

                        else: 
                            if len(cluster_add_sites)>0:
                                unclustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
                                # cluster_add_sites_set = set(cluster_add_sites) ###################
                                # mask = chrom_mask & np.array([pos in cluster_add_sites_set for pos in positions]) ##############
                                # unclustered_dms.append(mask) ############
                            last_add = i
                            #if chrom=='chr11':
                            #  print(chrom, cluster_sites, cluster_start, cluster_end, len(cluster_sites), sign_2 )
                        # Reset for new cluster
                        cluster_start = positions[i]
                        cluster_end = positions[i]
                        cluster_sites = [positions[i]]
                        cluster_add_sites = [] ##################################### Just added @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                        if last_add<=i: cluster_add_sites.append(positions[i])
                        cluster_diffs = [meth_diffs[i]]
                        sign_2 = [0, 0]
                        if meth_diffs[i]>0: sign_2[1] = sign_2[1] + 1;
                        else: sign_2[0] = sign_2[0] + 1;
    #  
            if len(cluster_sites) >= min_pos:
                trend = np.sign(np.mean(cluster_diffs))
                if (all(np.sign(cluster_diffs) == trend)):
                    cpg_diff = meth_df.loc[(chrom, slice(cluster_start, cluster_sites[-1])), "diff"]
                    if not cpg_diff.empty:
                        pos_trend_pct = (cpg_diff>neural_change_limit).sum()*100/cpg_diff.shape[0];
                        neg_trend_pct = (cpg_diff<-neural_change_limit).sum()*100/cpg_diff.shape[0];
                        neutral_pct = 100 - pos_trend_pct - neg_trend_pct
                        if neutral_pct < neurl_perc and ( (sign_2[1]>0 and neg_trend_pct<opposite_perc) or (sign_2[0]>0 and pos_trend_pct<opposite_perc)):
                            clustered_regions.append({
                                "chromosome": chrom,
                                "start": cluster_start,
                                "end": cluster_sites[-1],
                                "num_dms": len(cluster_sites),
                                "avg_sign_meth": np.mean(cluster_diffs),
                                "avg_meth": np.mean(cpg_diff),
                                "std_meth": np.std(cpg_diff),
                                "pos_trend_pct": pos_trend_pct,
                                "neg_trend_pct": neg_trend_pct,
                                "neutral_pct": neutral_pct
                            })
                            # cluster_add_sites_set = set(cluster_add_sites) ########
                            # mask = chrom_mask & np.array([pos in cluster_add_sites_set for pos in positions]) #####
                            # clustered_dms.append(mask) ######
                            clustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
                            # print("2", clustered_dms)
                            # 0, 1, 2, 3, 4; i=5;
                            if ((sign_2[1]>0 and sign_2[0]==0 and meth_diffs[i]>0) or (sign_2[0]>0 and sign_2[1]==0 and meth_diffs[i]<0) ):
                                last_add = i
                                i = i - (min_pos-1)
                        else:
                            if len(cluster_add_sites)>0:
                                # cluster_add_sites_set = set(cluster_add_sites) ##################
                                # mask = chrom_mask & np.array([pos in cluster_add_sites_set for pos in positions]) #########
                                # unclustered_dms.append(mask) ##############
                                unclustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
                            last_add = i;
            else:
                if len(cluster_add_sites)>0:
                    # cluster_add_sites_set = set(cluster_add_sites) ################
                    # mask = chrom_mask & np.array([pos in cluster_add_sites_set for pos in positions]) ##############
                    # unclustered_dms.append(mask) ##############
                    unclustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
        cluster_df = pd.DataFrame(clustered_regions)
        # TODO convert list of lists to dataframe
        unclustered_dms_df = pd.concat(unclustered_dms) if unclustered_dms else pd.DataFrame()
        # print(unclustered_dms)
        clustered_dms_df = pd.concat(clustered_dms) if clustered_dms else pd.DataFrame()
        return cluster_df, unclustered_dms_df.reset_index(), clustered_dms_df.reset_index()

    def map_positions_to_genes(self, positions: InputProcessor.data_container, gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], min_pos_diff=0, bed_file="CpG_gencodev42ccrenb_repeat_epic1v2hm450.bed", gtf_file="outfile_w_hm450.bed"):#"gencode.v41.chr_patch_hapl_scaff.annotation.gtf"):
        """
        
        Required columns: ["chrom", "chromStart", "chromEnd", "diff"]

        """
        
        assert isinstance(positions, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(positions,["chrom", "chromStart",  "diff"]) #["chrom", "chromStart", "strand", "diff"])

        if isinstance(gene_regions, str):
            gene_regions = [gene_regions]

        return position_based_gene(positions, gene_regions, min_pos_diff, bed_file, gtf_file)

    def __calculate_region_distance(self, start1, end1, start2, end2):
      if end1 >= start2 and end2 >= start1:
          return 0  # Regions overlap
      # Calculate distance if no overlap
      if end1 < start2:
          return start2 - end1  # Region 1 is before Region 2
      else:
          return start1 - end2  # Region 2 is before Region 1
    @staticmethod
    def parse_annotation(annotation):
      """
      Parse the annotation field to categorize into gene intron, gene exon, IG, and ENCODE types.
      """
      categories = set()
      encode_types = set()
      repeat = set()
      cpg_epic = set()
      enhance_list = dict();
      gene_list = dict();
      nb_suf = 'nb100'
      #EPIC not consider
      annotations = annotation.split('/')
      for entry in annotations:
          if ':' in entry:
              key, value = entry.split(':', 1)
              if key == 'IG':
                  categories.add('IG')
              elif key == 'ENCODE' or key == 'CCRE':
                  encode_type = value.split('@')[1] if '@' in value else value
                  if encode_type[-len(nb_suf):]==nb_suf:
                     encode_type = encode_type.split('_')[0]+'_nb'
                  encode_types.add(encode_type)
                  en_id = value.split('@')[0] if '@' in value else value
                  if en_id not in gene_list:
                     gene_list[ en_id ] = set()
                  gene_list[ en_id ].add(encode_type)
              elif key == 'REPEAT':
                  repeat.add( value);
              elif key=='EPIC' or key=='EPICv2' or key=='HM450':
                  cpg_epic.add( value)
              else:
                  if 'intron' in value.lower():
                      categories.add('Gene_Intron')
                      if 'intron' not in gene_list: gene_list[ 'intron'] = set()
                      gene_list[ 'intron'].add(key);
                  elif 'exon' in value.lower():
                      categories.add('Gene_Exon')
                      if 'exon' not in gene_list: gene_list[ 'exon'] = set()
                      gene_list[ 'exon'].add(key);
      # Ensure IG is not counted if others exist
      if ('Gene_Intron' in categories or 'Gene_Exon' in categories or encode_types):
          categories.discard('IG')
      if 'enhD' in encode_types or 'enhP' in encode_types or 'prom' in encode_types or 'K4m3' in encode_types:
         categories.discard('Gene_Intron')
      if 'enhD' in encode_types:
         encode_types.discard('enhD');
         encode_types.add('enh')
      if 'enhP' in encode_types:
         encode_types.discard('enhP');
         encode_types.add('enh')
      merg_enh = []
      for _et in encode_types:
         if '_' in _et and _et.split('_')[0] in ['enhD', 'enhP']:
            merg_enh.append( _et );
      for _et in merg_enh:
         encode_types.discard( _et );
         encode_types.add('enh'+'_'+_et.split('_')[1] )
      return categories, encode_types, repeat, cpg_epic, gene_list
    def process_regions(self, region_file: InputProcessor.data_container, annotation_file:str = "CpG_gencodev42ccrenb_repeat_epic1v2hm450.bed", gene_bed_file:str = "gencode.v42.chr_patch_hapl_scaff.annotation.genes.bed", ccre_file:str ="encodeCcreCombined.bed" ):
        """
        Reads region file and annotation file, finds matches, and counts occurrences.
        """
        assert isinstance(region_file, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(region_file, ['chrom', 'start', 'end'])
        # Read region file
        #regions_df = pd.read_csv(region_file, sep=',', header=0, names=['chr', 'start', 'end'], usecols=[0,1,2])
        regions_df = region_file
        if (('chromosome' in regions_df.columns) or ('chrom' in regions_df.columns)) and ('chr' not in regions_df.columns):
            regions_df = regions_df.rename(columns={'chromosome':'chr', 'chrom':"chr"})
        gene_bed = pd.read_csv(gene_bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'gene', 'strand'], usecols=[0,1,2,3,5])
        ccre_bed = pd.read_csv(ccre_file, sep='\t', header=None, names=['chr', 'start', 'end', 'ccre'], usecols=[0,1,2,13])
        # Read annotation file
        annotation_df = pd.read_csv(annotation_file, sep='\t', header=None, names=['chr', 'start', 'end', 'id', 'annotation'])
        exon_genes = defaultdict(int)
        intro_genes  = defaultdict(int)
        prom_genes = defaultdict(int)
        enhP_genes = defaultdict(int)
        enhD_genes = defaultdict(int)
        # Iterate through regions and match annotations
        for _, region in regions_df.iterrows():
            matched_annotations = annotation_df[
            (annotation_df['chr'] == region['chr']) &
            (annotation_df['start'] <= region['end']) &
            (annotation_df['end'] >= region['start'])
            ]
            for _, annotation in matched_annotations.iterrows():
                categories, encode_types, repeat, cpg_epic, gene_list = self.parse_annotation(annotation['annotation'])
                for _ge in gene_list:
                    for _de in gene_list[ _ge ]:
                        if _de=='enhD': t_dict = enhD_genes
                        elif _de=='enhP': t_dict = enhP_genes
                        elif _de=='prom': t_dict = prom_genes
                        else:
                            if _ge=='exon': t_dict = exon_genes
                            elif _ge=='intron': t_dict = intro_genes
                            else:
                                if not (len(_de)>3 and (_de[-3:]=='_nb' or _de=='K4m3') ):
                                    print("does not support", _ge, _de);
                                continue;
                        if _de in ['enhD', 'enhP', 'prom']:
                            t_dict[_ge] += 1;
                        else:
                            t_dict[_de ] += 1
        promF_genes = defaultdict(int)
        enhPF_genes = defaultdict(int)
        enhDF_genes = defaultdict(int)
        enhd_thr = 500000
        enhp_thr = 50000
        prom_thr = 2000
        thresholds = [enhd_thr, enhp_thr, prom_thr]
        gene_list = [enhD_genes, enhP_genes, prom_genes]
        res_list = [enhDF_genes, enhPF_genes, promF_genes]
        for _gi in range(len(thresholds)):
            t_thr = thresholds[ _gi ]
            t_ccre = gene_list[ _gi ]
            t_genes = res_list[ _gi ]
            for t_e in t_ccre:
                ccre_region = ccre_bed[ ccre_bed['ccre']==t_e ]
                if not ccre_region.shape[0]==1:
                    print("wrong id", t_e, ccre_region.shape, ccre_region)
                    continue;
                res_gene_def = gene_bed[ (gene_bed['chr']==ccre_region['chr'].iloc[0]) & ( gene_bed['start']-t_thr <=ccre_region['end'].iloc[0] ) & ( gene_bed['end']+t_thr>=ccre_region['start'].iloc[0]  )   ]
                if _gi==0: this_t = "enhD"
                elif _gi==1: this_t = "enhP"
                else: this_t = "Prom"
                outpt_list = [t_e+'\t'+this_t]
                if _gi in [0, 1]:
                    for _, g_row in res_gene_def.iterrows():
                        t_genes[ g_row['gene'] ] += 1;
                        outpt_list.append("\n\tANNOg:"+g_row['gene']+':'+str(self.__calculate_region_distance(ccre_region['start'].iloc[0], ccre_region['end'].iloc[0], g_row['start'],  g_row['end']) ))
                    print("ANNOe:"+''.join( outpt_list ) );
                else:
                    for _, g_row in res_gene_def.iterrows():
                        if g_row['strand']=='+':
                            if g_row['start']>ccre_region['start'].iloc[0]:
                                t_genes[ g_row['gene'] ] += 1;
                                outpt_list.append("\n\tANNOg:"+g_row['gene']+':'+str( g_row['start']-ccre_region['start'].iloc[0] ) )
                        else:
                            if g_row['end'] < ccre_region['start'].iloc[0]:
                                t_genes[ g_row['gene'] ] += 1;
                                outpt_list.append("\n\tANNOg:"+g_row['gene']+':'+str( g_row['end'] - ccre_region['start'].iloc[0] ) )
                    print("ANNOe:"+''.join( outpt_list ) );
        ccre_genes = defaultdict(int)
        for c_g_list in [prom_genes, enhP_genes, enhD_genes]:
            for _c_k in c_g_list:
                ccre_genes[ _c_k] += c_g_list[ _c_k ]
        final_list = set();
        rem_list = defaultdict(int)
        for c_g_list in [prom_genes, enhP_genes, enhD_genes, exon_genes, intro_genes]:
            for _c_k in c_g_list:
                if len(_c_k)>4 and _c_k[:4] in ['ENSG', 'LINC']:
                    rem_list[ _c_k] += 1;
                elif len(_c_k.split('-AS'))>1:
                    rem_list[ _c_k] += 1;
                else:
                    final_list.add( _c_k )
        rem_list = list(rem_list.keys())
        final_list = list(final_list)
        # rem_list = pd.DataFrame.from_dict(rem_list)
        # final_list = pd.Series(final_list)
        return rem_list, final_list
