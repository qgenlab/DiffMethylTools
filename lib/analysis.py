from collections import deque
import multiprocessing
from queue import Queue
from typing import Callable
import numpy as np
import pandas as pd

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
            col = "^" + col + "$"
            if not any(re.match(col, c) for c in columns):
                raise AssertionError(f"Required column pattern '{col}' not found in columns: {columns}")
    
    def assert_one_of_column_pairs(self, df, column_pairs):
        columns = df.columns
        for pair in column_pairs:
            if all(any(re.match(col, c) for c in columns) for col in pair):
                return
        raise AssertionError(f"None of the required column pairs {column_pairs} found in columns: {columns}")


    def merge_tables(self, min_cov = 10, cov_percentile = 1.0, min_samp_ctr = 4, min_samp_case = 3, case_data: InputProcessor.data_container = None, ctr_data: InputProcessor.data_container = None):
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
        d["ctr"] = [pl.from_pandas(x) for x in ctr_data]
        d["case"] = [pl.from_pandas(x) for x in case_data]

        for key in ["ctr", "case"]:
            for i in range(len(d[key])):
                df = d[key][i]
                for col in df.columns:
                    if "blockSizes" in col or "coverage" in col or "positive" in col or "negative" in col:
                        df = df.with_columns(pl.col(col).cast(pl.Float64, strict=False).alias(col))

                # for positive, negative setup, calcualte coverage and methylation percentage
                positive_col = next((col for col in df.columns if col.startswith("positive_")), None)
                negative_col = next((col for col in df.columns if col.startswith("negative_")), None)
                if positive_col and negative_col:
                    name = re.match(r"positive_(.*)", positive_col).group(1)
                    if f"negative_{name}" == negative_col:
                        df = df.with_columns((pl.col(positive_col) + pl.col(negative_col)).alias(f"coverage_{key}_{name}"))
                        df = df.with_columns((pl.col(positive_col) / (pl.col(positive_col) + pl.col(negative_col)) * 100).alias(f"blockSizes_{key}_{name}"))
                        df = df.drop([positive_col, negative_col])
                d[key][i] = df

                # select only necessary columns for merge
                d[key][i] = df.select(MERGE_LIST + [col for col in df.columns if "blockSizes" in col or "coverage" in col])
                
        ctr = d["ctr"][0]
        case = d["case"][0]
        del(d["ctr"][0])
        del(d["case"][0])
        for e in range(0,len(d["ctr"])):
            ctr = ctr.join(d["ctr"][e], on=MERGE_LIST, how='full', coalesce=True)
        del(d["ctr"])
        for e in range(0,len(d["case"])):
            case = case.join(d["case"][e], on=MERGE_LIST, how='full', coalesce=True)
        del(d)
        cov_columns = cs.starts_with("coverage_")
        case1 = case.with_columns(sum_case=pl.sum_horizontal(case.select((cov_columns > min_cov) & (cov_columns < cov_columns.replace(0,None).quantile(cov_percentile)))))
        ctr1 = ctr.with_columns(sum_ctr=pl.sum_horizontal(ctr.select((cov_columns > min_cov) & (cov_columns < cov_columns.replace(0,None).quantile(cov_percentile)))))
        sub = case1.select(MERGE_LIST + ["sum_case"]).join(ctr1.select(MERGE_LIST + ["sum_ctr"]), on=MERGE_LIST, how='inner', suffix ="_2", coalesce=True)
        sub = sub.filter((pl.col("sum_case") >= min_samp_case) & (pl.col("sum_ctr") >= min_samp_ctr))
        case = case.join(sub[tuple(MERGE_LIST)], on=MERGE_LIST, how='inner', coalesce=True)
        ctr = ctr.join(sub[tuple(MERGE_LIST)], on=MERGE_LIST, how='inner', coalesce=True)
        case = case.to_pandas()
        ctr = ctr.to_pandas()
        avg_case = case.filter(like="blockSizes").mean(axis=1, numeric_only=True)
        avg_ctr = ctr.filter(like="blockSizes").mean(axis=1, numeric_only=True)
        case["avg_case"] = avg_case
        ctr["avg_ctr"] = avg_ctr
        f = pd.merge(case,ctr,on=MERGE_LIST, how='inner')
        f["diff"] = f["avg_case"] - f["avg_ctr"]
        self._merged_result = f.sort_values(by=MERGE_LIST).reset_index(drop=True)
        return self._merged_result

    # helper windowed p-value function
    def helper_window_based(self, chrom, f, step, window, min_nbr, statistical_test):
        chrom_starts = f["chromStart"].values
        avg_ctr = f["avg_ctr"].values
        avg_case = f["avg_case"].values
        diffs = (f["avg_case"] - f["avg_ctr"]).values

        q = Queue()
        
        init = (chrom_starts[0] // step - 1) * step
        end = (chrom_starts[-1] // step + 1) * step
        print("Start", chrom)
        while init < end:
            mask = (chrom_starts >= init) & (chrom_starts < init + window)
            output_len = np.sum(mask)
            
            if output_len >= min_nbr:
                m = diffs[mask].mean()
                p_val = statistical_test(avg_ctr[mask], avg_case[mask], axis=0)[1]
                q.put((chrom, init, init + window, output_len, p_val, m))
            init += step


        print("End", chrom)
        q.put(None)
        return list(iter(q.get, None))
        
    # main windowed p-value function
    def window_based(self, data: InputProcessor.data_container, statistical_test : Callable, window = 1000, step = 500, min_nbr_per_win = 5, processes=12, min_std=0.1):
        """Required columns:

        ["chrom", "chromStart", "blockSizes_case.*", "blockSizes_ctr.*", "avg_case", "avg_ctr"]
        
        """
        assert isinstance(data, pd.DataFrame), "List input not acceptable for this function."
        assert hasattr(statistical_test, "__module__") and statistical_test.__module__.startswith("scipy.stats"), "statistical_test must be from scipy.stats package."
        self.assert_required_columns(data, ["chrom", "chromStart", "blockSizes_case.*", "blockSizes_ctr.*", "avg_case", "avg_ctr"])
        assert statistical_test in [mannwhitneyu, ttest_ind, ranksums, ks_2samp, median_test, brunnermunzel], "statistical_test must be one of the following scipy.stats functions: [mannwhitneyu, ttest_ind, ranksums, ks_2samp, median_test, brunnermunzel]"
        
        # filter standard deviation
        data = data[(data.filter(regex="blockSizes_case.*").std(axis=1) > min_std) & (data.filter(regex="blockSizes_ctr.*").std(axis=1) > min_std)]

        chroms = data["chrom"].unique()

        chroms = [x for x in chroms if "_" not in x]

        q_results = []

        pool = multiprocessing.Pool(processes=processes)

        for ch in chroms:
            f = data[data["chrom"] == ch]
            q = pool.apply_async(self.helper_window_based, args=(ch, f, step, window, min_nbr_per_win, statistical_test))
            q_results.extend([q])

        l = []

        for q in q_results:
            l.extend(q.get())

        del q_results
            
        pool.close()
        pool.join()

        return pd.DataFrame(l, columns = ["chrom", "start", "end", "nbr", "p-val", "diff"])
    
    def position_based(self, data: InputProcessor.data_container, method="limma", features=None, test_factor="Group", processes=22, model="eBayes", min_std=0.1, fill_na:bool=True):
        """
    
        Required columns: ["chrom", "chromStart", "blockSizes_case.*", "blockSizes_ctr.*"]

        Columns in the features csv must match the name given in the blockSizes columns (for example, use the name "test" for "blockSizes_ctr_test"). 

        """
        assert method in ["gamma", "limma"], "Select a valid method parameter. Options are: [""gamma"", ""limma""]"
        assert isinstance(data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(data, ["chrom", "chromStart", "blockSizes_case.*", "blockSizes_ctr.*"])

        # filter standard deviation
        data = data[(data.filter(regex="blockSizes_case.*").std(axis=1) > min_std) & (data.filter(regex="blockSizes_ctr.*").std(axis=1) > min_std)]
        # fill na with mean of blocksizes_case and blocksizes_ctr separately
        if fill_na:
            # Fill NA values for blockSizes_case columns by row average
            case_columns = data.filter(regex="blockSizes_case.*").columns
            data[case_columns] = data[case_columns].apply(lambda row: row.fillna(row.mean()), axis=1)

            # Fill NA values for blockSizes_ctr columns by row average
            ctr_columns = data.filter(regex="blockSizes_ctr.*").columns
            data[ctr_columns] = data[ctr_columns].apply(lambda row: row.fillna(row.mean()), axis=1)


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

            print("Start", ch)
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

    def generate_q_values(self, data: InputProcessor.data_container, method: str = "fdr_bh"):
        """Required columns:
        
        ["p-val"]

        """
        assert isinstance(data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(data, ["p-val"])

        _ = multitest.multipletests(data["p-val"], method=method)
        df = pl.DataFrame({'q-value':_[1].tolist()}).to_pandas()
        return pd.concat([data, df], axis=1)
        
    def filters(self, data: InputProcessor.data_container, max_q_value=0.05, abs_min_diff=25):
        """Required columns:

        ["q-value", "diff"]
        
        """
        assert isinstance(data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(data, ["q-value", "diff"])

        data = data[data["q-value"] < max_q_value]
        data = data[data["diff"].abs() >= abs_min_diff]

        return data
    
    def filter_diffs(self, min_diff = 25):
        if self.using_positional_pipeline:
            self._positional = self._positional[self._positional["diff"].abs() >= min_diff].reset_index(drop=True)
        else:
            self._res = self._res[self._res["diff"].abs() >= min_diff].reset_index(drop=True)

    
    def generate_DMR(self, significant_position_data: InputProcessor.data_container, position_data: InputProcessor.data_container, min_pos=3, neutral_change_limit=7.5, neutral_perc=30, opposite_perc=10):
        """
        
        significant_position_data:
            Required Columns: ["chrom", "chromStart", "diff"]
        position_data:
            Required Columns: ["chrom", "chromStart, "diff"]
        
        """
        assert isinstance(significant_position_data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(significant_position_data, ["chrom", "chromStart", "diff"])
        assert isinstance(position_data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(position_data, ["chrom", "chromStart", "diff"])
        

        meth_df = significant_position_data.sort_values(["chrom", "chromStart"])
        meth_df.set_index(["chrom", "chromStart"], inplace=True)
        dms_df = position_data.sort_values(["chrom", "chromStart"])
        dms_df.set_index(["chrom", "chromStart"], inplace=True)
        clustered_regions = []
        unclustered_dms = []
        clustered_dms = []
        for chrom, dms_chr in dms_df.groupby(level=0):  
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
            while i<len(positions)-1:
                i = i + 1;
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
                    cluster_end = cluster_sites[0] if len(cluster_sites)<min_pos else cluster_sites[-min_pos]
                else:
                    # Step 2: Process and validate cluster
                    if len(cluster_sites) >= min_pos:
                        trend = np.sign(np.mean(cluster_diffs))
                        if (all(np.sign(cluster_diffs) == trend)):  # Ensure same trend
                            cpg_diff = meth_df.loc[(chrom, slice(cluster_start, cluster_sites[-1])), "diff"]
                            if not cpg_diff.empty:
                                pos_trend_pct = (cpg_diff>neutral_change_limit).sum()*100/cpg_diff.shape[0];
                                neg_trend_pct = (cpg_diff<-neutral_change_limit).sum()*100/cpg_diff.shape[0];
                                neutral_pct = 100 - pos_trend_pct - neg_trend_pct
                                if neutral_pct < neutral_perc and ( (sign_2[1]>0 and neg_trend_pct<opposite_perc) or (sign_2[0]>0 and pos_trend_pct<opposite_perc)):
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
                                    clustered_dms.append(dms_chr.loc[(chrom, cluster_sites), :])
                                    # 0, 1, 2, 3, 4; i=5;
                                    # 
                                    if ((sign_2[1]>0 and sign_2[0]==0 and meth_diffs[i]>0) or (sign_2[0]>0 and sign_2[1]==0 and meth_diffs[i]<0) ):
                                        last_add = i
                                        i = i - (min_pos-1)
                                    #if chrom=='chr11':
                                    #   print('ADD', chrom, cluster_sites, cluster_start, cluster_end, len(cluster_sites), neutral_pct, sign_2, neg_trend_pct, pos_trend_pct)
                                else:
                                    if len(cluster_add_sites)>0:
                                        unclustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
                                    last_add = i;
                    else: 
                        if len(cluster_add_sites)>0:
                            unclustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
                        last_add = i;
                        #if chrom=='chr11':
                        #  print(chrom, cluster_sites, cluster_start, cluster_end, len(cluster_sites), sign_2 )
                    # Reset for new cluster
                    cluster_start = positions[i]
                    cluster_end = positions[i]
                    cluster_sites = [positions[i]]
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
                    pos_trend_pct = (cpg_diff>neutral_change_limit).sum()*100/cpg_diff.shape[0];
                    neg_trend_pct = (cpg_diff<-neutral_change_limit).sum()*100/cpg_diff.shape[0];
                    neutral_pct = 100 - pos_trend_pct - neg_trend_pct
                    if neutral_pct < neutral_perc and ( (sign_2[1]>0 and neg_trend_pct<opposite_perc) or (sign_2[0]>0 and pos_trend_pct<opposite_perc)):
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
                        clustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
                        # 0, 1, 2, 3, 4; i=5;
                        if ((sign_2[1]>0 and sign_2[0]==0 and meth_diffs[i]>0) or (sign_2[0]>0 and sign_2[1]==0 and meth_diffs[i]<0) ):
                            last_add = i
                            i = i - (min_pos-1)
                    else:
                        if len(cluster_add_sites)>0:
                            unclustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
                        last_add = i;
        else:
            if len(cluster_add_sites)>0:
                unclustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])

        cluster_df = pd.DataFrame(clustered_regions)
        unclustered_dms_df = pd.concat(unclustered_dms) if unclustered_dms else pd.DataFrame()
        clustered_dms_df = pd.concat(clustered_dms) if clustered_dms else pd.DataFrame()
        return cluster_df, unclustered_dms_df.reset_index(), clustered_dms_df.reset_index()

    def map_positions_to_genes(self, positions: InputProcessor.data_container, gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], min_pos_diff=0, bed_file="outfile_w_hm450.bed", gtf_file="gencode.v41.chr_patch_hapl_scaff.annotation.gtf"):
        """
        
        Required columns: ["chrom", "chromStart", "chromEnd", "diff"]

        """
        
        assert isinstance(positions, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(positions, ["chrom", "chromStart", "strand", "diff"])

        if isinstance(gene_regions, str):
            gene_regions = [gene_regions]

        return position_based_gene(positions, gene_regions, min_pos_diff)

    def map_windows_to_genes(self, windows: InputProcessor.data_container, gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], min_pos_diff=0, bed_file="outfile_w_hm450.bed", gtf_file="gencode.v41.chr_patch_hapl_scaff.annotation.gtf", enhd_thr = 500000, enhp_thr = 50000, prom_thr = 2000, processes=12):
        """
        
        Required columns: ["chrom", "start", "end", "diff"]
        
        """
        
        assert isinstance(windows, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(windows, ["chrom", "start", "end", "diff"])

        if isinstance(gene_regions, str):
            gene_regions = [gene_regions]

        return window_based_gene(windows, gene_regions, min_pos_diff)