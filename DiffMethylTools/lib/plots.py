from collections import deque
import math
from typing import Optional, Union
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
import pandas as pd
import polars as pl
import numpy as np
import seaborn
from .input_processor import InputProcessor
from .analysis import Analysis
from collections import defaultdict, Counter
import random

class Plots():
    def __init__(self):
        pass

    def assert_required_columns(self, df, required_columns):
        Analysis().assert_required_columns(df, required_columns)
    
    def assert_one_of_column_pairs(self, df, column_pairs):
        Analysis().assert_one_of_column_pairs(df, column_pairs)

    def volcano_plot(self, data: InputProcessor.data_container, name : str, threshold : Optional[float] = 0.05, line: Optional[float] = None, x_range: tuple[int, int] = (-100, 100), y_max: int = None, title: str = None, x_label: str = None, y_label: str = None):
        assert isinstance(data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(data, ["diff", "q-value"])
        
        fig, ax = plt.subplots(1,1, figsize=(8, 6))
        fig.subplots_adjust(hspace=0.3)
        df = data
        df["neg_log10_fdr"] = -1 * np.log10(df["q-value"])
        if threshold is not None:
            if line is not None:
                significant = df[(df["q-value"] < threshold) & (df["diff"].abs() > line)]
                insignificant = df[(df["q-value"] >= threshold) | (df["diff"].abs() <= line)]
            else:
                significant = df[df["q-value"] < threshold]
                insignificant = df[df["q-value"] >= threshold]
            ax.scatter(significant["diff"].to_numpy(), significant["neg_log10_fdr"].to_numpy(), s=5,zorder=1)
            ax.scatter(insignificant["diff"].to_numpy(), insignificant["neg_log10_fdr"].to_numpy(), s=5, c="gray", zorder=0)
        else:
            ax.scatter(df["diff"].to_numpy(), df["neg_log10_fdr"].to_numpy(), s=5)

        if x_range[0] is not None:
            left_edge = x_range[0]
        else:
            left_edge = df["diff"].min()
            
        if x_range[1] is not None:
            right_edge = x_range[1]
        else:
            right_edge = df["diff"].max()

        if y_max is not None:
            ax.set_ylim(0, y_max)

        ax.set_xlim(left_edge, right_edge)
        if threshold is not None:
            ax.axhline(y = np.log10([threshold])[0] * -1, color = 'maroon', linestyle = '--', label=f'FDR = {threshold}')
        if line is not None:
            ax.axvline(x = -line, color = 'maroon', linestyle = '--')
            ax.axvline(x = line, color = 'maroon', linestyle = '--')
        if title is not None:
            ax.set_title(title, fontsize=20)
        else:
            ax.set_title('Volcano Plot', fontsize=20)
        if x_label is not None:
            ax.set_xlabel(x_label, fontsize=15)
        else:
            ax.set_xlabel("Average case - Average control", fontsize=15)
        if y_label is not None:
            ax.set_ylabel(y_label, fontsize=15)
        else:
            ax.set_ylabel("- log10(FDR)", fontsize=15)

        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
        ax.legend(ax.get_legend_handles_labels()[0], ax.get_legend_handles_labels()[1], fontsize=13)
        plt.tight_layout()
        plt.savefig(name, dpi=300)
    
    def manhattan_plot(self, data: InputProcessor.data_container, name: str, threshold: Optional[float] = 0.05, title: str = None, x_label: str = None, y_label: str = None):
        assert isinstance(data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(data, ["chrom", "q-value"])
        self.assert_one_of_column_pairs(data, [("start",), ("chromStart",)])
        if "start" in data.columns:
            start_name = "start"
        else:
            start_name = "chromStart"

        fig, ax = plt.subplots(1, 1, figsize=(15, 10))
        df = data
        df["neg_log10_fdr"] = -1 * np.log10(df["q-value"])

        if "diff" in df.columns:
            print("Generating signed values")
            df["sign_neg_log10_fdr"] = df["neg_log10_fdr"] * np.sign(df["diff"])
        else:
            df["sign_neg_log10_fdr"] = df["neg_log10_fdr"]
        ############ if "diff" in list(df.columns): print("Generating signed values"); df["sign_neg_log10_fdr"] = df["neg_log10_fdr"] * np.sign(df['diff'])
        # Sort chromosomes correctly
        df["chrom"] = df["chrom"].str.extract(r'(\d+|X|Y)')[0]
        df["chrom"] = pd.Categorical(df["chrom"], categories=[str(i) for i in range(1, 23)] + ["X", "Y"], ordered=True)
        df = df.sort_values(["chrom", start_name])
        
        df_grouped = df.groupby("chrom", observed=True)

        colors = ['black', 'gray']
        x_labels = []
        x_labels_pos = []
        current_position = 0

        for num, (grouped_name, group) in enumerate(df_grouped):
            group = group.sort_values(start_name)
            group["plot_position"] = group[start_name] + current_position
            ax.scatter(group['plot_position'], group["sign_neg_log10_fdr"], color=colors[num % len(colors)], s=10)
            x_labels.append(grouped_name)
            x_labels_pos.append((group['plot_position'].iloc[-1] + group['plot_position'].iloc[0]) / 2)
            current_position = group['plot_position'].iloc[-1] + 1

        if threshold is not None:
            ax.axhline(y=-np.log10(threshold), color='maroon', linestyle='--')
            ax.axhline(y=np.log10(threshold), color='maroon', linestyle='--')
            ax.axhline(y=0, color='blue')

        ax.set_xticks(x_labels_pos)
        ax.set_xticklabels(x_labels)

        if title is not None:
            ax.set_title(title, fontsize=25)
        else:
            ax.set_title('Manhattan Plot', fontsize=25)
        if x_label is not None:
            ax.set_xlabel(x_label, fontsize=27) # 15
        else:
            ax.set_xlabel('Chromosome', fontsize=27) # 15
        # if y_label is not None:
        #     ax.set_ylabel(y_label, fontsize=17) # 15
        # else:
        #     ax.set_ylabel('- log10(FDR)', fontsize=17) # 15
        ax.set_ylabel("")
        ax.text(-0.05, 0.80, '-log10(FDR)', rotation=90, va='top', ha='center', fontsize=25, transform=ax.transAxes)
        ax.text(-0.05, 0.20, 'log10(FDR)', rotation=90, va='bottom', ha='center', fontsize=25, transform=ax.transAxes)
        max_abs = df["sign_neg_log10_fdr"].abs().max()
        ax.set_ylim([-max_abs, max_abs])

        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=22)
        plt.tight_layout()
        plt.savefig(name, dpi=600)

    def coverage_plot(self, ctr_data: Union[InputProcessor.data_container, list[InputProcessor.data_container]], case_data: Union[InputProcessor.data_container, list[InputProcessor.data_container]], name : str, cov_min : int = -1, cov_max : int = -1, cov_max_percentile : float = -1, bins:int = 20, title: str = None, x_label: str = None, y_label: str = None):
        # assert that data is a list of dataframes
        # TODO maybe change this to allow full file
        assert isinstance(ctr_data, list), "Input data must be from a list. (The same as merge_data)"
        for data in [ctr_data, case_data]:
            for df in data:
                self.assert_one_of_column_pairs(df, (["coverage.*"], ["positive_*", "negative_*"]))
            
        
        def __individual_plot(ax, data, title, cov_min=-1, cov_max=-1):
            filter = np.array([1] * len(data))
            if cov_max_percentile != -1:
                cov_max = np.percentile(data, cov_max_percentile)
            if cov_min != -1:
                filter &= (data > cov_min)
            if cov_max != -1:
                filter &= (data < cov_max)
            
            plot_data = data[np.where(filter)]
            ax.hist(plot_data, bins=bins)
            ax.set_title(title, fontsize=25)
            if x_label is not None:
                ax.set_xlabel(x_label, fontsize=20) # 15
            else:    
                ax.set_xlabel("Read Coverage", fontsize=20) # 15
            if y_label is not None:
                ax.set_ylabel(y_label, fontsize=20) # 15
            else:
                ax.set_ylabel("Count", fontsize=20) # 15
            ax.tick_params(axis='x', labelsize=10)
            ax.tick_params(axis='y', labelsize=10)
            # ax.set_aspect('auto')

        for i, df_list in enumerate([ctr_data, case_data]):
            ctr_case = "ctr" if i == 0 else "case"
            for df in df_list:
                # find if coverage_* is in columns already
                if any(df.columns.str.contains(f"^(?:coverage_{ctr_case}_.*)$", regex=True)):
                    continue
                # make coverage_ctr* = negavtive_ctr* + positive_ctr*
                print(df)
                tmp_name = df.filter(regex="^(negative)_.*").columns[0]
                tmp_name = tmp_name[tmp_name.find("_")+1:]
                df[f"coverage_{ctr_case}_{tmp_name}"] = df.filter(regex="^((negative)|(positive))_.*$").sum(axis=1)



        control_samples = [pl.from_pandas(df).select(pl.col(f"^coverage_ctr_.*$")).to_numpy().flatten() 
                           for df in ctr_data]
        case_samples = [pl.from_pandas(df).select(pl.col(f"^coverage_case_.*$")).to_numpy().flatten() 
                        for df in case_data]
        
        
        
        num_ctrl, num_case = len(control_samples), len(case_samples)
        max_samples = max(num_ctrl, num_case)
        n = int(math.sqrt(max_samples))
        m = math.ceil(max_samples / n) - n
        
        rows = n + m
        cols = n * 2

        subplot_size = 5  # Size of each subplot in inches

        # Calculate the aspect ratio of the grid
        aspect_ratio = cols / rows

        # Adjust the figure size to balance the aspect ratio
        if aspect_ratio > 1:  # More columns than rows
            fig_width = subplot_size * cols
            fig_height = subplot_size * rows * aspect_ratio/2
        else:  # More rows than columns
            fig_width = subplot_size * cols / aspect_ratio
            fig_height = subplot_size * rows
        print(rows, cols)
        print((fig_width, fig_height))
        fig, axes = plt.subplots(rows, cols, figsize=(fig_width, fig_height), squeeze=False)
        fig.subplots_adjust(hspace=0.4, wspace=0.4) 
    
        
        control_idx = 0
        case_idx = 0

        for i in range(n + m):
            for j in range(n*2):
                if j < n and control_idx < len(control_samples):
                    __individual_plot(axes[i, j], control_samples[control_idx], f"Control Sample {control_idx+1}", cov_min, cov_max)
                    control_idx += 1
                elif j >= n and case_idx < len(case_samples):
                    __individual_plot(axes[i, j], case_samples[case_idx], f"Case Sample {case_idx + 1}", cov_min, cov_max)
                    case_idx += 1
                else:
                    axes[i, j].axis('off')

        if title is not None:
            plt.suptitle(title, fontsize=30)
        else:
            plt.suptitle("Coverage Histograms", fontsize=30)

        plt.tight_layout()
        plt.savefig(name, dpi=300)

    def graph_gene_regions(self, gene_data: InputProcessor.data_container, ccre_data: InputProcessor.data_container, name: str, gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], intron_cutoff: int = -1, exon_cutoff: int = -1, upstream_cutoff: int = -1, CCRE_cutoff: int = -1,prom_cutoff:int = -1, title: str = None, x_label: str = None, intron_y_label: str = None, exon_y_label: str = None, upstream_y_label: str = None, CCRE_y_label: str = None, prom_y_label: str = None):
        # check if gene_Data and ccre_data are not list
        """

        Required Columns:

        gene_data: ["intron", "intron_diff", "exon", "exon_diff", "upstream", "upstream_diff"] depending on the gene_regions selected.
        ccre_data: ["CCRE", "CCRE_diff"] if "CCRE" is in gene_regions.

        """
        
        assert isinstance(gene_data, pd.DataFrame), "List input not acceptable for this function."
        assert isinstance(ccre_data, pd.DataFrame), "List input not acceptable for this function."


        if "intron" in gene_regions:
            self.assert_required_columns(gene_data, ["intron", "intron_diff"])
        if "exon" in gene_regions:
            self.assert_required_columns(gene_data, ["exon", "exon_diff"])
        if "upsteam" in gene_regions:
            self.assert_required_columns(gene_data, ["upsteam", "upsteam_diff"])
        if "CCRE" in gene_regions:
            self.assert_required_columns(ccre_data, ["CCRE", "CCRE_diff"])

        n_rows = len(gene_regions) // 2
        n_cols = len(gene_regions) - n_rows
        possible_gene_regions = ["intron", "exon", "upstream", "CCRE"]
        assert all([region in possible_gene_regions for region in gene_regions])


        ## code
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5*len(gene_regions), 4.5*len(gene_regions)))
        for j, type in enumerate(gene_regions):
            row, col = divmod(j, n_cols) # #########################
            ax = axes[row, col] #######################
            type = type.lower() if type != "CCRE" else type
            # get the variable {type}_cuttoff
            cutoff = eval(f"{type}_cutoff")
            data = pl.from_pandas(ccre_data if type == "CCRE" else gene_data)
            counts = data[type].to_numpy()
            if counts.size ==0: continue 
            counts_99_9 = np.percentile(counts, 99.9)
            counts = np.minimum(counts, counts_99_9)

            avg_diff = data[f"{type}_diff"].to_numpy()*100
            mask = (counts != 0)
            if cutoff != -1:
                mask &= (counts <= cutoff)
            counts = counts[mask]
            avg_diff = avg_diff[mask] 
            ################################ ax = axes[j]
            ax.scatter(avg_diff, counts, s=10)
            if len(counts) != 0:
                ax.set_ylim([0, counts.max()+2])
            if len(avg_diff) != 0:
                limit = -102
                ax.set_xlim([limit, -1 * limit])
            if x_label is not None:
                ax.set_xlabel(x_label, fontsize=30) # 13
            else:
                ax.set_xlabel("Average Difference (case - control)", fontsize=30) # 13
            y_label = eval(f"{type}_y_label")
            if y_label is not None:
                ax.set_ylabel(y_label, fontsize=17)
            else:
                ax.set_ylabel(f"# CpG sites of DMR ∩ {type[0].upper() + type[1:]}", fontsize=25)
            ax.tick_params(axis='x', labelsize=30)
            ax.tick_params(axis='y', labelsize=30)
            ax.yaxis.set_major_locator(MaxNLocator(steps=[1, 10], integer=True))
        if title is not None:
            fig.suptitle(title, fontsize=30, y=0.98)
        else:
            fig.suptitle(f"", fontsize=30, y=0.98)
        plt.tight_layout(rect=[0, 0, 1, 0.94], pad=3.0, w_pad=4.0, h_pad=4.0)
        plt.savefig(name, dpi=600)
    def __prepare_gene_methylation(self, position_data: InputProcessor.data_container,region_data: InputProcessor.data_container, left_distance: int = 1000, right_distance: int = 100, window_size: int = 100, hypermethylated: bool = True, gene_hypermethylated_min: int = 20, window_hypermethylated_min: int = 5, min_hypermethylated_windows: int = 5, hypomethylated: bool = True, gene_hypomethylated_max: int = -20, window_hypomethylated_max: int = -5, min_hypomethylated_windows: int = 5, position_count: int = 5, gtf_file: str= "gencode.v41.chr_patch_hapl_scaff.annotation.gtf", position_or_region: str = "region"):        
        """
            use the gtf file to get the upstream region for each gene
        """
        gtf = pl.read_csv(gtf_file, 
                          separator="\t",
                          skip_rows=5,
                          has_header=False,
                          columns=[0, 2, 3, 4, 6, 8],
                          new_columns=["chrom", "feature", "start", "end", "strand", "gene_info"]
                         )
        genes = gtf.filter(pl.col("feature") == "gene")
        genes = genes.with_columns(pl.col("gene_info").str.extract(r'gene_name "([^"]+)"').alias("gene"))
        upstream_regions = genes.with_columns([pl.when(pl.col("strand") == "+").then(pl.col("start") - left_distance).otherwise(pl.col("end")).alias("upstream_start"),
        pl.when(pl.col("strand") == "+").then(pl.col("start")).otherwise(pl.col("end") + right_distance).alias("upstream_end")])
        region_data = region_data.rename({"chromosome": "chrom"})
        regions = pl.from_pandas(region_data)
        if hypomethylated == False and hypermethylated == False:
            print("At least one of the hypo/hyper methylated parameters has to be True")
            return
        elif hypomethylated == True and hypermethylated == True:
            pass
        elif hypomethylated == False:
            regions = regions.filter((pl.col("diff") > 0))
        elif hypermethylated == False:
            regions = regions.filter((pl.col("diff") < 0))
        if position_or_region == "region":
            overlaps = (upstream_regions.join(regions, on="chrom", how="inner").filter((pl.col("start") < pl.col("end_right")) & (pl.col("end") > pl.col("start_right"))))
        else:
            overlaps = (regions.join(upstream_regions, on="chrom", how="inner").filter((pl.col("chromStart") >= pl.col("start")) & (pl.col("chromStart") <= pl.col("end"))))
            overlaps = overlaps.join(overlaps.group_by("gene").agg(pl.col("gene").count().alias("position_count")), on="gene", how="left")
            overlaps = overlaps.filter(pl.col("position_count") >= position_count)
        merged = overlaps.unique(subset="gene")[["gene", "chrom", "upstream_end", "strand", "gene_info"]]
        results = np.empty((0,(left_distance+right_distance)//window_size))
        names = deque()
        final_as_polars = pl.from_pandas(position_data)
        for gene, chr, gene_loc, strand, gene_info in merged.iter_rows():
            start = gene_loc - left_distance
            end = gene_loc + right_distance - 1
            filter = (pl.col("chrom") == chr) & (pl.col("chromStart") >= start) & (pl.col("chromStart") < end)
            loop_positions = final_as_polars.filter(filter)
            loop_positions = loop_positions.with_columns(pl.col("chromStart") - start)
            locs = loop_positions["chromStart"].to_numpy()
            diffs = (loop_positions["diff"] * 100).to_numpy()
            matrix = np.empty(left_distance+right_distance)
            matrix[:] = np.nan
            matrix[locs] = diffs
            if strand == "-":
                matrix = matrix[::-1]
            avg = np.array_split(matrix, len(matrix)//window_size)
            avg = [np.nanmean(y) if not np.isnan(y).all() else 0 for y in avg]
            hypomethylated_count = 0
            hypermethylated_count = 0
            for x in avg:
                if x <= window_hypomethylated_max:
                    hypomethylated_count += 1
                if x >= window_hypermethylated_min:
                    hypermethylated_count += 1
            if hypermethylated and hypomethylated:
                if ((hypermethylated and hypermethylated_count < min_hypermethylated_windows)
                    and (hypomethylated and hypomethylated_count < min_hypomethylated_windows)):
                    continue
            else:
                if ((hypermethylated and hypermethylated_count < min_hypermethylated_windows)
                    or (hypomethylated and hypomethylated_count < min_hypomethylated_windows)):
                    continue
            names.extend([gene])
            results = np.vstack([results, np.expand_dims(avg, axis=0)])
        return results, names
        
    def graph_upstream_gene_methylation(self, position_data: InputProcessor.data_container, region_data: InputProcessor.data_container, name: str = "upstream_methylation.png", csv_name: str = "upstream_methylation.csv", csv: Optional[str] = None, left_distance: int = 1000, right_distance: int = 1000, window_size: int = 100, hypermethylated: bool = True, gene_hypermethylated_min: int = 20, window_hypermethylated_min: int = 5, min_hypermethylated_windows: int = 5, hypomethylated: bool = True, gene_hypomethylated_max: int = -20, window_hypomethylated_max: int = -5, min_hypomethylated_windows: int = 5, position_count: int = 5, clamp_positive: int = 50, clamp_negative: int = -50, title:str = None, gtf_file: str= "gencode.v41.chr_patch_hapl_scaff.annotation.gtf", position_or_region:str = "region"):
        random.seed(42)
        np.random.seed(42)
        results, gene_names = self.__prepare_gene_methylation(position_data, region_data, left_distance, right_distance, window_size, hypermethylated, gene_hypermethylated_min, window_hypermethylated_min, min_hypermethylated_windows, hypomethylated, gene_hypomethylated_max, window_hypomethylated_max, min_hypomethylated_windows, position_count, gtf_file, position_or_region)
        print(results)
        # TODO move this into the normal function, no need for this..
        title = None
        for i, x in enumerate(results):
            for j, y in enumerate(x):
                if y > clamp_positive:
                    results[i][j] = clamp_positive
                elif y < clamp_negative:
                    results[i][j] = clamp_negative

        df = pd.DataFrame(results)
        df["index"] = list(np.array(gene_names))
        df = df.set_index("index")
        cmap=LinearSegmentedColormap.from_list('rg',["r", "w", "g"], N=256)

        print("results plots", results) ###########################

        cbar_limit = max(abs(results.min()), abs(results.max()))
        # cg = seaborn.clustermap(df.iloc[:,:df.shape[1]//2], col_cluster=False, cmap=cmap, vmin=-1*cbar_limit, center=0, vmax=cbar_limit, method='ward', figsize=(10,max(10, df.shape[0]*0.02))) #TODO add parameter

        left_bins = (left_distance)//window_size
        # cg = seaborn.clustermap(df.iloc[:,:left_bins], col_cluster=False, cmap=cmap, vmin=-1*cbar_limit, center=0, vmax=cbar_limit, method='ward', figsize=(10,max(10, df.shape[0]*-0.033 + 80)))
        fig_height = max(12, df.shape[0] * 0.1) # max(10, df.shape[0] * 0.035)
        # cg = seaborn.clustermap(df.iloc[:,:left_bins], col_cluster=False, cmap=cmap, vmin=-1*cbar_limit, center=0, vmax=cbar_limit, method='ward', figsize=(10,max(10, df.shape[0]*-0.012 + 49)))
        cg = seaborn.clustermap(df.iloc[:,:left_bins], col_cluster=False, cmap=cmap, vmin=-1*cbar_limit, center=0, vmax=cbar_limit, method='ward', figsize=(10,fig_height))
        cg.cax.set_visible(False)
        cg.ax_col_dendrogram.set_visible(False)
        ax = cg.ax_heatmap
        idxs = cg.dendrogram_row.reordered_ind
        df = df.iloc[idxs]
        seaborn.heatmap(df,ax=ax,cmap=cmap, vmin=-1*cbar_limit, center=0, vmax=cbar_limit)
        # ax.set_xticks([0] + [((left_distance+right_distance)//window_size)//2 + 0.5] + [((left_distance+right_distance)//window_size)], labels=[f'-{f"{left_distance//1000}k" if left_distance >= 1000 else left_distance}bp', 'Gene Start', f'{f"{right_distance//1000}k" if right_distance >= 1000 else right_distance}bp'])
        ###############################
        ax.set_xticks([0, (left_distance + right_distance)//window_size])
        ax.set_xticklabels([
        f'-{left_distance//1000}Kbp' if left_distance >= 1000 else f'-{left_distance}bp',
        f'{right_distance//1000}Kbp' if right_distance >= 1000 else f'{right_distance}bp'
        ])
        ax.tick_params(axis='x', labelrotation=0, labelsize=25)
        # ax.annotate("Gene Start",
        # xy=(left_bins + 0.5, -0.05),
        # xycoords=('data', 'axes fraction'),
        # ha='center', va='top',
        # fontsize=20)
        ###############################
        ax.tick_params(axis='x', labelrotation=0, labelsize=25) # no fontsize
        ax.axvline(x = left_bins + 0.5, color = 'black', linestyle = '-')
        ax.get_yaxis().set_visible(False)
        if title is not None:
            ax.set_title(title, fontsize=25)
        else:
            subtext = ""
            if hypermethylated and hypomethylated:
                subtext = "\n(Hypermethylated and Hypomethylated)"
            elif not hypermethylated and hypomethylated:
                subtext = "\n(Hypomethylated)"
            elif hypermethylated and not hypomethylated:
                subtext = "\n(Hypermethylated)"
            ax.set_title("", fontsize=25)

        plt.tight_layout()
        plt.savefig(name, dpi=600, bbox_inches='tight')
        df.to_csv(csv_name)

    def graph_upstream_UCSC(self, gene_name: str, position_data: InputProcessor.data_container, name: str="UCSC_graph.bedGraph", before_tss: int = 5000, gtf_file: str = ""):
        """
        
        Required Columns: ["chrom", "chromStart", "blockSizes_*"]
        
        """
        
        assert isinstance(position_data, pd.DataFrame), "List input not acceptable for this function."
        
        self.assert_required_columns(position_data, ["chrom", "chromStart", "blockSizes_.*"])


        gene_locations = pl.read_csv(gtf_file, 
                        skip_rows=5, 
                        separator="\t", 
                        has_header=False,
                        columns= [0,2,3,4,6,8],
                        new_columns=["chr", "gene_type", "start", "end", "strand", "gene"]
                        ) \
                        .filter((pl.col("gene_type") == "gene") & (pl.col("chr").str.contains("chr"))) \
                        .select(pl.exclude(["gene_type"])) 

        # get just the gene name and gene type from the gene info string
        gene_locations = gene_locations.with_columns(pl.col("gene").str.split(";"))

        gene_locations = gene_locations.with_columns(pl.col("gene").list.get(1).alias("gene_type")) \
            .with_columns(pl.col("gene_type").str.slice(12)) \
            .with_columns(pl.col("gene_type").str.head(-1))

        gene_locations = gene_locations.with_columns(pl.col("gene").list.get(2)) \
            .with_columns(pl.col("gene").str.slice(12)) \
            .with_columns(pl.col("gene").str.head(-1))

        gene_locations = gene_locations.with_columns(pl.col("start")-1)
        gene_locations = gene_locations.with_columns(pl.col("end")-1)

        
        chr, start, end, strand, gene, gene_type = gene_locations.filter(pl.col("gene") == gene_name)[0].to_numpy()[0]

        print(gene, start, end, strand)

        # get 5kbp upstream
        if strand == "+":
            end = start
            start -= before_tss
        elif strand == "-":
            start = end
            end += before_tss
        print(gene, start, end, strand)
            
        with open(name, "w") as file:
            file.write(f"browser position {chr}:{start}-{end}\n")
            data = position_data.filter(like="blockSizes")
            for key in data.columns:
                key_trim = key[key.find("_")+1:]
                case_ctr_str = key_trim[:key_trim.find("_")]
                color = "200,100,0" if case_ctr_str == "ctr" else "0,100,200" if case_ctr_str == "case" else "200,0,100"
                file.write(f"track type=bedGraph name=\"{gene} {key}\" visibility=full color={color} altColor=0,100,200\n")

                positions = position_data[["chrom", "chromStart", key]]
                
                positions = positions[(positions["chrom"] == chr) & (positions["chromStart"] >= start) & (positions["chromStart"] < end)]
                
                for _, row in positions.iterrows():
                    _, inner_start, percent = row
                    file.write(f"{chr}\t{inner_start}\t{inner_start+1}\t{percent}\n")

    def graph_full_gene(self, gene_name:str, position_data: InputProcessor.data_container, name="gene_methylation_graph.png", before_tss: int = 0, after_tss: Optional[int] = None, bin_size: int = 500, start_marker: bool = True, end_marker: bool = True, deviation_display: bool = True, aggregate_samples: bool = False, legend_size:int = 12, title: str = None, x_label:str = None, y_label:str=None, case_name: str = "Case", ctr_name: str = "Control", gtf_file: str = ""):        
        """
        
        Required Columns: ["chrom", "chromStart", "blockSizes_*"]
        
        """

        assert isinstance(position_data, pd.DataFrame), "List input not acceptable for this function."
        self.assert_required_columns(position_data, ["chrom", "chromStart", "blockSizes_.*"])
        
        
        gene_locations = pl.read_csv(gtf_file, 
                         skip_rows=5, 
                         separator="\t", 
                         has_header=False,
                         columns= [0,2,3,4,6,8],
                         new_columns=["chrom", "gene_type", "chromStart", "chromEnd", "strand", "gene"]
                        ) \
                        .filter((pl.col("gene_type") == "gene") & (pl.col("chrom").str.contains("chr"))) \
                        .select(pl.exclude(["gene_type"])) 
                # get just the gene name and gene type from the gene info string
        gene_locations = gene_locations.with_columns(pl.col("gene").str.split(";"))
        
        gene_locations = gene_locations.with_columns(pl.col("gene").list.get(1).alias("gene_type")) \
               .with_columns(pl.col("gene_type").str.slice(12)) \
               .with_columns(pl.col("gene_type").str.head(-1))
        
        gene_locations = gene_locations.with_columns(pl.col("gene").list.get(2)) \
               .with_columns(pl.col("gene").str.slice(12)) \
               .with_columns(pl.col("gene").str.head(-1))

        gene_locations = gene_locations.with_columns(pl.col("chromStart")-1)
        gene_locations = gene_locations.with_columns(pl.col("chromEnd")-1)
        
        gene = gene_locations.filter(pl.col("gene") == gene_name)[0]
        print(gene)
        chrom = gene[0]["chrom"].item()
        start = gene[0]["chromStart"].item()
        end = gene[0]["chromEnd"].item()
        strand = gene[0]["strand"].item()

        if strand == "-" and after_tss is not None:
            start, end = end, start
            before_tss, after_tss = after_tss, before_tss

        print(position_data, before_tss, after_tss, start, end)
        if after_tss is None:
            df = position_data[(position_data["chromStart"] >= start - before_tss) & (position_data["chromStart"] <= end) & (position_data["chrom"] == chrom)]
        else:
            df = position_data[(position_data["chromStart"] >= start - before_tss) & (position_data["chromStart"] <= start + after_tss) & (position_data["chrom"] == chrom)]

        if strand == "-" and after_tss is None:
            start, end = end, start
            
        # Create bin labels
        bin = []
        print(df)
        current = df["chromStart"].iloc[0]
        for _, row in df.iterrows():
            while current + bin_size < row["chromStart"]:
                current += bin_size

            bin.extend([current])
            
            # bin = (df['chromStart'] - 1) // bin_size * bin_size + 1
            # print(bin)
    
            # bin.iloc[0] = df.iloc[0]["chromStart"]
            # bin.iloc[-1] = df.iloc[-1]["chromStart"]
        # print(bin)
        bin[-1] =  df.iloc[-1]["chromStart"]
        df['position_bin'] = bin
            
        print(df)

        if aggregate_samples:
            # generate average of blockSizes_ctr* and blockSizes_case*
            df["avg_ctr"] = df.filter(like="blockSizes_ctr").mean(axis=1)
            df["avg_case"] = df.filter(like="blockSizes_case").mean(axis=1)
            print(df)

        cols = df.filter(like="blockSizes").columns if not aggregate_samples else ["avg_ctr", "avg_case"]

        # Aggregate by bins
        result = df.groupby('position_bin', as_index=False).agg({ x : [np.nanmean, np.nanstd] for x in cols})

        # print(result["blockSizes_case_0"])
        print(result)
        plt.clf()
        
        fig, ax = plt.subplots(figsize=(15, 4)) # figsize=(15, 4)
        
        for x in cols:

            # label = "WT" if "ctr" in x else "KO"
            label = ctr_name if "ctr" in x else case_name
            color = "orange" if "ctr" in x else "blue"
            # 2 4
            ax.plot(result["position_bin"], result[x]["nanmean"], 'o--', linewidth=3, ms=6, label=label, color = color) # ax.plot(result["position_bin"], result[x]["nanmean"]/100, 'o--', linewidth=3, ms=6, label=label, color = color)
            if deviation_display:
                ax.fill_between(result["position_bin"], 
                                (result[x]["nanmean"] - result[x]["nanstd"])/100, # /100 
                                (result[x]["nanmean"] + result[x]["nanstd"])/100, # /100
                                color=color, alpha=0.2)
        print(start)
        # vertical lines for gene start/end
        # TODO make this optional
        if start_marker:
            ax.axvline(start, color="green", linestyle=":", label="Gene Start")#, linewidth=5)

        if end_marker and (after_tss == None or (strand == "+" and df["chromStart"].iloc[-1] >= end) or (strand == "-" and df["chromStart"].iloc[0] <= end)):
            ax.axvline(end, color="red", linestyle=":", label="Gene End")
            
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), prop={'size':legend_size})

        if x_label is not None:
            ax.set_xlabel(x_label)
        else:
            ax.set_xlabel(f"Position (Chromosome {chrom[3:]})")
        
        if y_label is not None:
            ax.set_ylabel(y_label)
        else:
            ax.set_ylabel("Methylation Level")
        
        ax.set_ylim(0, 1)

        if title is not None:
            ax.set_title(title)
        else:
            plt.title("Methylation Levels in Gene " + gene_name)
        
        plt.tight_layout()

        loc, label = plt.xticks()
        loc = loc[1:-1]
        label = [format(int(x), ",") for x in loc]
        plt.xticks(loc, label)
        
        plt.savefig(name, dpi=600, bbox_inches='tight')
       
    def __sliding_window_avg(self, df, start_col, window_size, step_size, sample_start_ind):
        new_rows = []
        min_pos, max_pos = df[start_col].min(), df[start_col].max()
        for start in range(int(min_pos), int(max_pos) - int(window_size) + 1, int(step_size)):
            end = start + window_size
            window_data = df[(df[start_col] >= start) & (df[start_col] < end)]
            if not window_data.empty:
                # print("window_data.iloc[:, sample_start_ind:]:  ", window_data.iloc[:, sample_start_ind:])
                # avg_values = window_data.iloc[:, sample_start_ind:].mean().to_dict() 
                avg_values = window_data.filter(like="blockSizes_").astype(float).mean().to_dict()
                avg_values[start_col] = start + window_size // 2 
                new_rows.append(avg_values)
        return pd.DataFrame(new_rows)
    
    def plot_methylation_curve(self, region_data: InputProcessor.data_container, position_data: InputProcessor.data_container,  name:str, repeat_regions_df: str=None, enhancer_promoter_df: str = None, repeat_regions_columns:list[int] = [5,6,7,11], enhancer_promoter_columns:list[int] = [0,1,2,12,13], window_size:int = 50, step_size:int = 25, chr_filter:str = None, start_filter:int = None, end_filter:int = None, sample_start_ind:int = 3):
        """
        Generate smoothed methylation plots using sliding window averaging.
        - Left subplot: Individual sample curves.
        - Right subplot: Mean ± STD curves for two groups.
        region_data: BED format DataFrame ['chromosome', 'start', 'end'].
        position_data: DataFrame with ['chromosome', 'start'] + sample columns.
        """
        self.assert_required_columns(region_data, ["chrom", "start", "end"])
        self.assert_required_columns(position_data, ["chrom", "chromStart", "blockSizes_case_.*", "blockSizes_ctr_.*"])
        repeat_colors = {
            "SINE": "blue",
            "LINE": "teal",
            "LTR": "orchid",
            "Simple_repeat": "green",
            "Satellite": "pink" ,
            "Low_complexity": "brown",
            "DNA": "gray"
        }
        enhancer_colors = {
            "enhP": "red",
            "enhD": "orange",
            "prom": "gold"
        }
        annotation_height = 0.02
        add_y = 0.03
        min_std_area = 0.075;
        sep_area_dict = []
        repeat_regions_df = pl.read_csv(repeat_regions_df, separator='\t', has_header=False, columns=repeat_regions_columns, new_columns=['chromosome', 'start', 'end', 'labels'], infer_schema_length=10000).to_pandas()
        enhancer_promoter_df = pl.read_csv(enhancer_promoter_df, separator='\t', has_header=False, columns=enhancer_promoter_columns, new_columns=['chromosome', 'start', 'end', 'cat', 'id'], infer_schema_length=10000).to_pandas()
        if 'cat' in enhancer_promoter_df.columns and 'id' in enhancer_promoter_df.columns: enhancer_promoter_df['labels'] =  enhancer_promoter_df['cat'] + '_' + enhancer_promoter_df['id']
        group1_samples = position_data.filter(like="blockSizes_case").columns.tolist()
        group2_samples = position_data.filter(like="blockSizes_ctr").columns.tolist()
        if chr_filter != None: 
            region_data = region_data[region_data["chrom"] == chr_filter]
            position_data = position_data[position_data["chrom"] == chr_filter]
            if start_filter != None:
                if end_filter != None:
                    region_data = region_data[(region_data["start"] >= start_filter) & (region_data["end"] <= end_filter)]
                else:
                    region_data = region_data[(region_data["start"] >= start_filter)]
        print(region_data)
        for _, region in region_data.iterrows():
            chrom, start, end = region['chrom'], region['start'], region['end']
            base_fn = "{}-{}-{}".format(chrom, start, end)
            # Extract methylation data for the region
            region_data = position_data[
                (position_data['chrom'] == chrom) & 
                (position_data['chromStart'] >= start) & 
                (position_data['chromStart'] <= end)
            ].copy()
            if region_data.empty:
                print(f"No data for region: {chrom}:{start}-{end}")
                continue
            windowed_data = self.__sliding_window_avg(region_data, 'chromStart', window_size, step_size, sample_start_ind=3)
            if windowed_data.empty:
                print(f"No data for window_data: {chrom}:{start}-{end}")
                print( region_data )
                continue
            # Prepare figure with two subplots
            fig, axes = plt.subplots(2, 1, figsize=(20, 12)) ###################################################### 15, 5
            # === Left Subplot: Individual Sample Curves ===
            for sample in group1_samples:
                if sample in windowed_data.columns:
                    axes[0].plot(windowed_data['chromStart'], windowed_data[sample], label=sample, color='blue', alpha=0.5)
            for sample in group2_samples:
                if sample in windowed_data.columns:
                    axes[0].plot(windowed_data['chromStart'], windowed_data[sample], label=sample, color='red', alpha=0.5)
            axes[0].set_xlabel("Genomic Coordinate", fontsize=20) # fontsize=12 ; labelsize=12
            axes[0].set_ylabel("Methylation Level", fontsize=20)
            axes[0].tick_params(axis='both', labelsize=20)
            # === Right Subplot: Mean ± STD Curves for Groups ===
            print(windowed_data.columns)
            windowed_data['Group1_Mean'] = windowed_data[group1_samples].mean(axis=1)
            windowed_data['Group1_STD'] = windowed_data[group1_samples].std(axis=1)
            windowed_data['Group2_Mean'] = windowed_data[group2_samples].mean(axis=1)
            windowed_data['Group2_STD'] = windowed_data[group2_samples].std(axis=1)
            x_range = [windowed_data['chromStart'].iloc[0], windowed_data['chromStart'].iloc[-1] ]
            axes[1].plot(windowed_data['chromStart'], windowed_data['Group1_Mean'], label="Case", color='blue', linestyle='-')
            axes[1].fill_between(windowed_data['chromStart'], windowed_data['Group1_Mean'] - windowed_data['Group1_STD'], 
                                 windowed_data['Group1_Mean'] + windowed_data['Group1_STD'], color='blue', alpha=0.2)
            axes[1].plot(windowed_data['chromStart'], windowed_data['Group2_Mean'], label="Control", color='red', linestyle='-')
            axes[1].fill_between(windowed_data['chromStart'], windowed_data['Group2_Mean'] - windowed_data['Group2_STD'], 
                                 windowed_data['Group2_Mean'] + windowed_data['Group2_STD'], color='red', alpha=0.2)
            axes[1].set_xlabel("Genomic Coordinate", fontsize=20)
            axes[1].set_ylabel("Average Methylation ± STD", fontsize=20)
            axes[1].tick_params(axis='both', labelsize=20)
            def_loc = 'center right'
            low_part = 0.2;
            if windowed_data['Group1_Mean'].iloc[-1]-windowed_data['Group1_STD'].iloc[-1]>low_part and  windowed_data['Group2_Mean'].iloc[-1]-windowed_data['Group2_STD'].iloc[-1]>low_part:
               def_loc = 'lower right'
            elif windowed_data['Group1_Mean'].iloc[0]-windowed_data['Group1_STD'].iloc[0]>low_part and  windowed_data['Group2_Mean'].iloc[0]-windowed_data['Group2_STD'].iloc[0]>low_part:
               def_loc = 'lower left'
            elif windowed_data['Group1_Mean'].iloc[0]+windowed_data['Group1_STD'].iloc[0]<100-low_part and  windowed_data['Group2_Mean'].iloc[0]+windowed_data['Group2_STD'].iloc[0]<100-low_part:
               def_loc = 'upper left'
            elif windowed_data['Group1_Mean'].iloc[-1]+windowed_data['Group1_STD'].iloc[-1]<100-low_part and  windowed_data['Group2_Mean'].iloc[-1]+windowed_data['Group2_STD'].iloc[-1]<100-low_part:
               def_loc = 'upper right'
            axes[1].legend(fontsize=15, loc=def_loc, ncol=2) # fontsize=10
            xmin, xmax = max(axes[0].get_xlim()[0], axes[1].get_xlim()[0]), min(axes[0].get_xlim()[1], axes[1].get_xlim()[1])
            if not (repeat_regions_df is None):
               repeat_in_region = repeat_regions_df[
                      (repeat_regions_df['chromosome'] == chrom) &
                      (repeat_regions_df['start'] <= end) &
                      (repeat_regions_df['end'] >= start)
               ]
               print(base_fn, repeat_in_region.shape)
               for ax_i in range(2):
                  ymin, ymax = axes[ax_i].get_ylim()
                  annotation_y = ymax + add_y  
                  for _, repeat in repeat_in_region.iterrows():
                      repeat_type = repeat.get('labels', "Other")  
                      color = repeat_colors.get(repeat_type, "gray")  
                      plot_start = xmin if repeat['start']<xmin else repeat['start']
                      plot_end = xmax if repeat['end']>xmax else repeat['end']
                      axes[ax_i].add_patch(plt.Rectangle((plot_start, annotation_y), plot_end - plot_start, annotation_height, 
                                                      color=color, alpha=0.7, label=repeat['labels']))
                      axes[ax_i].text( (plot_end+plot_start) / 2, 
                           annotation_y + annotation_height,
                           repeat['labels'],
                           ha='center', va='center', fontsize=12, color='black', fontweight='bold'
                       )# fontsize=7
                  if repeat_in_region.shape[0]>0:
                      axes[ax_i].set_ylim(0, annotation_y + add_y)
            if not (enhancer_promoter_df is None):
               enhancer_in_region = enhancer_promoter_df[
                      (enhancer_promoter_df['chromosome'] == chrom) &
                      (enhancer_promoter_df['start'] <= end) &
                      (enhancer_promoter_df['end'] >= start)
               ]
               print( base_fn, enhancer_in_region.shape)
               for ax_i in range(2):
                  ymin, ymax = axes[ax_i].get_ylim()
                  annotation_y = ymax + add_y 
                  for _, enhancer in enhancer_in_region.iterrows():
                      enhancer_type = enhancer.get('cat', "Other")  
                      color = enhancer_colors.get(enhancer_type, "black") 
                      plot_start = xmin if enhancer['start']<xmin else enhancer['start']
                      plot_end = xmax if enhancer['end']>xmax else enhancer['end']
                      axes[ax_i].add_patch(plt.Rectangle((plot_start, annotation_y), plot_end - plot_start, annotation_height, 
                                                      color=color, alpha=0.7, label=enhancer['labels']))
                      axes[ax_i].text((plot_end+plot_start) / 2,  
                           annotation_y + annotation_height,  
                           enhancer['labels'],
                           ha='center', va='center', fontsize=12, color='black', fontweight='bold'
                       ) # fontsize=7
                  if enhancer_in_region.shape[0]>0:
                     axes[ax_i].set_ylim(0, annotation_y + add_y)
            plt.tight_layout()
            plt.savefig(name+'/long_region_plot_'+base_fn+'.png', dpi=600)
            plt.close('all')
            t_n_ear = 0;
            for _i_sa in range(len(windowed_data['Group1_Mean'])):
               t_n_ear += (windowed_data['Group1_Mean'].iloc[_i_sa] - windowed_data['Group2_Mean'].iloc[_i_sa])/( (windowed_data['Group1_STD'].iloc[_i_sa] if windowed_data['Group1_STD'].iloc[_i_sa]>min_std_area else min_std_area)+ (windowed_data['Group2_STD'].iloc[_i_sa] if windowed_data['Group2_STD'].iloc[_i_sa]>min_std_area else min_std_area) )
            t_n_ear = t_n_ear / (len(windowed_data['Group1_Mean']) if len(windowed_data['Group1_Mean']) >=10 else 10)
            sep_area_dict.append([abs(t_n_ear), t_n_ear, base_fn])
        return sorted(sep_area_dict)
    def match_region_annotation(self, regions_df: InputProcessor.data_container, bed_file: InputProcessor.data_container, name: str = "match_region_annotation", annotation_or_region: str = "region", show_counts: bool = False):
        """
        Reads region file and annotation file, finds matches, and counts occurrences.
        """
        print(bed_file)
        if bed_file == "CpG_gencodev42ccrenb_repeat_epic1v2hm450.bed": bed_file = Path(__file__).resolve().parent.parent / bed_file
        annotation_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'id', 'annotation'])
        total_counts_1 = defaultdict(int)
        total_counts_2 = defaultdict(int)
        total_counts_3 = defaultdict(int)
        total_counts_4 = defaultdict(int)
        total_counts_p4 = defaultdict(int)
        total_gene = {}
        special_types = ['enhD', 'enhP', 'prom', 'K4m3', 'enh']
        add_to_name = ""
        # print(regions_df)
        chroms = regions_df['chrom'].unique()
        # print(chroms)
        for chrom in list(chroms):
            print(chrom)
            regions_df_chr = regions_df[regions_df["chrom"] == chrom]
            annotation_df_chr = annotation_df[annotation_df['chr'] == chrom]
            for _, region in regions_df_chr.iterrows():
                matched_annotations = annotation_df_chr[
                (annotation_df_chr['start'] <= region['end']) &
                (annotation_df_chr['end'] >= region['start'])
                ]
                this_counts_1 = defaultdict(int)
                this_counts_2 = defaultdict(int)
                en_list = defaultdict(int)
                this_counts_3 = defaultdict(int)
                this_counts_4 = defaultdict(int)
                gend_counts = defaultdict(int)
                if annotation_or_region == "annotation":
                    add_to_name = "annotation_based"
                    all_categories = []
                    all_encode_types = []
                    all_repeat = []
                    all_cpg_epic = []
                    all_gene_list = {}
                    for _, annotation in matched_annotations.iterrows():
                        try:
                            categories, encode_types, repeat, cpg_epic, gene_list = Analysis.parse_annotation(annotation['annotation'])
                        except AttributeError:
                            continue
                        all_categories.extend(categories)
                        all_encode_types.extend(encode_types)
                        all_repeat.extend(repeat)
                        # print(all_repeat)
                        all_cpg_epic.extend(cpg_epic)
                        for g in gene_list:
                            if g not in all_gene_list:
                                all_gene_list[g] = set(gene_list[g])
                            else:
                            #for k in gene_list[g]:
                                all_gene_list[g] |= set(gene_list[g])
                                # all_gene_list[g][k] = all_gene_list[g].get(k, 0) + gene_list[g][k]
                    for cat in set(all_categories):
                        if cat in special_types:
                            this_counts_1[cat] += 1
                            key = 'Gene_Body' if 'Gene_' in cat else cat
                            this_counts_2[key] += 1
                        else:
                            count_cat = all_categories.count(cat)
                            this_counts_1[cat] += 1 ############################################# count_cat
                            key = 'Gene_Body' if 'Gene_' in cat else cat
                            this_counts_2[key] += 1 ############################################### count_cat
                    if len(all_repeat) == 0:
                        this_counts_3['No_repeat'] += 1
                    else:
                        for rep_c in set(all_repeat):
                            this_counts_3[rep_c] += 1 ################################################all_repeat.count(rep_c)
                    if len(all_cpg_epic) == 0:
                        this_counts_4['non-EPIC'] += 1
                    else:
                        this_counts_4['EPIC'] += 1
                    for encode_type in set(all_encode_types):
                        this_counts_1[encode_type] += 1
                        this_counts_2[encode_type] += 1
                        en_list[encode_type] += 1
                    for _ge in all_gene_list:
                        for _de in all_gene_list[_ge]:
                            key = _ge + ':' + _de
                            gend_counts[key] += 1
                        #if key not in gend_counts:
                        #gend_counts[key] = all_gene_list[_ge][_de]
                        # else:
                        #    gend_counts[key] += all_gene_list[_ge][_de]
                else:
                    add_to_name = "region_based"
                    for _, annotation in matched_annotations.iterrows():
                        try:
                            categories, encode_types, repeat, cpg_epic, gene_list = Analysis.parse_annotation(annotation['annotation'])
                        except AttributeError:
                            continue
                        for _ge in gene_list:
                            for _de in gene_list[_ge]:
                                key = _ge + ':' + _de
                                gend_counts[key] = gend_counts.get(key, 0) + 1
                        for cat in categories:
                            this_counts_1[cat] += 1
                            this_counts_2['Gene_Body' if 'Gene_' in cat else cat] += 1
                            # print("s1 - 1: ", this_counts_1)
                            # print("s1 - 2: ", this_counts_2)
                        if len(repeat) == 0:
                            this_counts_3['No_repeat'] += 1
                        else:
                            for rep_c in repeat:
                                this_counts_3[rep_c] += 1
                        if len(cpg_epic) == 0:
                            this_counts_4['non-EPIC'] += 1
                        else:
                            this_counts_4['EPIC'] += 1
                        for encode_type in encode_types:
                            this_counts_1[encode_type] += 1
                            this_counts_2[encode_type] += 1
                            en_list[encode_type] += 1
                # print("s2 - 1: ", this_counts_1)
                # print("s2 - 2: ", this_counts_2)
                de_gede = []
                for _gede in gend_counts:
                    if _gede.endswith('_nb') and _gede[:-3] in gend_counts:
                        gend_counts[_gede[:-3]] += gend_counts[_gede]
                        de_gede.append(_gede)
                for _gede in de_gede:
                    del gend_counts[_gede]
                total_gene = Counter(total_gene) + Counter(dict(gend_counts))
                if 'Gene_Intron' in this_counts_1 or 'Gene_Exon' in this_counts_1 or len(en_list) > 0:
                    this_counts_1.pop('IG', None)
                    this_counts_2.pop('IG', None)
            # print("s3 - 1: ", this_counts_1)
            # print("s3 - 2: ", this_counts_2)
                if any(et in en_list for et in ['enhD', 'enhP', 'prom', 'K4m3', 'enh']):
                    # this_counts_1.pop('Gene_Intron', None)
                    this_counts_2['Gene_Body'] -= this_counts_1['Gene_Intron']
                    this_counts_1.pop('Gene_Intron', None)
            # print("s4 - 1: ", this_counts_1)
            # print("s4 - 2: ", this_counts_2)
                if annotation_or_region == "region" and 'No_repeat' in this_counts_3 and this_counts_3['No_repeat'] < matched_annotations.shape[0] - this_counts_3['No_repeat']:
                    this_counts_3.pop('No_repeat', None)
                en_dis = []
                en_cod_type = list(this_counts_1.keys())
                for _c in en_cod_type:
                    if '_' in _c and _c.split('_')[0] in special_types and _c.split('_')[0] in en_cod_type:
                        en_dis.append(_c)
                for _c in en_dis:
                    this_counts_1.pop(_c, None)
                    this_counts_2.pop(_c, None)
            # print("s5 - 1: ", this_counts_1)
            # print("s5 - 2: ", this_counts_2)
                for _c in this_counts_1:
                    total_counts_1[_c] += this_counts_1[_c]
                for _c in this_counts_2:
                    total_counts_2[_c] += this_counts_2[_c]
                for _c in this_counts_3:
                    total_counts_3[_c] += this_counts_3[_c]
                for _c in this_counts_4:
                    total_counts_4[_c] += this_counts_4[_c]
                for _c in this_counts_4:
                    total_counts_p4[_c] += this_counts_4[_c]
        self.__annotation_chart(total_counts_1, "Fraction of Occurrences: Gene Intron, Gene Exon, IG, and ENCODE Types", add_to_name + "_" + name+'_inex2.png')
        self.__annotation_chart(total_counts_2, "Fraction of Occurrences: Gene Body, IG, and ENCODE Types", add_to_name + "_" + name+'_body2.png')
        self.__annotation_chart(total_counts_3, "", add_to_name + "_" + name+'_repeat2.png')
        self.__annotation_chart(total_counts_4, "", add_to_name + "_" + name+'_epic2.png')
        self.__annotation_chart(total_counts_p4, "", add_to_name + "_" + name+'_epicP2.png')
        return [pd.DataFrame.from_dict(list(dict(total_counts_1).items())), pd.DataFrame.from_dict(list(dict(total_counts_2).items())), pd.DataFrame.from_dict(list(dict(total_counts_3).items())), pd.DataFrame.from_dict(list(dict(total_counts_4).items())), pd.DataFrame.from_dict(list(dict(total_counts_p4).items())), pd.DataFrame.from_dict(list(dict(total_gene).items()))]
    
    def match_position_annotation(self, regions_df: InputProcessor.data_container, bed_file: InputProcessor.data_container, name:str="match_position_annotation", show_counts = False):
        """
        Reads position file and annotation file, finds matches, and counts occurrences.
        """
        annotation_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'id', 'annotation'])
        total_counts_1 = defaultdict(int)
        total_counts_2 = defaultdict(int)
        total_counts_3 = defaultdict(int)
        total_counts_4 = defaultdict(int)
        total_counts_p4 = defaultdict(int)
        total_gene = {}
        # regions_df = pd.read_csv(regions_df)
        merged_df = pd.merge(annotation_df, regions_df, left_on=["chr", "start"], right_on=["chrom", "chromStart"], how="inner")
        this_counts_1 = defaultdict(int)
        this_counts_2 = defaultdict(int)
        en_list = defaultdict(int)
        this_counts_3 = defaultdict(int)
        this_counts_4 = defaultdict(int)
        gend_counts = defaultdict(int)
        merged_df = merged_df.dropna(subset=["annotation"])
        for _, annotation in merged_df.iterrows():
            if annotation['annotation'] == "nan" : continue
            # print(annotation['annotation'])
            categories, encode_types, repeat, cpg_epic, gene_list = Analysis.parse_annotation(annotation['annotation'])
            for _ge in gene_list:
               for _de in gene_list[ _ge ]:
                  if (_ge + ':'+_de) not in gend_counts: gend_counts[  _ge + ':'+_de ] = 1;
                  else: gend_counts[  _ge + ':'+_de ] += 1
            for cat in categories:
                this_counts_1[cat] += 1
                this_counts_2['Gene_Body' if 'Gene_' in cat else cat] += 1
            if len(repeat)==0: this_counts_3['No_repeat'] += 1;
            else:
                for rep_c in repeat:
                    this_counts_3[rep_c] += 1
            if len(cpg_epic)==0:
                this_counts_4['non-EPIC'] += 1;
            else:
                this_counts_4['EPIC'] += 1;
            for encode_type in encode_types:
                this_counts_1[encode_type] += 1
                this_counts_2[encode_type] += 1
                en_list[encode_type] += 1
        de_gede = []
        for _gede in gend_counts:
           if _gede[-3:]=='_nb' and _gede[:-3] in gend_counts:
              gend_counts[ _gede[:-3] ] += gend_counts[ _gede ]
              de_gede.append( _gede)
        for _gede in de_gede:
           del gend_counts[ _gede ]
        total_gene = Counter(total_gene) + Counter(dict(gend_counts))
        # if 'Gene_Intron' in this_counts_1 or 'Gene_Exon' in this_counts_1 or len(en_list)>0:
        #   if 'IG' in this_counts_1: del this_counts_1['IG']
        #   if 'IG' in this_counts_2: del this_counts_2['IG']
        #if 'enhD' in en_list or 'enhP' in en_list or 'prom' in en_list or 'K4m3' in en_list or 'enh' in en_list:
        #   if 'Gene_Intron' in this_counts_1: del this_counts_1['Gene_Intron']
        # if 'No_repeat' in this_counts_3 and this_counts_3['No_repeat']< matched_annotations.shape[0] - this_counts_3['No_repeat']:
        #    del this_counts_3['No_repeat']
        # if 'Gene_Exon' in this_counts_1:
        #    print("Exon", region)
        en_dis = []
        print(this_counts_1, this_counts_2, this_counts_3, this_counts_4) #, gend_counts)
        en_cod_type = list(this_counts_1.keys())
        for _c in en_cod_type:
           if '_' in _c and _c.split('_')[0] in ['enhD', 'enhP', 'prom', 'K4m3', 'enh'] and _c.split('_')[0] in en_cod_type:
              en_dis.append( _c )
        for _c in en_dis:
           if _c in this_counts_1: del this_counts_1[_c]
           if _c in this_counts_2: del this_counts_2[_c]
        for _c in this_counts_1:
           total_counts_1[ _c ] = 1 + total_counts_1[ _c ]
        for _c in this_counts_2:
           total_counts_2[ _c ] = 1 + total_counts_2[ _c ]
        for _c in this_counts_3:
           total_counts_3[ _c ] = 1 + total_counts_3[ _c ]
        for _c in this_counts_4:
           total_counts_4[ _c ] = 1 + total_counts_4[ _c ]
        for _c in this_counts_4:
           total_counts_p4[ _c ] = this_counts_4[ _c ] + total_counts_p4[ _c ]
        self.__annotation_chart(total_counts_1, "Fraction of Occurrences: Gene Intron, Gene Exon, IG, and ENCODE Types", name+'_inex2.png', nb = show_counts)
        self.__annotation_chart(total_counts_2, "Fraction of Occurrences: Gene Body, IG, and ENCODE Types", name+'_body2.png', nb = show_counts)
        self.__annotation_chart(total_counts_3, "", name+'_repeat2.png', nb = show_counts)
        self.__annotation_chart(total_counts_4, "", name+'_epic2.png', nb = show_counts)
        self.__annotation_chart(total_counts_p4, "", name+'_epicP2.png', nb = show_counts)
        return [pd.DataFrame.from_dict(list(dict(total_counts_1).items())), pd.DataFrame.from_dict(list(dict(total_counts_2).items())), pd.DataFrame.from_dict(list(dict(total_counts_3).items())), pd.DataFrame.from_dict(list(dict(total_counts_4).items())), pd.DataFrame.from_dict(list(dict(total_counts_p4).items())), pd.DataFrame.from_dict(list(dict(total_gene).items()))]
    def __annotation_chart(self, data, title, fig_name, keep_small=False, nb = False):
        """
        Plot and save a pie chart.
        """
        print(data)
        labels = sorted(list(data.keys()))
        sizes = [data[_lb_] for _lb_ in labels]
        total = sum(sizes)
        filtered_labels = []
        filtered_sizes = []
        var = 1 if not keep_small else 0
        for label, size in zip(labels, sizes):
            percent = 100 * size / total
            if percent >= var:
                filtered_labels.append(label)
                filtered_sizes.append(size)
        def make_autopct(nb):
            def autopct_func(pct):
                absolute = int(round(pct * total / 100.0))
                if nb:
                    return f'{pct:.1f}%\n({absolute})' if pct >= 1 else ''
                else:
                    return f'{pct:.1f}%' if pct >= 1 else ''
            return autopct_func
        colors = plt.get_cmap('Blues')(np.linspace(0.2, 0.7, len(filtered_sizes)))
        plt.figure(figsize=(7, 7))
        plt.pie(
            filtered_sizes,
            labels=filtered_labels,
            autopct=make_autopct(nb),
            startangle=140,
            colors=colors,
            textprops={'fontsize': 18},
            labeldistance=1.1
        )
        plt.title(title, fontsize=16)
        plt.savefig(fig_name, dpi=600, bbox_inches='tight')
        plt.close('all')
    def __match_to_gene(self, positions_df, bed_file, name: str = "match_region_annotation"):
        """
            Reads region file and annotation file, finds matches, and returns gene annotation matrix with average diff.
        """
        annotation_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'id', 'annotation'])
        # regions_df = pd.read_csv(regions_df)
        merged_df = pd.merge(annotation_df, positions_df, left_on=["chr", "start"], right_on=["chrom", "chromStart"], how="inner").dropna(subset=["annotation"])
        print(merged_df)
        gend_counts = defaultdict(int)
        gend_counts_diff = defaultdict(list)
        for _, annotation in merged_df.iterrows():
            categories, encode_types, repeat, cpg_epic, gene_list = Analysis.parse_annotation(annotation['annotation'])
            for _ge in gene_list:
                for _de in gene_list[_ge]:
                    key = f"{_ge}:{_de}"
                    gend_counts[key] += 1
                    gend_counts_diff[key + "_diff"].append(annotation["diff"])
        cols = {'CDS', 'CTCF', 'CTCF_nb', 'K4m3', 'K4m3_nb', 'enhD', 'enhD_nb', 'enhP', 'enhP_nb', 'exon', 'intron', 'prom', 'prom_nb'}
        d = {}
        for key, val in gend_counts.items():
            k = key.split(":")
            if len(k) != 2: continue
            col, gene = (k[1], k[0]) if k[0] not in cols else (k[0], k[1])
            if col not in cols: continue
            if gene not in d:
                d[gene] = {}
            d[gene][col] = d[gene].get(col, 0) + val
            d[gene][col + "_diff"] = sum(gend_counts_diff[key + "_diff"]) / len(gend_counts_diff[key + "_diff"])
        return pd.DataFrame.from_dict(d, orient='index').fillna(0)
