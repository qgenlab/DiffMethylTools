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


class Plots():
    def __init__(self):
        pass

    def assert_required_columns(self, df, required_columns):
        Analysis().assert_required_columns(df, required_columns)
    
    def assert_one_of_column_pairs(self, df, column_pairs):
        Analysis().assert_one_of_column_pairs(df, column_pairs)

    def volcano_plot(self, data: InputProcessor.data_container, name : str, threshold : Optional[float] = 0.05, line: Optional[float] = 15, x_range: tuple[int, int] = (-100, 100), y_max: int = None, title: str = None, x_label: str = None, y_label: str = None):
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
        self.assert_one_of_column_pairs(data, [("start",), ("chromStart")])
        if "start" in data.columns:
            start_name = "start"
        else:
            start_name = "chromStart"

        fig, ax = plt.subplots(1, 1, figsize=(15, 6))
        df = data
        df["neg_log10_fdr"] = -1 * np.log10(df["q-value"])
        
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
            ax.scatter(group['plot_position'], group['neg_log10_fdr'], color=colors[num % len(colors)], s=10)
            x_labels.append(grouped_name)
            x_labels_pos.append((group['plot_position'].iloc[-1] + group['plot_position'].iloc[0]) / 2)
            current_position = group['plot_position'].iloc[-1] + 1

        if threshold is not None:
            ax.axhline(y=-np.log10(threshold), color='maroon', linestyle='--')

        ax.set_xticks(x_labels_pos)
        ax.set_xticklabels(x_labels)

        if title is not None:
            ax.set_title(title, fontsize=20)
        else:
            ax.set_title('Manhattan Plot', fontsize=20)
        if x_label is not None:
            ax.set_xlabel(x_label, fontsize=15)
        else:
            ax.set_xlabel('Chromosome', fontsize=15)
        if y_label is not None:
            ax.set_ylabel(y_label, fontsize=15)
        else:
            ax.set_ylabel('-log10(p-value)', fontsize=15)

        ax.tick_params(axis='x', labelsize=11)
        ax.tick_params(axis='y', labelsize=12)
        plt.tight_layout()
        plt.savefig(name, dpi=300)

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
                ax.set_xlabel(x_label, fontsize=15)
            else:    
                ax.set_xlabel("Read Coverage", fontsize=15)
            if y_label is not None:
                ax.set_ylabel(y_label, fontsize=15)
            else:
                ax.set_ylabel("Count", fontsize=15)
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

    def graph_gene_regions(self, gene_data: InputProcessor.data_container, ccre_data: InputProcessor.data_container, name: str, gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], intron_cutoff: int = -1, exon_cutoff: int = -1, upstream_cutoff: int = -1, CCRE_cutoff: int = -1, title: str = None, x_label: str = None, intron_y_label: str = None, exon_y_label: str = None, upstream_y_label: str = None, CCRE_y_label: str = None):
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
        if "upstream" in gene_regions:
            self.assert_required_columns(gene_data, ["upstream", "upstream_diff"])
        if "CCRE" in gene_regions:
            self.assert_required_columns(ccre_data, ["CCRE", "CCRE_diff"])


        possible_gene_regions = ["intron", "exon", "upstream", "CCRE"]
        assert all([region in possible_gene_regions for region in gene_regions])

        ## code
        fig, axes = plt.subplots(1,len(gene_regions), figsize=(4.5*len(gene_regions), 7))
        for j, type in enumerate(gene_regions):
            type = type.lower() if type != "CCRE" else type
            # get the variable {type}_cuttoff
            cutoff = eval(f"{type}_cutoff")
            data = pl.from_pandas(ccre_data if type == "CCRE" else gene_data)
            counts = data[type].to_numpy()
            avg_diff = data[f"{type}_diff"].to_numpy()
            mask = (counts != 0)
            if cutoff != -1:
                mask &= (counts <= cutoff)
            counts = counts[mask]
            avg_diff = avg_diff[mask] 
            ax = axes[j]
            ax.scatter(avg_diff, counts, s=10)
            if len(counts) != 0:
                ax.set_ylim([0, counts.max()+2])
            if len(avg_diff) != 0:
                limit = min(avg_diff.min(), -1 * avg_diff.max()) - 2
                ax.set_xlim([limit, -1 * limit])
            if x_label is not None:
                ax.set_xlabel(x_label, fontsize=13)
            else:
                ax.set_xlabel("Average Difference (case - control)", fontsize=13)
            y_label = eval(f"{type}_y_label")
            if y_label is not None:
                ax.set_ylabel(y_label, fontsize=17)
            else:
                ax.set_ylabel(f"# DMP in {type[0].upper() + type[1:]}", fontsize=17)
            ax.tick_params(axis='x', labelsize=12)
            ax.tick_params(axis='y', labelsize=12)
            ax.yaxis.set_major_locator(MaxNLocator(steps=[1, 10], integer=True))
        if title is not None:
            fig.suptitle(title, fontsize=20)
        else:
            fig.suptitle(f"Gene Locations vs. Average Methylation Difference", fontsize=20)
        plt.tight_layout()
        plt.savefig(name, dpi=300)

    def __prepare_gene_methylation(self, position_data: InputProcessor.data_container, gene_data: InputProcessor.data_container, gene_region: str = "upstream", left_distance: int = 5000, right_distance: int = 5000, window_size: int = 100, hypermethylated: bool = True, gene_hypermethylated_min: int = 20, window_hypermethylated_min: int = 5, min_hypermethylated_windows: int = 5, hypomethylated: bool = True, gene_hypomethylated_max: int = -20, window_hypomethylated_max: int = -5, min_hypomethylated_windows: int = 5, position_count: int = 5, gtf_file: str= "gencode.v41.chr_patch_hapl_scaff.annotation.gtf"):
        """ 
        
        Required Columns: 
        
        position_data: ["chrom", "chromStart", "diff"]
        gene_data: ONE of the following: [("intron", "intron_diff"), ("exon", "exon_diff"), ("upstream", "upstream_diff")]
        
        """
        
        assert isinstance(gene_region, str) and gene_region in ["intron", "exon", "upstream"], "Invalid gene region. Must be in [""intron"", ""exon"", ""upstream""]"
        self.assert_required_columns(position_data, ["chrom", "chromStart", "diff"]), "Invalid position data columns. Must be: [""chrom"", ""chromStart"", ""diff""]"
        self.assert_one_of_column_pairs(gene_data, [("intron", "intron_diff"), ("exon", "exon_diff"), ("upstream", "upstream_diff")]), "Invalid gene data columns. Must be one of: [(""intron"", ""intron_diff""), (""exon"", ""exon_diff""), (""upstream"", ""upstream_diff"")]"
        assert isinstance(gene_data, pd.DataFrame), "List input not acceptable for this function."
        assert isinstance(position_data, pd.DataFrame), "List input not acceptable for this function."

        gene_locations = pl.read_csv(gtf_file, 
                            skip_rows=5, 
                            separator="\t", 
                            has_header=False,
                            columns= [0,2,3,6,8],
                            new_columns=["chrom", "gene_type", "chromStart", "strand", "gene_info"]
                        ) \
                        .filter((pl.col("gene_type") == "gene") & (pl.col("chrom").str.contains("chr"))) \
                        .select(pl.exclude(["gene_type"])) 
        
        # keep gene_type
        
        # get just the gene name and gene type from the gene info string
        gene_locations = gene_locations.with_columns(pl.col("gene_info").str.split(";"))
        

        gene_locations = gene_locations.with_columns(gene_type=pl.col("gene_info").list.get(1)) \
               .with_columns(pl.col("gene_type").str.slice(12)) \
               .with_columns(pl.col("gene_type").str.head(-1))
        

        gene_locations = gene_locations.with_columns(gene=pl.col("gene_info").list.get(2)) \
               .with_columns(pl.col("gene").str.slice(12)) \
               .with_columns(pl.col("gene").str.head(-1))
        
        gene_locations = gene_locations.with_columns(pl.col("chromStart")-1)

        
        filter = (pl.col(gene_region) >= position_count)
        if hypermethylated and hypomethylated:
            filter &= ((pl.col(gene_region+"_diff") > gene_hypermethylated_min) | (pl.col(gene_region+"_diff") < gene_hypomethylated_max))
        elif hypermethylated:
            filter &= (pl.col(gene_region+"_diff") > gene_hypermethylated_min)
        elif hypomethylated:
            filter &= (pl.col(gene_region+"_diff") < gene_hypomethylated_max)
        significant_genes = pl.from_pandas(gene_data.reset_index()).filter(filter).sort(gene_region+"_diff", descending=True)
        merged = significant_genes.join(gene_locations, how="left", on="gene").drop_nulls().select(pl.col(['gene', 'chrom', 'chromStart', 'strand', 'gene_info']))
        # return merged 
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
            diffs = loop_positions["diff"].to_numpy()

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
    
    def graph_upstream_gene_methylation(self, position_data: InputProcessor.data_container, gene_data: InputProcessor.data_container, gene_region: str = "upstream", png_name: str = "upstream_methylation.png", csv_name: str = "upstream_methylation.csv", csv: Optional[str] = None, left_distance: int = 5000, right_distance: int = 5000, window_size: int = 100, hypermethylated: bool = True, gene_hypermethylated_min: int = 20, window_hypermethylated_min: int = 5, min_hypermethylated_windows: int = 5, hypomethylated: bool = True, gene_hypomethylated_max: int = -20, window_hypomethylated_max: int = -5, min_hypomethylated_windows: int = 5, position_count: int = 5, clamp_positive: int = 50, clamp_negative: int = -50, title:str = None, gtf_file: str= "gencode.v41.chr_patch_hapl_scaff.annotation.gtf"):
        """ 

        Required Columns: 

        position_data: ["chrom", "chromStart", "diff"]
        gene_data: ONE of the following: [("intron", "intron_diff"), ("exon", "exon_diff"), ("upstream", "upstream_diff")]

        """

        if position_data is not None:
            assert isinstance(position_data, pd.DataFrame), "List input not acceptable for this function."
        if gene_data is not None:
            assert isinstance(gene_data, pd.DataFrame), "List input not acceptable for this function."
        assert (position_data is not None and gene_data is not None) or (position_data is None and gene_data is None and csv is not None), "Either provide position_data and gene_data or a csv file."
        # column checks are done in __prepare_gene_methylation

        if csv is not None:
            df = pd.read_csv(csv, index_col=0)
            gene_names = list(df.index.values)
            results = df.to_numpy()
        else:
            results, gene_names = self.__prepare_gene_methylation(position_data, gene_data, gene_region, left_distance, right_distance, window_size, hypermethylated, gene_hypermethylated_min, window_hypermethylated_min, min_hypermethylated_windows, hypomethylated, gene_hypomethylated_max, window_hypomethylated_max, min_hypomethylated_windows, position_count, gtf_file)
        
        
        # TODO move this into the normal function, no need for this.. 
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
        cbar_limit = max(abs(results.min()), abs(results.max()))
        # cg = seaborn.clustermap(df.iloc[:,:df.shape[1]//2], col_cluster=False, cmap=cmap, vmin=-1*cbar_limit, center=0, vmax=cbar_limit, method='ward', figsize=(10,max(10, df.shape[0]*0.02))) #TODO add parameter
        
        left_bins = (left_distance)//window_size
        # cg = seaborn.clustermap(df.iloc[:,:left_bins], col_cluster=False, cmap=cmap, vmin=-1*cbar_limit, center=0, vmax=cbar_limit, method='ward', figsize=(10,max(10, df.shape[0]*-0.033 + 80)))
        
        # cg = seaborn.clustermap(df.iloc[:,:left_bins], col_cluster=False, cmap=cmap, vmin=-1*cbar_limit, center=0, vmax=cbar_limit, method='ward', figsize=(10,max(10, df.shape[0]*-0.012 + 49)))
        cg = seaborn.clustermap(df.iloc[:,:left_bins], col_cluster=False, cmap=cmap, vmin=-1*cbar_limit, center=0, vmax=cbar_limit, method='ward', figsize=(10,max(10, df.shape[0]*0.035)))
        cg.cax.set_visible(False)
        ax = cg.ax_heatmap
        idxs = cg.dendrogram_row.reordered_ind
        df = df.iloc[idxs]
        seaborn.heatmap(df,ax=ax,cmap=cmap)
        ax.set_xticks([0] + [((left_distance+right_distance)//window_size)//2 + 0.5] + [((left_distance+right_distance)//window_size)], labels=[f'-{f"{left_distance//1000}k" if left_distance >= 1000 else left_distance}bp', 'Gene Start', f'{f"{right_distance//1000}k" if right_distance >= 1000 else right_distance}bp'])
        ax.tick_params(axis='x', labelrotation=0)
        # num_bins = (left_distance+right_distance)//window_size
        ax.axvline(x = left_bins + 0.5, color = 'black', linestyle = '-') 
        ax.get_yaxis().set_visible(False)
        if title is not None:
            ax.set_title(title, fontsize=15)
        else:
            subtext = ""
            if hypermethylated and hypomethylated:
                subtext = "\n(Hypermethylated and Hypomethylated)"
            elif not hypermethylated and hypomethylated:
                subtext = "\n(Hypomethylated)"
            elif hypermethylated and not hypomethylated:
                subtext = "\n(Hypermethylated)"
            ax.set_title(f"Average Methylation Difference around Significant Gene Start Sites{subtext}", fontsize=15)
        plt.tight_layout()
        
        plt.savefig(png_name, dpi=300, bbox_inches='tight')
        df.to_csv(csv_name)

    def graph_average_upstream_gene_methylation(self, position_data: InputProcessor.data_container, gene_data: InputProcessor.data_container, gene_region: str = "upstream", png_name="average_upstream_methylation.png", csv_name: str = "average_upstream_methylation.csv", csv: Optional[str] = None, left_distance: int = 5000, right_distance: int = 5000, window_size: int = 100, hypermethylated: bool = True, gene_hypermethylated_min: int = 20, window_hypermethylated_min: int = 5, min_hypermethylated_windows: int = 5, hypomethylated: bool = True, gene_hypomethylated_max: int = -20, window_hypomethylated_max: int = -5, min_hypomethylated_windows: int = 5, position_count: int = 5, clamp_positive: int = 50, clamp_negative: int = -50, title:str = None, y_label:str = None, gtf_file: str= "gencode.v41.chr_patch_hapl_scaff.annotation.gtf"):        
        if csv is not None:
            df = pl.read_csv(csv)
        else:
            results, gene_names = self.__prepare_gene_methylation(position_data, gene_data, gene_region, left_distance, right_distance, window_size, hypermethylated, gene_hypermethylated_min, window_hypermethylated_min, min_hypermethylated_windows, hypomethylated, gene_hypomethylated_max, window_hypomethylated_max, min_hypomethylated_windows, position_count, gtf_file)
        
            results = results.tolist()

            df = pl.DataFrame(results, orient="row")

            df = df.with_columns(pl.Series(name="index", values=gene_names))

        means = np.array(df.select(pl.exclude("index").mean()).rows()[0])
        deviations = np.array(df.select(pl.exclude("index").std()).rows()[0])

        plt.clf()

        if title is None:
            text = "Hypomethylated" if hypomethylated else ("Hypermethylated" if hypermethylated else "")
            text = "Hypomethylated and Hypermethylated" if hypermethylated and hypomethylated else text
            plt.title(f"Average of {text}\nGene Regions")
        else:
            plt.title(title)
        plt.ylim(-100,100)
        plt.xticks([])
        plt.errorbar(range(len(means)), means, yerr=deviations, fmt='o--', ms=4, color = "orange")

        if y_label is not None:
            plt.ylabel(y_label)
        else:
            plt.ylabel("AVG(case - control)%")

        plt.savefig(png_name, dpi=300, bbox_inches='tight')

    def graph_upstream_UCSC(self, gene_name: str, position_data: InputProcessor.data_container, name: str="UCSC_graph.bedGraph", before_tss: int = 5000, gtf_file: str = "gencode.v41.chr_patch_hapl_scaff.annotation.gtf"):
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
                file.write(f"track type=bedGraph name=\"{gene} {key}\" visibility=full color={"200,100,0" if case_ctr_str == "ctr" else "0,100,200" if case_ctr_str == "case" else "200,0,100"} altColor=0,100,200\n")

                positions = position_data[["chrom", "chromStart", key]]
                
                positions = positions[(positions["chrom"] == chr) & (positions["chromStart"] >= start) & (positions["chromStart"] < end)]
                
                for _, row in positions.iterrows():
                    _, inner_start, percent = row
                    file.write(f"{chr}\t{inner_start}\t{inner_start+1}\t{percent}\n")

    def graph_full_gene(self, gene_name:str, position_data: InputProcessor.data_container, name="gene_methylation_graph.png", before_tss: int = 0, after_tss: Optional[int] = None, bin_size: int = 500, start_marker: bool = True, end_marker: bool = True, deviation_display: bool = True, legend_size:int = 12, title: str = None, x_label:str = None, y_label:str=None, case_name: str = "Case", ctr_name: str = "Control", gtf_file: str = "gencode.v41.chr_patch_hapl_scaff.annotation.gtf"):        
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

        cols = df.filter(like="blockSizes").columns

        # Aggregate by bins
        result = df.groupby('position_bin', as_index=False).agg({ x : ['mean', 'std'] for x in cols})

        # print(result["blockSizes_case_0"])
        print(result)
        plt.clf()
        
        fig, ax = plt.subplots(figsize=(15, 4))
        
        for x in cols:

            # label = "WT" if "ctr" in x else "KO"
            label = ctr_name if "ctr" in x else case_name
            color = "orange" if "ctr" in x else "blue"
            # 2 4
            ax.plot(result["position_bin"], result[x]["mean"]/100, 'o--', linewidth=3, ms=6, label=label, color = color)
            if deviation_display:
                ax.fill_between(result["position_bin"], 
                                (result[x]["mean"] - result[x]["std"])/100, 
                                (result[x]["mean"] + result[x]["std"])/100, 
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
        
        plt.savefig(name, dpi=300, bbox_inches='tight')
       
    def pie_chart(self, gene_data: InputProcessor.data_container, CCRE_data: InputProcessor.data_container, regions: list[str]|str = ["intron", "exon", "CCRE", "upstream"], name: str = "pie_chart.png", hypermethylated:bool = True, hypomethylated:bool = True, hypermethylated_title:str = None, hypermehylated_min:float = 20, hypomethylated_max:float = -20, hypomethylated_title:str = None, title:str = None):
        assert isinstance(gene_data, pd.DataFrame), "List input not acceptable for this function."
        
        if "intron" in regions:
            self.assert_required_columns(gene_data, ["intron", "intron_diff"])
        if "exon" in regions:
            self.assert_required_columns(gene_data, ["exon", "exon_diff"])
        if "upstream" in regions:
            self.assert_required_columns(gene_data, ["upstream", "upstream_diff"])
        if "CCRE" in regions:
            self.assert_required_columns(CCRE_data, ["CCRE", "CCRE_diff", "tag"])

        def __gen_label(gene_count, pie_chart):
            res_labels=[x[0].upper() + x[1:] for x in gene_count.keys()]
            sum = np.sum(list(gene_count.values()))
            percents=[x/sum for x in gene_count.values()]
            labels = [(x[0] + " - " + f"{x[1]:.1%}") for x in zip(res_labels, percents)]
            labels, slices, _ = zip(*sorted(zip(labels, pie_chart[0], percents), key=lambda x: x[2], reverse=True))
            return labels, slices

        if type(regions) == str:
            regions = [regions]
        hypermethylated_gene_count = {}
        hypomethylated_gene_count = {}
        for region in regions:

            if region == "CCRE":
                df = CCRE_data

                hypermethylated_counts = df[df["CCRE_diff"] > hypermehylated_min]["tag"]
                hypomethylated_counts = df[df["CCRE_diff"] < hypomethylated_max]["tag"]

                for x in ["enhancer", "promoter", "other"]:
                    hypermethylated_gene_count[x] = 0
                    hypomethylated_gene_count[x] = 0
                    if x == "other":
                        hypermethylated_gene_count[x] = hypermethylated_counts[~hypermethylated_counts.str.contains("enhancer|promoter")].count()
                        hypomethylated_gene_count[x] = hypomethylated_counts[~hypomethylated_counts.str.contains("enhancer|promoter")].count()
                    else:
                        hypermethylated_gene_count[x] = hypermethylated_counts.str.contains(x).sum()
                        hypomethylated_gene_count[x] = hypomethylated_counts.str.contains(x).sum()
            else:
                hypermethylated_gene_count[region] = 0
                hypomethylated_gene_count[region] = 0
                df = gene_data
                hypermethylated_counts = df[df[region + "_diff"] > hypermehylated_min][region]
                hypomethylated_counts = df[df[region + "_diff"] < hypomethylated_max][region]

                hypermethylated_gene_count[region] = hypermethylated_counts.sum()
                hypomethylated_gene_count[region] = hypomethylated_counts.sum()

        print(hypermethylated_gene_count)
        print(hypomethylated_gene_count)

        if hypermethylated and hypomethylated:
            fig, ax = plt.subplots(1, 2, figsize=(8, 5))

            # Hypermethylated pie chart
            res_hyper = ax[0].pie(hypermethylated_gene_count.values(), colors=plt.cm.Paired.colors[1::2])

            labels, slices = __gen_label(hypermethylated_gene_count, res_hyper)
            ax[0].legend(slices, labels, loc='lower center', bbox_to_anchor=(0.5,-0.1), ncol=2, fontsize=8)
            
            if hypermethylated_title is not None:
                ax[0].set_title(hypermethylated_title)
            else:
                ax[0].set_title("Hypermethylated")

            # Hypomethylated pie chart
            res_hypo = ax[1].pie(hypomethylated_gene_count.values(), colors=plt.cm.Paired.colors[::2])

            labels, slices = __gen_label(hypomethylated_gene_count, res_hypo)
            ax[1].legend(slices, labels, loc='lower center', bbox_to_anchor=(0.5,-0.1), ncol=2, fontsize=8)
            
            if hypomethylated_title is not None:
                ax[1].set_title(hypomethylated_title)
            else:
                ax[1].set_title("Hypomethylated")
        else:
            if hypermethylated:
                data = hypermethylated_gene_count
                if hypermethylated_title is not None:
                    title = hypermethylated_title
                else:
                    title = "Hypermethylated"
            elif hypomethylated:
                data = hypomethylated_gene_count
                if hypomethylated_title is not None:
                    title = hypomethylated_title
                else:
                    title = "Hypomethylated"

            res = plt.pie(data.values(), colors=plt.cm.Paired.colors[1::2])
            
            labels, slices = __gen_label(data, res)
            plt.legend(slices, labels, loc='lower center', bbox_to_anchor=(0.5,-0.1), ncol=3, fontsize=8)

            plt.title(title)

        if title is not None:
            plt.suptitle(title)
        else:
            plt.suptitle("Distribution of Significant Positions\nin Each Gene Region",fontsize=16)

        plt.tight_layout()
        
        plt.savefig(name, dpi=300, bbox_inches='tight')

