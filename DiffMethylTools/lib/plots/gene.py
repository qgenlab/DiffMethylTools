import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import Optional, Union
import polars as pl
import random
from collections import deque
from matplotlib.colors import LinearSegmentedColormap
import seaborn
from matplotlib.ticker import MaxNLocator


from ..input_processor import InputProcessor

class GeneMixin():
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
    def _prepare_gene_methylation(self, position_data: InputProcessor.data_container,region_data: InputProcessor.data_container, left_distance: int = 1000, right_distance: int = 100, window_size: int = 100, hypermethylated: bool = True, gene_hypermethylated_min: int = 20, window_hypermethylated_min: int = 5, min_hypermethylated_windows: int = 5, hypomethylated: bool = True, gene_hypomethylated_max: int = -20, window_hypomethylated_max: int = -5, min_hypomethylated_windows: int = 5, position_count: int = 5, gtf_file: str= "gencode.v41.chr_patch_hapl_scaff.annotation.gtf", position_or_region: str = "region"):        
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
        results, gene_names = self._prepare_gene_methylation(position_data, region_data, left_distance, right_distance, window_size, hypermethylated, gene_hypermethylated_min, window_hypermethylated_min, min_hypermethylated_windows, hypomethylated, gene_hypomethylated_max, window_hypomethylated_max, min_hypomethylated_windows, position_count, gtf_file, position_or_region)
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
        left_bins = (left_distance)//window_size
        fig_height = max(12, df.shape[0] * 0.1) # max(10, df.shape[0] * 0.035)
        cg = seaborn.clustermap(df.iloc[:,:left_bins], col_cluster=False, cmap=cmap, vmin=-1*cbar_limit, center=0, vmax=cbar_limit, method='ward', figsize=(10,fig_height))
        cg.cax.set_visible(False)
        cg.ax_col_dendrogram.set_visible(False)
        ax = cg.ax_heatmap
        idxs = cg.dendrogram_row.reordered_ind
        df = df.iloc[idxs]
        seaborn.heatmap(df,ax=ax,cmap=cmap, vmin=-1*cbar_limit, center=0, vmax=cbar_limit)
        ax.set_xticks([0, (left_distance + right_distance)//window_size])
        ax.set_xticklabels([
        f'-{left_distance//1000}Kbp' if left_distance >= 1000 else f'-{left_distance}bp',
        f'{right_distance//1000}Kbp' if right_distance >= 1000 else f'{right_distance}bp'
        ])
        ax.tick_params(axis='x', labelrotation=0, labelsize=25)
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


