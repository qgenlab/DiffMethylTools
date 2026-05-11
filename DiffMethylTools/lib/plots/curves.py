import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import Optional, Union
import polars as pl


from ..input_processor import InputProcessor


class CurvesMixin():
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
        bin[-1] =  df.iloc[-1]["chromStart"]
        df['position_bin'] = bin
            
        print(df)

        if aggregate_samples:
            df["avg_ctr"] = df.filter(like="blockSizes_ctr").mean(axis=1)
            df["avg_case"] = df.filter(like="blockSizes_case").mean(axis=1)
            print(df)

        cols = df.filter(like="blockSizes").columns if not aggregate_samples else ["avg_ctr", "avg_case"]
        result = df.groupby('position_bin', as_index=False).agg({ x : [np.nanmean, np.nanstd] for x in cols})
        print(result)
        plt.clf()
        
        fig, ax = plt.subplots(figsize=(15, 4)) # figsize=(15, 4)
        
        for x in cols:
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
        try:
            repeat_regions_df = pl.read_csv(repeat_regions_df, separator='\t', has_header=False, columns=repeat_regions_columns, new_columns=['chromosome', 'start', 'end', 'labels'], infer_schema_length=10000).to_pandas()
        except pl.exceptions.NoDataError:
            repeat_regions_df = pd.DataFrame(columns=['chromosome', 'start', 'end', 'labels'])
        try:
            enhancer_promoter_df = pl.read_csv(enhancer_promoter_df, separator='\t', has_header=False, columns=enhancer_promoter_columns, new_columns=['chromosome', 'start', 'end', 'cat', 'id'], infer_schema_length=10000).to_pandas()
        except pl.exceptions.NoDataError:
            enhancer_promoter_df = pd.DataFrame(columns=['chromosome', 'start', 'end', 'cat', 'id'])
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
