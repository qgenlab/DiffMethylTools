import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import Optional, Union


from ..input_processor import InputProcessor


class CoverageMixin():
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
