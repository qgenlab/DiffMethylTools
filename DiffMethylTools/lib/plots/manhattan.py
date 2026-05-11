import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import Optional, Union

# Since plots/ is a subfolder, you need to go up one level (..) to import InputProcessor for your type hints
from ..input_processor import InputProcessor

class ManhattanMixin():
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
