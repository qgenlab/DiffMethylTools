import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import Optional, Union

# Since plots/ is a subfolder, you need to go up one level (..) to import InputProcessor for your type hints
from ..input_processor import InputProcessor



class VolcanoMixin():
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
