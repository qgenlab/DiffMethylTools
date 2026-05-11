import seaborn
from ..input_processor import InputProcessor
from ..analysis import Analysis
from collections import defaultdict, Counter
import random
from matplotlib.figure import Figure
from pathlib import Path
import matplotlib.pyplot as plt

from .volcano import VolcanoMixin
from .annotation import AnnotationMixin
from .coverage import CoverageMixin
from .curves import CurvesMixin
from .gene import GeneMixin
from .manhattan import ManhattanMixin

# class Plots(VolcanoMixin, AnnotationMixin, CoverageMixin, CurvesMixin, GeneMixin, ManhattanMixin):
#     def __init__(self):
#         pass
# 
#     def assert_required_columns(self, df, required_columns):
#         Analysis().assert_required_columns(df, required_columns)
#     
#     def assert_one_of_column_pairs(self, df, column_pairs):
#         Analysis().assert_one_of_column_pairs(df, column_pairs)


class Plots(VolcanoMixin, AnnotationMixin, CoverageMixin, CurvesMixin, GeneMixin, ManhattanMixin):
    def __init__(self, output_dir="plots"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._patch_matplotlib()
    def _patch_matplotlib(self):
        outdir = self.output_dir
        if hasattr(plt, "_DiffMethylTools_patched"):
            return
        original_plt_savefig = plt.savefig
        original_fig_savefig = Figure.savefig
        def custom_plt_savefig(fname, *args, **kwargs):
            fname = Path(fname)
            if fname.parts and fname.parts[0] == outdir.name: path = fname
            else: path = outdir / fname
            path.parent.mkdir(parents=True, exist_ok=True)
            print(f"Saving figure to: {path}")
            return original_plt_savefig(path, *args, **kwargs)
        def custom_fig_savefig(fig, fname, *args, **kwargs):
            fname = Path(fname)
            if fname.parts and fname.parts[0] == outdir.name: path = fname
            else: path = outdir / fname
            path.parent.mkdir(parents=True, exist_ok=True)
            print(f"Saving figure to: {path}")
            return original_fig_savefig(fig, path, *args, **kwargs)
        plt.savefig = custom_plt_savefig
        Figure.savefig = custom_fig_savefig
        plt._DiffMethylTools_patched = True
    def assert_required_columns(self, df, required_columns):
        Analysis().assert_required_columns(df, required_columns)
    def assert_one_of_column_pairs(self, df, column_pairs):
        Analysis().assert_one_of_column_pairs(df, column_pairs)
