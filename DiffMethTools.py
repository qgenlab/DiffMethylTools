from typing import Callable, Optional, Type, get_args, get_origin, Union
from types import UnionType
from lib import InputProcessor, Analysis, Plots, FormatDefinition
from scipy.stats import mannwhitneyu
import inspect
import argparse
import pandas as pd

class DiffMethTools():
    def __init__(self, pipeline=True):
        self.pipeline = pipeline
        self.obj = Analysis()
        self.plots = Plots()
        self.saved_results = {}

    def __prepare_parameters(self, parameters, **kwargs):
        parameters.pop("self")
        for key in parameters:
            if key in kwargs.keys():
                parameters[key] = kwargs[key]
        
        return parameters

    MERGE_TABLES_REQUIRED_COLUMNS = {
        "case_data": ["chromosome", "position_start", "coverage", "methylation_percentage", "positive_methylation_count", "negative_methylation_count"],
        "ctr_data": ["chromosome", "position_start", "coverage", "methylation_percentage", "positive_methylation_count", "negative_methylation_count"]
    }
    def merge_tables(self, case_data: InputProcessor, ctr_data: InputProcessor, min_cov = 10, cov_percentile = 1.0, min_samp_ctr = 2, min_samp_case = 2) -> pd.DataFrame:
        """Merge case and control data tables.

        .. note::
            ``case_data`` and ``ctr_data`` must be from a list of samples, where each contains one of the following column formats:
                - ``["chromosome", "position_start", "coverage", "methylation_percentage"]``
                - ``["chromosome", "position_start", "positive_methylation_count", "negative_methylation_count"]``
                - Additionally, ``"position_end"`` and ``"strand"`` columns can be used. These will be considered during the join process. 

        :param case_data: The case data to be merged.
        :type case_data: InputProcessor
        :param ctr_data: The control data to be merged.
        :type ctr_data: InputProcessor
        :param min_cov: Minimum coverage filter, defaults to 10
        :type min_cov: int, optional
        :param cov_percentile: Maximum coverage filter (percentile of sample coverage), defaults to 1.0
        :type cov_percentile: float, optional
        :param min_samp_ctr: Minimum samples in control, defaults to 2
        :type min_samp_ctr: int, optional
        :param min_samp_case: Minimum samples in case, defaults to 2
        :type min_samp_case: int, optional
        :return: Merged data table
        :rtype: pd.DataFrame
        """
        parameters = locals().copy()

        case_data = case_data.copy()
        ctr_data = ctr_data.copy()

        case_data.process()
        ctr_data.process()

        parameters = self.__prepare_parameters(parameters, case_data=case_data.data_container, ctr_data=ctr_data.data_container)

        res = self.obj.merge_tables(**parameters)

        if self.pipeline:
            self.saved_results[self.merge_tables] = res

        return res
    
    WINDOW_BASED_REQUIRED_COLUMNS = {
        "data": ["chromosome", "position_start", "methylation_percentage*", "avg_case", "avg_ctr"]
    }
    def window_based(self, data: Optional[InputProcessor] = None, statistical_test : Callable = mannwhitneyu, window = 1000, step = 500, min_nbr_per_win = 5, processes=12, min_std=0.1) -> pd.DataFrame:
        """Perform window-based analysis.
        
        .. note::
            ``data`` must contain the following column format:
                - ``["chromosome", "position_start", "methylation_percentage*", "avg_case", "avg_ctr"]``
                - There must be a ``methylation_percentage`` column for each sample to test.

        .. note::
            ``statistical_test`` must be a function from the `scipy.stats <https://docs.scipy.org/doc/scipy/reference/stats.html>`_ package. Accepted functions include:
                - ``mannwhitneyu``
                - ``ttest_ind``
                - ``ranksums``
                - ``ks_2samp``
                - ``median_test``
                - ``brunnermunzel``

        :param data: Input data. Not necessary if the pipeline is in use, defaults to None
        :type data: InputProcessor, optional
        :param statistical_test: Statistical test to use, must be from `scipy.stats <https://docs.scipy.org/doc/scipy/reference/stats.html>`_ package, defaults to ``mannwhitneyu``
        :type statistical_test: Callable, optional
        :param window: Sample window size (in bp), defaults to 1000
        :type window: int, optional
        :param step: Window step size, defaults to 500
        :type step: int, optional
        :param min_nbr_per_win: Minimum number of sample positions per window, defaults to 5
        :type min_nbr_per_win: int, optional
        :param processes: Number of CPU processes, defaults to 12
        :type processes: int, optional
        :param min_std: Minimum standard deviation filter, defaults to 0.1
        :type min_std: float, optional
        :return: Window-based analysis results
        :rtype: pd.DataFrame
        """
        assert (not self.pipeline and data is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 

        parameters = locals().copy()

        if data is not None:
            data = data.copy()
            data.process()
            data = data.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            data = self.saved_results[self.merge_tables]
            
        parameters = self.__prepare_parameters(parameters, data=data)

        res = self.obj.window_based(**parameters)

        if self.pipeline:
            self.saved_results[self.window_based] = res

        return res
    
    POSITION_BASED_REQUIRED_COLUMNS = {
        "data": ["chromosome", "position_start", "methylation_percentage*"]
    }
    def position_based(self, data: Optional[InputProcessor] = None , method="limma", features=None, test_factor="Group", processes=12, model="eBayes", min_std=0.1, fill_na:bool=True) -> pd.DataFrame:
        """Perform position-based analysis. Has options for using the gamma function, or the limma R package.

        .. note::
            ``data`` must contain the following column format:
                - ``["chromosome", "position_start", "methylation_percentage*"]``
                - There must be a ``methylation_percentage`` column for each sample to test.

        :param data: Input data. Not necessary if the pipeline is in use, defaults to None
        :type data: InputProcessor, optional
        :param method: Position-based method to use, options are ``["gamma", "limma"]``, defaults to "limma"
        :type method: str, optional
        :param features: Features .csv file used with the "limma" method, has indicies with the names of samples in ``data`` and columns for the features, defaults to None
        :type features: str, optional
        :param test_factor: Test factor used with the "limma" method, defaults to "Group"
        :type test_factor: str, optional
        :param processes: Number of CPU processes, defaults to 12
        :type processes: int, optional
        :param model: Limma model to use, options are ``["eBayes", "treat"]`` defaults to "eBayes"
        :type model: str, optional
        :param min_std: Minimum standard deviation filter, defaults to 0.1
        :type min_std: float, optional
        :param fill_na: Fill NA values with row group average, defaults to True
        :return: Position-based analysis results
        :rtype: pd.DataFrame
        """
        assert (not self.pipeline and data is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided."

        parameters = locals().copy()

        if data is not None:
            data = data.copy()
            data.process()
            data = data.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            data = self.saved_results[self.merge_tables]
            
        parameters = self.__prepare_parameters(parameters, data=data)

        res = self.obj.position_based(**parameters)

        if self.pipeline:
            self.saved_results[self.position_based] = res

        return res

    MAP_WIN_2_POS_REQUIRED_COLUMNS = {
        "window_data": ["chromosome", "region_start", "region_end"],
        "position_data": ["chromosome", "position_start", "avg_case", "avg_ctr"]
    }
    def map_win_2_pos(self, window_data: Optional[InputProcessor] = None, position_data: Optional[InputProcessor] = None, processes=12, sub_window_size = 100, sub_window_step = 100, sub_window_min_diff=0, pipeline_window_result: str = "auto") -> pd.DataFrame:
        """Map windows to positions.

        .. note::
            Required columns for ``window_data``:
                - ``["chromosome", "region_start", "region_end"]``

            Required columns for ``position_data``:
                - ``["chromosome", "position_start", "avg_case", "avg_ctr"]``

        :param window_data: Window data to map with. Not necessary if the pipeline is in use, defaults to None
        :type window_data: InputProcessor, optional
        :param position_data: Position data to map the windows to. Not necessary if the pipeline is in use, defaults to None
        :type position_data: InputProcessor, optional
        :param processes: Number of CPU processes, defaults to 12
        :type processes: int, optional
        :param sub_window_size: Sub-window size for a deeper difference filtering, defaults to 100
        :type sub_window_size: int, optional
        :param sub_window_step: Sub-window step size, defaults to 100
        :type sub_window_step: int, optional
        :param sub_window_min_diff: Sub-window minimum difference, defaults to 0
        :type sub_window_min_diff: int, optional
        :param pipeline_window_result: The function results to use as input if DiffMethTools is pipelined and no data is provided. Options are ``["auto", "filters", "generate_q_values", "window_based"]``, defaults to "auto"
        :type pipeline_window_result: str, optional
        :return: Mapped windows to positions
        :rtype: pd.DataFrame
        """
        assert (not self.pipeline and window_data is not None and position_data is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert pipeline_window_result in ["auto", "filters", "generate_q_values", "window_based"], "Invalid parameter for pipeline_window_result. Options are: [""auto"", ""filters"", ""generate_q_values"", ""window_based""]"

        parameters = locals().copy()

        if window_data is not None:
            window_data = window_data.copy()
            window_data.process()
            window_data = window_data.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            # TODO after adding the DMR function, add option for which data to use
            if pipeline_window_result == "auto":
                if self.saved_results.get(self.filters) is not None:
                    window_data = self.saved_results[self.window_based]
                elif self.saved_results.get(self.generate_q_values) is not None:
                    window_data = self.saved_results[self.generate_q_values]
                elif self.saved_results.get(self.window_based) is not None:
                    window_data = self.saved_results[self.window_based]
            elif pipeline_window_result == "filters":
                window_data = self.saved_results[self.filters]
            elif pipeline_window_result == "generate_q_values":
                window_data = self.saved_results[self.generate_q_values]
            elif pipeline_window_result == "window_based":
                window_data = self.saved_results[self.window_based]
        if position_data is not None:
            position_data = position_data.copy()
            position_data.process()
            position_data = position_data.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            position_data = self.saved_results[self.merge_tables]
        
        parameters.pop("pipeline_window_result")
        parameters = self.__prepare_parameters(parameters, window_data=window_data, position_data=position_data)

        res = self.obj.map_win_2_pos(**parameters)

        if self.pipeline:
            self.saved_results[self.map_win_2_pos] = res

        return res

    GENERATE_Q_VALUES_REQUIRED_COLUMNS = {
        "data": ["p-val"]
    }
    def generate_q_values(self, data: Optional[InputProcessor] = None, method: str = "fdr_bh", position_or_window: str = "auto") -> pd.DataFrame:
        """Generate q-values (also referred to as adjusted p-values or FDR).

        .. note::
            Required columns for ``data``:
                - ``["p-val"]``

        :param data: Input p-value data. Not necessary if the pipeline is in use, defaults to None
        :type data: InputProcessor, optional
        :param method: P-value adjustment method to use from `statsmodels.stats.multitest.multipletests <https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html#statsmodels.stats.multitest.multipletests-parameters>`_, defaults to "fdr_bh".
        :type method: str, optional
        :param position_or_window: The position-based or window-based results to use as input if DiffMethTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto"
        :type position_or_window: str, optional
        :return: Data with q-values (adjusted p-values or FDR)
        :rtype: pd.DataFrame
        """
        assert (not self.pipeline and data is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert position_or_window in ["auto", "position", "window"], "Invalid parameter for position_or_window. Options are: [""auto"", ""position"", ""window""]"

        parameters = locals().copy()

        if data is not None:
            data = data.copy()
            data.process()
            data = data.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            if position_or_window == "auto":
                if self.saved_results.get(self.window_based) is not None:
                    data = self.saved_results[self.window_based]
                    position_or_window = "window"
                else:
                    data = self.saved_results[self.position_based]
                    position_or_window = "position"
            elif position_or_window == "position":
                data = self.saved_results[self.position_based]
            else:
                data = self.saved_results[self.window_based]
            
            
        parameters.pop("position_or_window")
        parameters = self.__prepare_parameters(parameters, data=data)

        res = self.obj.generate_q_values(**parameters)

        if self.pipeline:
            self.saved_results[(self.generate_q_values, position_or_window)] = res

        return res

    FILTERS_REQUIRED_COLUMNS = {
        "data": ["q-value", "diff"]
    }
    def filters(self, data: Optional[InputProcessor] = None, max_q_value=0.05, abs_min_diff=25, position_or_window: str = "auto") -> pd.DataFrame:
        """Filter data by q-value and minimum difference.

        .. note::
            Required columns for ``data``:
                - ``["q-value", "diff"]``

        :param data: Input data. Not necessary if the pipeline is in use, defaults to None
        :type data: InputProcessor, optional
        :param max_q_value: Maximum q-value filter, defaults to 0.05
        :type max_q_value: float, optional
        :param abs_min_diff: Absolute minimum difference filter, defaults to 25
        :type abs_min_diff: int, optional
        :param position_or_window: The position-based or window-based results to use as input if DiffMethTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto"
        :type position_or_window: str, optional
        :return: Filtered input data
        :rtype: pd.DataFrame
        """
        assert (not self.pipeline and data is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert position_or_window in ["auto", "position", "window"], "Invalid parameter for position_or_window. Options are: [""auto"", ""position"", ""window""]"

        parameters = locals().copy()

        if data is not None:
            data = data.copy()
            data.process()
            data = data.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            if position_or_window == "auto":
                if self.saved_results.get((self.generate_q_values, "window")) is not None:
                    data = self.saved_results[(self.generate_q_values, "window")]
                    position_or_window = "window"
                elif self.saved_results.get((self.generate_q_values, "position")) is not None:
                    data = self.saved_results[(self.generate_q_values, "position")]
                    position_or_window = "position"
                else:
                    # TODO invalid
                    pass
        
        
        parameters.pop("position_or_window")
        parameters = self.__prepare_parameters(parameters, data=data)

        res = self.obj.filters(**parameters)

        if self.pipeline:
            self.saved_results[(self.filters, position_or_window)] = res

        return res

    GENERATE_DMR_REQUIRED_COLUMNS = {
        "significant_position_data": ["chromosome", "position_start", "diff"],
        "position_data": ["chromosome", "position_start", "diff"]
    }
    def generate_DMR(self, significant_position_data: Optional[InputProcessor] = None, position_data: Optional[InputProcessor] = None, min_pos=3, neutral_change_limit=7.5, neutral_perc=30, opposite_perc=10, significant_position_pipeline: str = "auto") -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        # TODO ask chris what to write here for documentation
        """Generate Differentially Methylated Regions (DMRs).
        
        .. note::
            Required columns for ``significant_position_data``:
                - ``["chromosome", "position_start", "diff"]``
            Required columns for ``position_data``:
                - ``["chromosome", "position_start", "diff"]``

        :param significant_position_data: Significant position data. Not necessary if the pipeline is in use, defaults to None
        :type significant_position_data: InputProcessor, optional
        :param position_data: All position data. Not necessary if the pipeline is in use, defaults to None
        :type position_data: InputProcessor, optional
        :param min_pos: Minimum positions, defaults to 3
        :type min_pos: int, optional
        :param neutral_change_limit: Neutral change limit, defaults to 7.5
        :type neutral_change_limit: float, optional
        :param neutral_perc: Neutral percentage, defaults to 30
        :type neutral_perc: int, optional
        :param opposite_perc: Opposite percentage, defaults to 10
        :type opposite_perc: int, optional
        :param significant_position_pipeline: The significant position-based or window-based results to use as input if DiffMethTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto"
        :type significant_position_pipeline: str, optional
        :return: DMR results
        :rtype: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        """
        assert (not self.pipeline and significant_position_data is not None and position_data is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert significant_position_pipeline in ["auto", "position", "window"], "Invalid parameter for significant_position_pipeline. Options are: [""auto"", ""position"", ""window""]"

        parameters = locals().copy()

        if significant_position_data is not None:
            significant_position_data = significant_position_data.copy()
            significant_position_data.process()
            significant_position_data = significant_position_data.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            # can be either auto, position or window
            if significant_position_pipeline == "auto":
                if self.saved_results.get((self.filters, "position")) is not None:
                    significant_position_data = self.saved_results[(self.filters, "position")]
                elif self.saved_results.get(self.map_win_2_pos) is not None:
                    significant_position_data = self.saved_results[self.map_win_2_pos]
                else:
                    # TODO invalid
                    pass
            elif significant_position_pipeline == "position":
                significant_position_data = self.saved_results[(self.filters, "position")]
            else:
                significant_position_data = self.saved_results[self.map_win_2_pos]

        if position_data is not None:
            position_data = position_data.copy()
            position_data.process()
            position_data = position_data.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            position_data = self.saved_results[self.merge_tables]
        
        
        parameters.pop("significant_position_pipeline")
        parameters = self.__prepare_parameters(parameters, significant_position_data=significant_position_data, position_data=position_data)

        res = self.obj.generate_DMR(**parameters)

        if self.pipeline:
            self.saved_results[(self.generate_DMR, "cluster_df")] = res[0]
            self.saved_results[(self.generate_DMR, "unclustered_dms_df")] = res[1]
            self.saved_results[(self.generate_DMR, "clustered_dms_df")] = res[2]
        return res
    
    MAP_POSITIONS_TO_GENES_REQUIRED_COLUMNS = {
        "positions": ["chromosome", "position_start", "diff"]
    }
    def map_positions_to_genes(self, positions: Optional[InputProcessor] = None, gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], min_pos_diff=0, bed_file="outfile_w_hm450.bed", gtf_file="gencode.v41.chr_patch_hapl_scaff.annotation.gtf", pipeline_input_source = "auto") -> tuple[pd.DataFrame, pd.DataFrame]:
        """Map positions to genes.

        .. note::
            Required columns for ``positions``:
                - ``["chromosome", "position_start", "diff"]``

        :param positions: Position data. Not necessary if the pipeline is in use, defaults to None
        :type positions: InputProcessor, optional
        :param gene_regions: Gene regions to map to. Options are any combination of ``["intron", "exon", "upstream", "CCRE"]``, defaults to ``["intron", "exon", "upstream", "CCRE"]``
        :type gene_regions: list[str] | str, optional
        :param min_pos_diff: Minimum position difference for mapping, defaults to 0
        :type min_pos_diff: int, optional
        :param bed_file: BED annotation file with unflexible input, defaults to "outfile_w_hm450.bed"
        :type bed_file: str, optional
        :param gtf_file: GTF annotation file with unflexible input, defaults to "gencode.v41.chr_patch_hapl_scaff.annotation.gtf"
        :type gtf_file: str, optional
        :param pipeline_input_source: Pipeline input source for pipelining, options are ["auto", "map_win_2_pos", "generate_q_values", "filters"], defaults to "auto"
        :type pipeline_input_source: str, optional
        :return: Mapped positions to genes in a list: ``[gene mapping dataframe, CCRE mapping dataframe]``
        :rtype: list[pd.DataFrame] 
        """
        assert (not self.pipeline and positions is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert pipeline_input_source in ["auto", "map_win_2_pos", "generate_q_values", "filters"], "Invalid parameter for pipeline_input_source. Options are: [""auto"", ""map_win_2_pos"", ""generate_q_values"", ""filters""]"

        parameters = locals().copy()

        if positions is not None:
            positions = positions.copy()
            positions.process()
            positions = positions.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            if pipeline_input_source == "auto":
                if self.saved_results.get((self.filters, "position")) is not None:
                    positions = self.saved_results[(self.filters, "position")]
                elif self.saved_results.get((self.generate_q_values, "position")) is not None:
                    positions = self.saved_results[(self.generate_q_values, "position")]
                elif self.saved_results.get(self.map_win_2_pos) is not None:
                    positions = self.saved_results[self.map_win_2_pos]
            elif pipeline_input_source == "map_win_2_pos":
                positions = self.saved_results[self.map_win_2_pos]
            elif pipeline_input_source == "generate_q_values":
                positions = self.saved_results[(self.generate_q_values, "position")]
            else:
                positions = self.saved_results[(self.filters, "position")]
        
        
        parameters.pop("pipeline_input_source")
        parameters = self.__prepare_parameters(parameters, positions=positions)

        res = self.obj.map_positions_to_genes(**parameters)

        if self.pipeline:
            self.saved_results[(self.map_positions_to_genes, "gene")] = res[0]
            self.saved_results[(self.map_positions_to_genes, "CCRE")] = res[1]

        return tuple(x.reset_index() for x in res)
    
    MAP_WINDOWS_TO_GENES_REQUIRED_COLUMNS = {    
        "windows": ["chromosome", "region_start", "region_end", "diff"]
    }
    def map_windows_to_genes(self, windows: Optional[InputProcessor] = None, gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], min_pos_diff=0, bed_file="outfile_w_hm450.bed", gtf_file="gencode.v41.chr_patch_hapl_scaff.annotation.gtf", enhd_thr = 500000, enhp_thr = 50000, prom_thr = 2000, processes=12, pipeline_input_source = "auto") -> tuple[pd.DataFrame, pd.DataFrame]:
        """Map windows to genes.

        .. note::
            Required columns for ``windows``:
                - ``["chromosome", "region_start", "region_end", "diff"]``

        :param windows: Window data. Not necessary if the pipeline is in use, defaults to None
        :type windows: InputProcessor, optional
        :param gene_regions: Gene regions to map to. Options are any combination of ``["intron", "exon", "upstream", "CCRE"]``, defaults to ``["intron", "exon", "upstream", "CCRE"]``
        :type gene_regions: list[str] | str, optional
        :param min_pos_diff: Minimum position difference for mapping, defaults to 0
        :type min_pos_diff: int, optional
        :param bed_file: BED annotation file with unflexible input, defaults to "outfile_w_hm450.bed"
        :type bed_file: str, optional
        :param gtf_file: GTF annotation file with unflexible input, defaults to "gencode.v41.chr_patch_hapl_scaff.annotation.gtf"
        :type gtf_file: str, optional
        :param enhd_thr: Distal enhancer distance threshold for finding nearest gene to CCRE, defaults to 500000
        :type enhd_thr: int, optional
        :param enhp_thr: Proximal enhancer distance threshold for finding nearest gene to CCRE, defaults to 50000
        :type enhp_thr: int, optional
        :param prom_thr: Promoter enhancer distance threshold for finding nearest gene to CCRE, defaults to 2000
        :type prom_thr: int, optional
        :param processes: Number of CPU processes, defaults to 12
        :type processes: int, optional
        :param pipeline_input_source: Pipeline input source for pipelining, options are ["auto", "generate_DMR", "generate_q_values", "filters"], defaults to "auto"
        :type pipeline_input_source: str, optional
        :return: Mapped windows to genes in a list: ``[gene mapping dataframe, CCRE mapping dataframe]``
        :rtype: list[pd.DataFrame]
        """
        assert (not self.pipeline and windows is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert pipeline_input_source in ["auto", "generate_DMR", "generate_q_values", "filters"], "Invalid parameter for pipeline_input_source. Options are: [""auto"", ""generate_DMR"", ""generate_q_values"", ""filters""]"

        parameters = locals().copy()

        if windows is not None:
            windows = windows.copy()
            windows.process()
            windows = windows.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            if pipeline_input_source == "auto":
                if self.saved_results.get(self.generate_DMR) is not None:
                    windows = self.saved_results[(self.generate_DMR, "clustered_dms_df")]
                elif self.saved_results.get((self.filters, "window")) is not None:
                    windows = self.saved_results[(self.filters, "window")]
                elif self.saved_results.get((self.generate_q_values, "window")) is not None:
                    windows = self.saved_results[(self.generate_q_values, "window")]
            elif pipeline_input_source == "generate_DMR":
                windows = self.saved_results[(self.generate_DMR, "clustered_dms_df")]
            elif pipeline_input_source == "generate_q_values":
                windows = self.saved_results[(self.generate_q_values, "window")]
            else:
                windows = self.saved_results[(self.filters, "window")]
        
        
        parameters.pop("pipeline_input_source")

        parameters = self.__prepare_parameters(parameters, windows=windows)

        res = self.obj.map_windows_to_genes(**parameters)

        if self.pipeline:
            self.saved_results[(self.map_windows_to_genes, "gene")] = res[0]
            self.saved_results[(self.map_windows_to_genes, "CCRE")] = res[1]

        return tuple(x.reset_index() for x in res)
    
    VOLCANO_PLOT_REQUIRED_COLUMNS = {
        "data": ["q-value", "diff"]
    }
    def volcano_plot(self, data: Optional[InputProcessor] = None, name : str = "volcano_plot.png", threshold : Optional[float] = 0.05, line: Optional[float] = 15, x_range: tuple[int, int] = (-100, 100), y_max: int = None, title: str = None, x_label: str = None, y_label: str = None, position_or_window: str = "auto") -> None:
        """Generate a volcano plot.
        
        .. note::
            Required columns for ``data``:
                - ``["q-value", "diff"]``

        :param data: Input data. Not necessary if the pipeline is in use, defaults to None
        :type data: InputProcessor, optional
        :param name: Output file name, defaults to "volcano_plot.png"
        :type name: str, optional
        :param threshold: Q-value threshold for horizontal line. Set to None to have no line, defaults to 0.05
        :type threshold: float, optional
        :param line: Vertical line threshold for the `abs(line)` vertical line. Set to NOne to have no lines, defaults to 15
        :type line: float, optional
        :param x_range: X-axis range, defaults to (-100, 100)
        :type x_range: tuple[int, int], optional
        :param y_max: Y-axis maximum, defaults to None
        :type y_max: int, optional
        :param title: Plot title, defaults to None for a generic title
        :type title: str, optional
        :param x_label: X-axis label, defaults to None for a generic label
        :type x_label: str, optional
        :param y_label: Y-axis label, defaults to None for a generic label
        :type y_label: str, optional
        :param position_or_window: The position-based or window-based results to use as input if DiffMethTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto"
        :type position_or_window: str, optional
        """
        assert (not self.pipeline and data is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert position_or_window in ["auto", "position", "window"], "Invalid parameter for position_or_window. Options are: [""auto"", ""position"", ""window""]"

        parameters = locals().copy()

        if data is not None:
            data = data.copy()
            data.process()
            data = data.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            if position_or_window == "auto":
                if self.saved_results.get((self.generate_q_values, "window")) is not None:
                    data = self.saved_results[(self.generate_q_values, "window")]
                elif self.saved_results.get((self.generate_q_values, "position")) is not None:
                    data = self.saved_results[(self.generate_q_values, "position")]
                else:
                    # TODO invalid
                    pass
            elif position_or_window == "position":
                data = self.saved_results[(self.generate_q_values, "position")]
            else:
                data = self.saved_results[(self.generate_q_values, "window")]
            
        parameters.pop("position_or_window")
        parameters = self.__prepare_parameters(parameters, data=data)

        self.plots.volcano_plot(**parameters)

    MANHATTAN_PLOT_REQUIRED_COLUMNS = {
        "data": ["chromosome", "q-value", "position_start", "region_start"]
    }
    def manhattan_plot(self, data: Optional[InputProcessor] = None, name: str = "manhattan_plot.png", threshold: Optional[float] = 0.05, title: str = None, x_label: str = None, y_label: str = None, position_or_window: str = "auto") -> None:
        """Generate a manhattan plot.
        
        .. note::
            Required columns for ``data`` are in one of the following column formats:
                - ``["chromosome", "q-value", "position_start"]``
                - ``["chromosome", "q-value", "region_start"]``

        :param data: Input data. Not necessary if the pipeline is in use, defaults to None
        :type data: InputProcessor, optional
        :param name: Output file name, defaults to "manhattan_plot.png"
        :type name: str
        :param threshold: Q-value threshold for horizontal line. Set to None to have no line, defaults to 0.05
        :type threshold: float, optional
        :param title: Plot title, defaults to None for a generic title
        :type title: str, optional
        :param x_label: X-axis label, defaults to None for a generic label
        :type x_label: str, optional
        :param y_label: Y-axis label, defaults to None for a generic label
        :type y_label: str, optional
        :param position_or_window: The position-based or window-based results to use as input if DiffMethTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto"
        :type position_or_window: str, optional
        """        
        assert (not self.pipeline and data is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert position_or_window in ["auto", "position", "window"], "Invalid parameter for position_or_window. Options are: [""auto"", ""position"", ""window""]"


        parameters = locals().copy()

        if data is not None:
            data = data.copy()
            data.process()
            data = data.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            if position_or_window == "auto":
                if self.saved_results.get((self.generate_q_values, "window")) is not None:
                    data = self.saved_results[(self.generate_q_values, "window")]
                elif self.saved_results.get((self.generate_q_values, "position")) is not None:
                    data = self.saved_results[(self.generate_q_values, "position")]
                else:
                    # TODO invalid
                    pass
            elif position_or_window == "position":
                data = self.saved_results[(self.generate_q_values, "position")]
            else:
                data = self.saved_results[(self.generate_q_values, "window")]
                
        parameters.pop("position_or_window")
        parameters = self.__prepare_parameters(parameters, data=data)
        
        self.plots.manhattan_plot(**parameters)

    COVERAGE_PLOT_REQUIRED_COLUMNS = {
        "case_data": ["coverage", "positive_methylation_count", "negative_methylation_count"],
        "ctr_data": ["coverage", "positive_methylation_count", "negative_methylation_count"]
    }
    def coverage_plot(self, case_data: InputProcessor, ctr_data: InputProcessor, name : str = "coverage_plot.png", cov_min : int = -1, cov_max : int = -1, cov_max_percentile : float = -1, bins:int = 20, title: str = None, x_label: str = None, y_label: str = None) -> None:
        """Generate a coverage plot.

        .. note::
            ``case_data`` and ``ctr_data`` must be from a list of samples, where each contains one of the following column formats:
                - ``["coverage"]``
                - ``["positive_methylation_count", "negative_methylation_count"]``

        :param case_data: The case data to plot.
        :type case_data: InputProcessor
        :param ctr_data: The control data to plot.
        :type ctr_data: InputProcessor
        :param name: Output file name, defaults to "coverage_plot.png"
        :type name: str, optional
        :param cov_min: Minimum coverage display, defaults to -1
        :type cov_min: int, optional
        :param cov_max: Maximum coverage display, defaults to -1
        :type cov_max: int, optional
        :param cov_max_percentile: Maximum coverage percentile display, defaults to -1. Overrides ``cov_max`` if set.
        :type cov_max_percentile: float, optional
        :param bins: Number of bins, defaults to 20
        :type bins: int, optional
        :param title: Plot title, defaults to None for a generic title
        :type title: str, optional
        :param x_label: X-axis label, defaults to None for a generic label
        :type x_label: str, optional
        :param y_label: Y-axis label, defaults to None for a generic label
        :type y_label: str, optional
        """
        parameters = locals().copy()

        case_data = case_data.copy()
        ctr_data = ctr_data.copy()

        case_data.process()
        ctr_data.process()

        parameters = self.__prepare_parameters(parameters, case_data=case_data.data_container, ctr_data=ctr_data.data_container)

        self.plots.coverage_plot(**parameters)

    GRAPH_GENE_REGIONS_REQUIRED_COLUMNS = {
        "gene_data": ["intron", "intron_diff", "exon", "exon_diff", "upstream", "upstream_diff"],
        "ccre_data": ["CCRE", "CCRE_diff"]
    }
    def graph_gene_regions(self, gene_data: Optional[InputProcessor] = None, ccre_data: Optional[InputProcessor] = None, name: str ="gene_regions.png", gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], intron_cutoff: int = -1, exon_cutoff: int = -1, upstream_cutoff: int = -1, CCRE_cutoff: int = -1, title: str = None, x_label: str = None, intron_y_label: str = None, exon_y_label: str = None, upstream_y_label: str = None, CCRE_y_label: str = None, position_or_window: str = "auto") -> None:
        """Generate a graph of gene regions.
        
        .. note::
            Required columns for ``gene_data``:
                - If ``"intron"`` is in ``gene_regions``: ``["intron", "intron_diff"]``
                - If ``"exon"`` is in ``gene_regions``: ``["exon", "exon_diff"]``
                - If ``"upstream"`` is in ``gene_regions``: ``["upstream", "upstream_diff"]``
            Required columns for ``ccre_data``:
                - If ``"CCRE"`` is in ``gene_regions``: ``["CCRE", "CCRE_diff"]`` 

        :param gene_data: Gene data. Not necessary if the pipeline is in use, defaults to None
        :type gene_data: InputProcessor, optional
        :param ccre_data: CCRE data. Not necessary if the pipeline is in use, defaults to None
        :type ccre_data: InputProcessor, optional
        :param name: Output file name, defaults to "gene_regions.png"
        :type name: str, optional
        :param gene_regions: Gene regions to map to. Options are any combination of ``["intron", "exon", "upstream", "CCRE"]``, defaults to ``["intron", "exon", "upstream", "CCRE"]``
        :type gene_regions: list[str] | str, optional
        :param intron_cutoff: Intron count vertical display cutoff, defaults to -1 (for no cutoff)
        :type intron_cutoff: int, optional
        :param exon_cutoff: Exon count vertical display cutoff, defaults to -1 (for no cutoff)
        :type exon_cutoff: int, optional
        :param upstream_cutoff: Upstream count vertical display cutoff, defaults to -1 (for no cutoff)
        :type upstream_cutoff: int, optional
        :param CCRE_cutoff: CCRE count vertical display cutoff, defaults to -1 (for no cutoff)
        :type CCRE_cutoff: int, optional
        :param title: Plot title, defaults to None for a generic title
        :type title: str, optional
        :param x_label: X-axis label, defaults to None for a generic label
        :type x_label: str, optional
        :param intron_y_label: Intron Y-axis label, defaults to None for a generic label
        :type intron_y_label: str, optional
        :param exon_y_label: Exon Y-axis label, defaults to None for a generic label
        :type exon_y_label: str, optional
        :param upstream_y_label: Upstream Y-axis label, defaults to None for a generic label
        :type upstream_y_label: str, optional
        :param CCRE_y_label: CCRE Y-axis label, defaults to None for a generic label
        :type CCRE_y_label: str, optional
        :param position_or_window: The position-based or window-based results to use as input if DiffMethTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto"
        :type position_or_window: str, optional
        """
        assert (not self.pipeline and (gene_data is not None or ccre_data is not None)) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert position_or_window in ["auto", "position", "window"], "Invalid parameter for position_or_window. Options are: [""auto"", ""position"", ""window""]"

        parameters = locals().copy()

        if gene_data is not None:
            gene_data = gene_data.copy()
            gene_data.process()
            gene_data = gene_data.data_container
        else:
            if position_or_window == "auto":
                if self.saved_results.get((self.map_positions_to_genes, "gene")) is not None:
                    gene_data = self.saved_results[(self.map_positions_to_genes, "gene")]
                elif self.saved_results.get((self.map_windows_to_genes, "gene")) is not None:
                    gene_data = self.saved_results[(self.map_windows_to_genes, "gene")]
                else:
                    # TODO invalid
                    pass
            elif position_or_window == "position":
                gene_data = self.saved_results[(self.map_positions_to_genes, "gene")]
            else:
                gene_data = self.saved_results[(self.map_windows_to_genes, "gene")]

        if ccre_data is not None:
            ccre_data = ccre_data.copy()
            ccre_data.process()
            ccre_data = ccre_data.data_container
        else:
            if position_or_window == "auto":
                if self.saved_results.get((self.map_positions_to_genes, "CCRE")) is not None:
                    ccre_data = self.saved_results[(self.map_positions_to_genes, "CCRE")]
                elif self.saved_results.get((self.map_windows_to_genes, "CCRE")) is not None:
                    ccre_data = self.saved_results[(self.map_windows_to_genes, "CCRE")]
                else:
                    # TODO invalid
                    pass
            elif position_or_window == "position":
                ccre_data = self.saved_results[(self.map_positions_to_genes, "CCRE")]
            else:
                ccre_data = self.saved_results[(self.map_windows_to_genes, "CCRE")]
        
        parameters.pop("position_or_window")
        parameters = self.__prepare_parameters(parameters, gene_data=gene_data, ccre_data=ccre_data)

        res = self.plots.graph_gene_regions(**parameters)

    GRAPH_UPSTREAM_GENE_METHYLATION_REQUIRED_COLUMNS = {
        "position_data": ["chromosome", "position_start", "diff"],
        "gene_data": ["intron", "intron_diff", "exon", "exon_diff", "upstream", "upstream_diff"]
    }
    def graph_upstream_gene_methylation(self, position_data: Optional[InputProcessor] = None, gene_data: Optional[InputProcessor] = None, gene_region: str = "upstream", png_name: str = "upstream_methylation.png", csv_name: str = "upstream_methylation.csv", csv: Optional[str] = None, left_distance: int = 5000, right_distance: int = 5000, window_size: int = 100, hypermethylated: bool = True, gene_hypermethylated_min: int = 20, window_hypermethylated_min: int = 5, min_hypermethylated_windows: int = 5, hypomethylated: bool = True, gene_hypomethylated_max: int = -20, window_hypomethylated_max: int = -5, min_hypomethylated_windows: int = 5, position_count: int = 5, clamp_positive: int = 50, clamp_negative: int = -50, title:str = None, gtf_file: str= "gencode.v41.chr_patch_hapl_scaff.annotation.gtf", position_or_window: str = "auto") -> None:
        """Generate a graph of upstream gene methylation.
        
        .. note::
            Required columns for ``position_data``:
                - ``["chromosome", "position_start", "diff"]``
            Required columns for ``gene_data``:
                - If ``"intron"`` is in ``gene_regions``: ``["intron", "intron_diff"]``
                - If ``"exon"`` is in ``gene_regions``: ``["exon", "exon_diff"]``
                - If ``"upstream"`` is in ``gene_regions``: ``["upstream", "upstream_diff"]``

        :param position_data: Position data. Not necessary if the pipeline is in use, defaults to None
        :type position_data: InputProcessor, optional
        :param gene_data: Gene data. Not necessary if the pipeline is in use, defaults to None
        :type gene_data: InputProcessor, optional
        :param gene_region: Gene region, defaults to "upstream"
        :type gene_region: str, optional
        :param png_name: Output PNG file name, defaults to "upstream_methylation.png"
        :type png_name: str, optional
        :param csv_name: Output CSV file name, defaults to "upstream_methylation.csv"
        :type csv_name: str, optional
        :param csv: Input CSV file, defaults to None
        :type csv: str, optional
        :param left_distance: Left distance, defaults to 5000
        :type left_distance: int, optional
        :param right_distance: Right distance, defaults to 5000
        :type right_distance: int, optional
        :param window_size: Window size, defaults to 100
        :type window_size: int, optional
        :param hypermethylated: Include hypermethylated regions, defaults to True
        :type hypermethylated: bool, optional
        :param gene_hypermethylated_min: Minimum hypermethylation for genes, defaults to 20
        :type gene_hypermethylated_min: int, optional
        :param window_hypermethylated_min: Minimum hypermethylation for windows, defaults to 5
        :type window_hypermethylated_min: int, optional
        :param min_hypermethylated_windows: Minimum hypermethylated windows, defaults to 5
        :type min_hypermethylated_windows: int, optional
        :param hypomethylated: Include hypomethylated regions, defaults to True
        :type hypomethylated: bool, optional
        :param gene_hypomethylated_max: Maximum hypomethylation for genes, defaults to -20
        :type gene_hypomethylated_max: int, optional
        :param window_hypomethylated_max: Maximum hypomethylation for windows, defaults to -5
        :type window_hypomethylated_max: int, optional
        :param min_hypomethylated_windows: Minimum hypomethylated windows, defaults to 5
        :type min_hypomethylated_windows: int, optional
        :param position_count: Position count, defaults to 5
        :type position_count: int, optional
        :param clamp_positive: Positive clamp value, defaults to 50
        :type clamp_positive: int, optional
        :param clamp_negative: Negative clamp value, defaults to -50
        :type clamp_negative: int, optional
        :param title: Plot title, defaults to None for a generic title
        :type title: str, optional
        :param gtf_file: GTF file, defaults to "gencode.v41.chr_patch_hapl_scaff.annotation.gtf"
        :type gtf_file: str, optional
        :param position_or_window: Position or window, options are ["auto", "position", "window"], defaults to "auto"
        :type position_or_window: str, optional
        """
        assert (not self.pipeline and gene_data is not None and position_data is not None) or (not self.pipeline and csv is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert position_or_window in ["auto", "position", "window"], "Invalid parameter for position_or_window. Options are: [""auto"", ""position"", ""window""]"

        parameters = locals().copy()

        if gene_data is not None:
            gene_data = gene_data.copy()
            gene_data.process()
            gene_data = gene_data.data_container
        elif csv is None:
            if position_or_window == "auto":
                if self.saved_results.get((self.map_positions_to_genes, "gene")) is not None:
                    gene_data = self.saved_results[(self.map_positions_to_genes, "gene")]
                elif self.saved_results.get((self.map_windows_to_genes, "gene")) is not None:
                    gene_data = self.saved_results[(self.map_windows_to_genes, "gene")]
                else:
                    # TODO invalid
                    pass
            elif position_or_window == "position":
                gene_data = self.saved_results[(self.map_positions_to_genes, "gene")]
            else:
                gene_data = self.saved_results[(self.map_windows_to_genes, "gene")]

        if position_data is not None:
            position_data = position_data.copy()
            position_data.process()
            position_data = position_data.data_container
        elif csv is None:
            position_data = self.saved_results[self.merge_tables]
        
        parameters.pop("position_or_window")
        parameters = self.__prepare_parameters(parameters, gene_data=gene_data, position_data=position_data)

        self.plots.graph_upstream_gene_methylation(**parameters)

    GRAPH_AVERAGE_UPSTREAM_GENE_METHYLATION_REQUIRED_COLUMNS = {
        "position_data": ["chromosome", "position_start", "diff"],
        "gene_data": ["intron", "intron_diff", "exon", "exon_diff", "upstream", "upstream_diff"]
    }
    def graph_average_upstream_gene_methylation(self, position_data: Optional[InputProcessor] = None, gene_data: Optional[InputProcessor] = None, gene_region: str = "upstream", png_name="average_upstream_methylation.png", csv_name: str = "average_upstream_methylation.csv", csv: Optional[str] = None, left_distance: int = 5000, right_distance: int = 5000, window_size: int = 100, hypermethylated: bool = True, gene_hypermethylated_min: int = 20, window_hypermethylated_min: int = 5, min_hypermethylated_windows: int = 5, hypomethylated: bool = True, gene_hypomethylated_max: int = -20, window_hypomethylated_max: int = -5, min_hypomethylated_windows: int = 5, position_count: int = 5, clamp_positive: int = 50, clamp_negative: int = -50, title:str = None, y_label:str=None, gtf_file: str= "gencode.v41.chr_patch_hapl_scaff.annotation.gtf", position_or_window="auto") -> None:   
        """Generate a graph of average upstream gene methylation.

        .. note::
            Required columns for ``position_data``:
                - ``["chromosome", "position_start", "diff"]``
            Required columns for ``gene_data``:
                - If ``"intron"`` is in ``gene_regions``: ``["intron", "intron_diff"]``
                - If ``"exon"`` is in ``gene_regions``: ``["exon", "exon_diff"]``
                - If ``"upstream"`` is in ``gene_regions``: ``["upstream", "upstream_diff"]``

        :param position_data: Position data. Not necessary if the pipeline is in use, defaults to None
        :type position_data: InputProcessor, optional
        :param gene_data: Gene data. Not necessary if the pipeline is in use, defaults to None
        :type gene_data: InputProcessor, optional
        :param gene_region: Gene region, defaults to "upstream"
        :type gene_region: str, optional
        :param png_name: Output PNG file name, defaults to "upstream_methylation.png"
        :type png_name: str, optional
        :param csv_name: Output CSV file name, defaults to "upstream_methylation.csv"
        :type csv_name: str, optional
        :param csv: Input CSV file, defaults to None
        :type csv: str, optional
        :param left_distance: Left distance, defaults to 5000
        :type left_distance: int, optional
        :param right_distance: Right distance, defaults to 5000
        :type right_distance: int, optional
        :param window_size: Window size, defaults to 100
        :type window_size: int, optional
        :param hypermethylated: Include hypermethylated regions, defaults to True
        :type hypermethylated: bool, optional
        :param gene_hypermethylated_min: Minimum hypermethylation for genes, defaults to 20
        :type gene_hypermethylated_min: int, optional
        :param window_hypermethylated_min: Minimum hypermethylation for windows, defaults to 5
        :type window_hypermethylated_min: int, optional
        :param min_hypermethylated_windows: Minimum hypermethylated windows, defaults to 5
        :type min_hypermethylated_windows: int, optional
        :param hypomethylated: Include hypomethylated regions, defaults to True
        :type hypomethylated: bool, optional
        :param gene_hypomethylated_max: Maximum hypomethylation for genes, defaults to -20
        :type gene_hypomethylated_max: int, optional
        :param window_hypomethylated_max: Maximum hypomethylation for windows, defaults to -5
        :type window_hypomethylated_max: int, optional
        :param min_hypomethylated_windows: Minimum hypomethylated windows, defaults to 5
        :type min_hypomethylated_windows: int, optional
        :param position_count: Position count, defaults to 5
        :type position_count: int, optional
        :param clamp_positive: Positive clamp value, defaults to 50
        :type clamp_positive: int, optional
        :param clamp_negative: Negative clamp value, defaults to -50
        :type clamp_negative: int, optional
        :param title: Plot title, defaults to None for a generic title
        :type title: str, optional
        :param y_label: Y-axis label, defaults to None for a generic label
        :type y_label: str, optional
        :param gtf_file: GTF file, defaults to "gencode.v41.chr_patch_hapl_scaff.annotation.gtf"
        :type gtf_file: str, optional
        :param position_or_window: Position or window, defaults to "auto"
        :type position_or_window: str, optional
        """
        assert (not self.pipeline and gene_data is not None and position_data is not None) or (not self.pipeline and csv is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert position_or_window in ["auto", "position", "window"], "Invalid parameter for position_or_window. Options are: [""auto"", ""position"", ""window""]"

        parameters = locals().copy()

        if gene_data is not None:
            gene_data = gene_data.copy()
            gene_data.process()
            gene_data = gene_data.data_container
        elif csv is None:
            if position_or_window == "auto":
                if self.saved_results.get((self.map_positions_to_genes, "gene")) is not None:
                    gene_data = self.saved_results[(self.map_positions_to_genes, "gene")]
                elif self.saved_results.get((self.map_windows_to_genes, "gene")) is not None:
                    gene_data = self.saved_results[(self.map_windows_to_genes, "gene")]
                else:
                    # TODO invalid
                    pass
            elif position_or_window == "position":
                gene_data = self.saved_results[(self.map_positions_to_genes, "gene")]
            else:
                gene_data = self.saved_results[(self.map_windows_to_genes, "gene")]

        if position_data is not None:
            position_data = position_data.copy()
            position_data.process()
            position_data = position_data.data_container
        elif csv is None:
            position_data = self.saved_results[self.merge_tables]
        
        parameters.pop("position_or_window")
        parameters = self.__prepare_parameters(parameters, gene_data=gene_data, position_data=position_data)

        self.plots.graph_average_upstream_gene_methylation(**parameters)

    GRAPH_UPSTREAM_UCSC_REQUIRED_COLUMNS = {
        "position_data": ["chromosome", "position_start", "methylation_percentage*"]
    }
    def graph_upstream_UCSC(self, gene_name: str, position_data: Optional[InputProcessor] = None, name: str="UCSC_graph.bedGraph", before_tss: int = 5000, gtf_file: str = "gencode.v41.chr_patch_hapl_scaff.annotation.gtf") -> None:
        """Generate a UCSC graph of upstream gene methylation.

        .. note::
            Required columns for ``position_data``:
                - ``["chromosome", "position_start", "methylation_percentage*"]``
                - There must be a ``methylation_percentage`` column for each sample to plot.

        :param gene_name: Gene name
        :type gene_name: str
        :param position_data: Position data
        :type position_data: InputProcessor, optional
        :param name: Output BEDGraph file name, defaults to "UCSC_graph.bedGraph"
        :type name: str, optional
        :param before_tss: Distance before transcrption start site (TSS), defaults to 5000
        :type before_tss: int, optional
        :param gtf_file: GTF file, defaults to "gencode.v41.chr_patch_hapl_scaff.annotation.gtf"
        :type gtf_file: str, optional
        """
        assert (not self.pipeline and position_data is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        parameters = locals().copy()

        if position_data is not None:
            position_data = position_data.copy()
            position_data.process()
            position_data = position_data.data_container
        else:
            position_data = self.saved_results[self.merge_tables]
        
        parameters = self.__prepare_parameters(parameters, position_data=position_data)

        self.plots.graph_upstream_UCSC(**parameters)

    GRAPH_FULL_GENE_REQUIRED_COLUMNS = {
        "position_data": ["chromosome", "position_start", "methylation_percentage*"]
    }
    def graph_full_gene(self, gene_name:str, position_data: Optional[InputProcessor] = None, name="gene_methylation_graph.png", before_tss: int = 0, after_tss: Optional[int] = None, bin_size: int = 500, start_marker: bool = True, end_marker: bool = True, deviation_display: bool = True, legend_size:int = 12, title: str = None, x_label:str = None, y_label:str=None, case_name: str = "Case", ctr_name: str = "Control",  gtf_file: str = "gencode.v41.chr_patch_hapl_scaff.annotation.gtf") -> None:
        """Generate a graph of full gene methylation.

        .. note::
            Required columns for ``position_data``:
                - ``["chromosome", "position_start", "methylation_percentage*"]``
                - There must be a ``methylation_percentage`` column for each sample to plot.

        :param gene_name: The name of the gene to graph.
        :type gene_name: str
        :param position_data: Position data. Not necessary if the pipeline is in use, defaults to None
        :type position_data: InputProcessor, optional
        :param name: Output PNG file name, defaults to "gene_methylation_graph.png"
        :type name: str, optional
        :param before_tss: Distance before transcription start site, defaults to 0
        :type before_tss: int, optional
        :param after_tss: Distance after transcription start site, defaults to None
        :type after_tss: int, optional
        :param bin_size: Bin size, defaults to 500
        :type bin_size: int, optional
        :param start_marker: Include start marker, defaults to True
        :type start_marker: bool, optional
        :param end_marker: Include end marker, defaults to True
        :type end_marker: bool, optional
        :param deviation_display: Display standard deviation regions, defaults to True
        :type deviation_display: bool, optional
        :param legend_size: Legend size, defaults to 12
        :type legend_size: int, optional
        :param title: Plot title, defaults to None for a generic title
        :type title: str, optional
        :param x_label: X-axis label, defaults to None for a generic label
        :type x_label: str, optional
        :param y_label: Y-axis label, defaults to None for a generic label
        :type y_label: str, optional
        :param case_name: Case name in the legend, defaults to "Case"
        :type case_name: str, optional
        :param ctr_name: Control name in the legend, defaults to "Control"
        :type ctr_name: str, optional
        :param gtf_file: GTF file, defaults to "gencode.v41.chr_patch_hapl_scaff.annotation.gtf"
        :type gtf_file: str, optional
        """
        assert (not self.pipeline and position_data is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        parameters = locals().copy()

        if position_data is not None:
            position_data = position_data.copy()
            position_data.process()
            position_data = position_data.data_container
        else:
            position_data = self.saved_results[self.merge_tables]
        
        parameters = self.__prepare_parameters(parameters, position_data=position_data)

        self.plots.graph_full_gene(**parameters)

    PIE_CHART_REQUIRED_COLUMNS = {
        "gene_data": ["intron", "intron_diff", "exon", "exon_diff", "upstream", "upstream_diff"],
        "ccre_data": ["CCRE", "CCRE_diff"]
    }
    def pie_chart(self, gene_data: Optional[InputProcessor] = None, CCRE_data: Optional[InputProcessor] = None, regions: list[str]|str = ["intron", "exon", "CCRE", "upstream"], name: str = "pie_chart.png", hypermethylated:bool = True, hypomethylated:bool = True, hypermehylated_min:float = 20, hypomethylated_max:float = -20, hypermethylated_title:str = None, hypomethylated_title:str = None, title:str = None, position_or_window: str = "auto") -> None:
        """Generate a pie chart of gene regions.

        .. note::
            Required columns for ``gene_data``:
                - If ``"intron"`` is in ``gene_regions``: ``["intron"]``
                - If ``"exon"`` is in ``gene_regions``: ``["exon"]``
                - If ``"upstream"`` is in ``gene_regions``: ``["upstream"]``
            Required columns for ``CCRE_data``:
                - If ``"CCRE"`` is in ``gene_regions``: ``["CCRE"]``

        :param gene_data: Gene data. Not necessary if the pipeline is in use, defaults to None
        :type gene_data: InputProcessor, optional
        :param CCRE_data: CCRE data. Not necessary if the pipeline is in use, defaults to None
        :type CCRE_data: InputProcessor, optional
        :param regions: Gene regions, defaults to ["intron", "exon", "CCRE", "upstream"]
        :type regions: list[str] | str, optional
        :param name: Output PNG file name, defaults to "pie_chart.png"
        :type name: str, optional
        :param hypermethylated: Include hypermethylated regions, defaults to True
        :type hypermethylated: bool, optional
        :param hypomethylated: Include hypomethylated regions, defaults to True
        :type hypomethylated: bool, optional
        :param hypermehylated_min: Minimum hypermethylation, defaults to 20
        :type hypermehylated_min: float, optional
        :param hypomethylated_max: Maximum hypomethylation, defaults to -20
        :type hypomethylated_max: float, optional
        :param hypermethylated_title: Hypermethylated title, defaults to None for a generic title
        :type hypermethylated_title: str, optional
        :param hypomethylated_title: Hypomethylated title, defaults to None for a generic title
        :type hypomethylated_title: str, optional
        :param title: Plot title, defaults to None for a generic title
        :type title: str, optional
        :param position_or_window: Position or window, options are ["auto", "position", "window"], defaults to "auto"
        :type position_or_window: str, optional
        """
        assert (not self.pipeline and (gene_data is not None or CCRE_data is not None)) or (self.pipeline), "If the pipeline isn't in use, data must be provided."
        assert position_or_window in ["auto", "position", "window"], "Invalid parameter for position_or_window. Options are: [""auto"", ""position"", ""window""]"


        parameters = locals().copy()

        if gene_data is not None:
            gene_data = gene_data.copy()
            gene_data.process()
            gene_data = gene_data.data_container
        else:
            if position_or_window == "auto":
                if self.saved_results.get((self.map_positions_to_genes, "gene")) is not None:
                    gene_data = self.saved_results[(self.map_positions_to_genes, "gene")]
                elif self.saved_results.get((self.map_windows_to_genes, "gene")) is not None:
                    gene_data = self.saved_results[(self.map_windows_to_genes, "gene")]
                else:
                    # TODO invalid
                    pass
            elif position_or_window == "position":
                gene_data = self.saved_results[(self.map_positions_to_genes, "gene")]
            else:
                gene_data = self.saved_results[(self.map_windows_to_genes, "gene")]

        if CCRE_data is not None:
            CCRE_data = CCRE_data.copy()
            CCRE_data.process()
            CCRE_data = CCRE_data.data_container
        else:
            if position_or_window == "auto":
                if self.saved_results.get((self.map_positions_to_genes, "CCRE")) is not None:
                    CCRE_data = self.saved_results[(self.map_positions_to_genes, "CCRE")]
                elif self.saved_results.get((self.map_windows_to_genes, "CCRE")) is not None:
                    CCRE_data = self.saved_results[(self.map_windows_to_genes, "CCRE")]
                else:
                    # TODO invalid
                    pass
            elif position_or_window == "position":
                CCRE_data = self.saved_results[(self.map_positions_to_genes, "CCRE")]
            else:
                CCRE_data = self.saved_results[(self.map_windows_to_genes, "CCRE")]
        
        parameters.pop("position_or_window")
        parameters = self.__prepare_parameters(parameters, gene_data=gene_data, CCRE_data=CCRE_data)

        self.plots.pie_chart(**parameters)

    ALL_ANALYSIS_REQUIRED_COLUMNS = {
        "case_data": ["chromosome", "position_start", "coverage", "methylation_percentage", "positive_methylation_count", "negative_methylation_count"],
        "ctr_data": ["chromosome", "position_start", "coverage", "methylation_percentage", "positive_methylation_count", "negative_methylation_count"]
    }

    def all_analysis(self, case_data: InputProcessor, ctr_data: InputProcessor, window_based=True, prefix=".", min_cov = 10, cov_percentile = 1.0, min_samp_ctr = 2, min_samp_case = 2) -> pd.DataFrame:
        """Run all analysis methods.

        .. note::
            ``case_data`` and ``ctr_data`` must be from a list of samples, where each contains one of the following column formats:
                - ``["chromosome", "position_start", "coverage", "methylation_percentage"]``
                - ``["chromosome", "position_start", "positive_methylation_count", "negative_methylation_count"]``

        .. note::
            if ``window_based`` is True, the following methods will be run:
                - ``merge_tables``
                - ``window_based``
                - ``generate_q_values``
                - ``filters``
                - ``map_win_2_pos``
            if ``window_based`` is False, the following methods will be run:
                - ``merge_tables``
                - ``position_based``
                - ``generate_q_values``
                - ``filters``

        :param case_data: Case data (list of files)
        :type case_data: InputProcessor
        :param ctr_data: Control data (list of files)
        :type ctr_data: InputProcessor
        :param window_based: Window-based analysis, defaults to True
        :type window_based: bool, optional
        :param prefix: Output file directory, defaults to "."
        :type prefix: str, optional
        :param min_cov: Minimum coverage, defaults to 10
        :type min_cov: int, optional
        :param cov_percentile: Coverage percentile, defaults to 1.0
        :type cov_percentile: float, optional
        :param min_samp_ctr: Minimum control samples, defaults to 2
        :type min_samp_ctr: int, optional
        :param min_samp_case: Minimum case samples, defaults to 2
        :type min_samp_case: int, optional
        :return: Final significant positions
        :rtype: pd.DataFrame
        """
        self.pipeline = False # True may take a lot of memory
        print("merging")
        merged = self.merge_tables(case_data, ctr_data, min_cov=min_cov, cov_percentile=cov_percentile, min_samp_ctr=min_samp_ctr, min_samp_case=min_samp_case)
        del case_data, ctr_data
        merged.to_csv(f"{prefix}/merged.csv", index=False)
        if window_based:
            print("windows")
            res = self.window_based(InputProcessor(merged))
            res.to_csv(f"{prefix}/windows.csv", index=False)
            print("p-val")
            res = self.generate_q_values(InputProcessor(res))
            res.to_csv(f"{prefix}/q_values.csv", index=False)
            print("filter")
            res = self.filters(InputProcessor(res))
            res.to_csv(f"{prefix}/significant_windows.csv", index=False)
            print("mapping")
            res = self.map_win_2_pos(InputProcessor(res), InputProcessor(merged))
            res.to_csv(f"{prefix}/mapped_positions.csv", index=False)
        else:
            res = self.position_based(InputProcessor(merged))
            del merged
            res.to_csv(f"{prefix}/positions.csv", index=False)
            res = self.generate_q_values(InputProcessor(res))
            res.to_csv(f"{prefix}/q_values.csv", index=False)
            res = self.filters(InputProcessor(res))
            res.to_csv(f"{prefix}/significant_positions.csv", index=False)
        return res
    
    ALL_PLOTS_REQUIRED_COLUMNS = {
        "data": ["chromosome", "position_start", "region_start", "q-value", "diff", "methylation_percentage*", "coverage"],
        "gene_data": ["intron", "intron_diff", "exon", "exon_diff", "upstream", "upstream_diff"],
        "ccre_data": ["CCRE", "CCRE_diff"]
    }
    def all_plots(self, data: InputProcessor, gene_data: InputProcessor, ccre_data: InputProcessor) -> None:
        """Run all plot methods."""
        self.pipeline = False
        self.volcano_plot(data)
        self.manhattan_plot(data)
        self.graph_gene_regions(gene_data, ccre_data)
        self.graph_upstream_gene_methylation(data, gene_data)
        self.graph_average_upstream_gene_methylation(data, gene_data)
        self.pie_chart(gene_data, ccre_data)
    



def parse_arguments():
    parser = argparse.ArgumentParser(description="CLI for DiffMethTools")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    obj = DiffMethTools(pipeline=False)
    # Dynamically add subcommands for each method in DiffMethTools
    for name, method in inspect.getmembers(obj, predicate=inspect.ismethod):
        if name.startswith("_"):  # Skip private methods
            continue

        names_argument_added = []

        method_parser = subparsers.add_parser(name, help=f"Run the {name} method")
        sig = inspect.signature(method)
        for param_name, param in sig.parameters.items():
            if param_name == "self":
                continue
            arg_type = param.annotation if param.annotation != inspect.Parameter.empty else str
            if type(arg_type) == UnionType:
                arg_type = Union

            if arg_type != InputProcessor and arg_type != Optional[InputProcessor]:
                default = param.default if param.default != inspect.Parameter.empty else None
                if param.kind != param.POSITIONAL_ONLY:
                    if get_origin(arg_type) is Optional:
                        arg_type = param.annotation.__args__[0]
                    method_parser.add_argument(f"--{param_name}", type=arg_type, default=default, help=f"(default: {default})")
                else:
                    method_parser.add_argument(f"--{param_name}", type=arg_type, required=True, help="(required)")
            else:
                method_parser.add_argument(f"--{param_name}_file", type=str, nargs="+", required=True, help="(required, can be a list of files)")
                method_parser.add_argument(f"--{param_name}_has_header", action="store_true", required=False, help="(use this flag if there is a header.)")
                method_parser.add_argument(f"--{param_name}_separator", type=str, required=False, default=',', help="(default: ,)")
                required_columns = getattr(obj, f"{name.upper()}_REQUIRED_COLUMNS")
                for data_name, columns in required_columns.items():
                    if param_name == data_name:
                        for column in columns:
                            if "methylation_percentage" in column or "methylation_count" in column or "coverage" in column:
                                # remove * from column if exists
                                column = column[:-1] if column[-1] == "*" else column
                                # allow list input
                                method_parser.add_argument(f"--{param_name}_{column}_column_index", type=int, nargs="+", required=False, help="(default: None, assumes all column are named as expected)")
                                if param_name not in names_argument_added:
                                    method_parser.add_argument(f"--{param_name}_column_names", type=str, nargs="+", required=False, help="(default: None, assumes all column are named as expected)")
                                    if "methylation_count" not in column and "ctr" not in param_name and "case" not in param_name:
                                        method_parser.add_argument(f"--{param_name}_column_control_case", type=str, nargs="+", required=False, help="(default: None, assumes column names are not provided)")
                                    names_argument_added.extend([param_name])
                            else:
                                method_parser.add_argument(f"--{param_name}_{column}_column_index", type=int, required=False, help="(default: None, assumes all columns are named as expected)")

        # add output parameter if return type is not none
        return_annotation = sig.return_annotation
        if return_annotation != inspect.Parameter.empty and return_annotation is not None and return_annotation != type(None):
            # Determine if it's a DataFrame or collection of DataFrames
            is_dataframe = return_annotation == pd.DataFrame
            is_collection = False
            
            if get_origin(return_annotation) in (tuple, list):
                args = get_args(return_annotation)
                is_collection = any(arg == pd.DataFrame for arg in args)
            
            if is_dataframe:
                default_output = f"{name}.csv"

                method_parser.add_argument("--output", type=str, default=default_output, 
                                        help=f"Output filename (default: {default_output})")
            elif is_collection:
                method_parser.add_argument("--output_prefix", type=str, default=name, 
                                        help=f"Prefix for output filenames (default: {name})")

    return parser

def main():
    parser = parse_arguments()
    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        return

    # Initialize DiffMethTools instance
    tool = DiffMethTools(pipeline=False)

    # Get the selected method
    method = getattr(tool, args.command)

    # Prepare arguments for the method
    sig = inspect.signature(method)
    method_args = {}
    output = method.__name__ + ".csv"
    for param_name, param in sig.parameters.items():
        if param_name == "self":
            continue
        arg_type = param.annotation if param.annotation != inspect.Parameter.empty else str

        if param_name == "output" or param_name == "output_prefix":
            output = getattr(args, param_name, "out.csv" if param_name == "output" else "out")
            continue
        
        if arg_type == InputProcessor or arg_type == Optional[InputProcessor]:
            data_indicies = {k[k.find(param_name)+len(param_name)+1:k.rfind("_column_index")] : v for k, v in vars(args).items() if "index" in k and param_name in k and v is not None}
            
            data_indicies = {k if k != "coverage" else "coverage_KEY": v for k, v in data_indicies.items()}
            # remove trailing * if exists
            # data_indicies = {k[:-1] if k[-1] == "*" else k: v for k, v in data_indicies.items()}
            file_name = getattr(args, param_name+"_file", None)
            has_header = getattr(args, param_name+"_has_header", True)
            separator = getattr(args, param_name+"_separator", ",")

            if separator == "t":
                separator = "\t"
                
            names = getattr(args, param_name+"_column_names", None)

            if "ctr" in param_name:
                ctr_case = "ctr"
            elif "case" in param_name:
                ctr_case = "case"
            else:
                ctr_case = getattr(args, param_name+"_column_control_case", None)
            
            format = FormatDefinition(column_mapping=data_indicies, sep=separator)

            # Convert CSV file paths to pandas DataFrames
            if len(file_name) == 1:
                file_name = file_name[0]

            if ctr_case is not None:
                if names is not None:
                    value = InputProcessor(file_name, format, has_header=has_header, names=names, ctr_case=ctr_case)
                else:
                    value = InputProcessor(file_name, format, has_header=has_header, ctr_case=ctr_case)
            else:
                if names is not None:
                    value = InputProcessor(file_name, format, has_header=has_header, names=names)
                else:
                    value = InputProcessor(file_name, format, has_header=has_header)
        else:
            value = getattr(args, param_name, None)
        method_args[param_name] = value

    # Call the method and print the result
    result = method(**method_args)
    # TODO send to csv if dataframe output. Remember there can be multiple dataframes in output. Make default name <function>.csv
    print(output)
    if type(result) == pd.DataFrame:
        result.to_csv(output, index=False)
    elif type(result) == list:
        for i, df in enumerate(result):
            df.to_csv(f"{output}_{i}.csv", index=False)

if __name__ == "__main__":
    main()

