from typing import Callable, Optional, Type, get_args, get_origin, Union
from types import UnionType
from lib import InputProcessor, Plots, FormatDefinition, Analysis#################Analysis
from scipy.stats import mannwhitneyu
from functools import wraps
import inspect
import argparse
import pandas as pd
import os
import yaml
import inspect
from pathlib import Path

class DiffMethylTools():
    def __init__(self, pipeline=True, results_path=None):
        self.pipeline = pipeline
        if results_path is None:
            results_path = os.getcwd()
        if results_path is not None and results_path[-1] == "/":
            results_path = results_path[:-1]
        self.results_path = results_path
        self.obj = Analysis()
        self.plots = Plots()
        self.saved_results = {}

    def __prepare_parameters(self, parameters, **kwargs):
        parameters.pop("self")
        if "rerun" in parameters.keys():
            parameters.pop("rerun")
        
        for key in parameters:
            if key in kwargs.keys():
                parameters[key] = kwargs[key]
        
        return parameters
    
    def analysis_function(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            if os.path.exists(self.results_path + "/data/state.yaml"):
                with open(self.results_path + "/data/state.yaml", "r") as file:
                    state = yaml.safe_load(file)

                
                func_state = state.get(func.__name__, False)
                # if function is "True" in state.yml and rerun is "False", 
                if func_state and not kwargs.get("rerun", False):
                    
                    is_gene = "to_genes" in func.__name__
                    genes = False
                    CCRE = False

                    if is_gene:
                        genes = os.path.exists(self.results_path + "/data/" + func.__name__ + "_genes.csv")
                        CCRE = os.path.exists(self.results_path + "/data/" + func.__name__ + "_CCRE.csv")
                
                    if genes or CCRE or os.path.exists(self.results_path + "/data/" + func.__name__ + ".csv"):     
                        print(f"Skipping {func.__name__} as it has already been run.")
                        if not genes and not CCRE:
                            res = InputProcessor(self.results_path + "/data/" + func.__name__ + ".csv", format=FormatDefinition(sep=","))
                            res.process()
                            result = res.data_container
                        else:
                            res = []
                            if genes:
                                res.append(InputProcessor(self.results_path + "/data/" + func.__name__ + "_genes.csv", format=FormatDefinition(sep=",")))
                                res[-1].process()
                                res[-1] = res[-1].data_container
                            if CCRE:
                                res.append(InputProcessor(self.results_path + "/data/" + func.__name__ + "_CCRE.csv", format=FormatDefinition(sep=",")))
                                res[-1].process()
                                res[-1] = res[-1].data_container
                            print(res)
                            result = tuple(res)
                        return result
                        
                        # turn chromstart chromend into ints
                        # result["chr"] = result["chr"].astype(str)
                        # result["chromStart"] = result["chromStart"].astype(int)
                        # result["chromEnd"] = result["chromEnd"].astype(int)
                        return result
                    else:
                        state[func.__name__] = False
                        with open(self.results_path + "/data/state.yaml", "w") as file:
                            yaml.dump(state, file)

            # make folders
            if self.results_path is not None:  
                os.makedirs(self.results_path, exist_ok=True)
                os.makedirs(self.results_path + "/plots", exist_ok=True)
                os.makedirs(self.results_path + "/data", exist_ok=True)

                
            png_name = kwargs.get("name", None)
            csv_name = kwargs.get("csv_name", None)

            if png_name is not None:
                if "/" not in png_name:
                    kwargs["name"] = self.results_path + "/plots/" + png_name

            if csv_name is not None:
                if "/" not in csv_name:
                    kwargs["csv_name"] = self.results_path + "/plots/" + csv_name

            result = func(self, *args, **kwargs)

            if self.results_path is not None:
                
                try:
                    with open(self.results_path + "/data/state.yaml", "r") as file:
                        state = yaml.safe_load(file)
                except FileNotFoundError:
                    state = {}

                state[func.__name__] = True

                with open(self.results_path + "/data/state.yaml", "w") as file:
                    yaml.dump(state, file)

                if result is not None:
                    if isinstance(result, tuple):
                        result = list(result)
                    elif isinstance(result, pd.DataFrame):
                        result = [result]
                    
                    result_len = len(result)
                    for i, dataframe in enumerate(result):
                        if result_len > 1:
                            if func.__name__ == "map_positions_to_genes" or func.__name__ == "map_windows_to_genes":
                                if i == 0:
                                    dataframe.to_csv(self.results_path + "/data/" + func.__name__ + "_genes.csv", index=False)
                                else:
                                    dataframe.to_csv(self.results_path + "/data/" + func.__name__ + "_CCRE.csv", index=False)
                            else:
                                dataframe.to_csv(self.results_path + "/data/" + func.__name__ + "_" + str(i) + ".csv", index=False)
                        else:
                            if func.__name__ == "map_positions_to_genes" or func.__name__ == "map_windows_to_genes":
                                for j, df in enumerate(dataframe):
                                    if df is not None:
                                        if isinstance(df, pd.DataFrame):
                                            df.to_csv(self.results_path + "/data/" + func.__name__ + "_" + ("genes" if j == 0 else "CCRE") + ".csv", index=False)
                                        else:
                                            raise ValueError(f"Invalid dataframe type: {type(df)}")
                            else:
                                dataframe.to_csv(self.results_path + "/data/" + func.__name__ + ".csv", index=False)
                
            return result
        return wrapper

    MERGE_TABLES_REQUIRED_COLUMNS = {
        "case_data": ["chromosome", "position_start", "coverage", "methylation_percentage", "positive_methylation_count", "negative_methylation_counta", "strand"],
        "ctr_data": ["chromosome", "position_start", "coverage", "methylation_percentage", "positive_methylation_count", "negative_methylation_count", "strand"]
    }
    @analysis_function
    def merge_tables(self, case_data: InputProcessor, ctr_data: InputProcessor, min_cov_individual = 10, min_cov_group = 15, filter_samples_ratio=0.6, meth_group_threshold=0.2, cov_percentile = 100.0, min_samp_ctr = 2, min_samp_case = 2, rerun=False, small_mean = 1) -> pd.DataFrame:
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
        :param min_cov_individual: Minimum coverage filter (individual), defaults to 10
        :type min_cov_individual: int, optional
        :type min_cov: int, optional
        :param min_cov_group: Minimum coverage filter (group), defaults to 15
        :param filter_samples_ratio: Minimum sample ratio filter. Used with min_cov_group, defaults to 0.6
        :type filter_samples_ratio: float, optional
        :param meth_group_threshold: Methylation group threshold. Used with min_cov_group, defaults to 0.2
        :type meth_group_threshold: float, optional
        :type min_cov_group: int, optional
        :param cov_percentile: Maximum coverage filter (percentile of sample coverage). Ranges from 0.0-100.0, defaults to 100.0
        :type cov_percentile: float, optional
        :param min_samp_ctr: Minimum samples in control, defaults to 2
        :type min_samp_ctr: int, optional
        :param min_samp_case: Minimum samples in case, defaults to 2
        :type min_samp_case: int, optional
        :param rerun: Rerun the analysis. If False, load previous output. Defaults to False.
        :type rerun: bool, optional
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
            self.saved_results[self.merge_tables.__name__] = res

        return res
    POSITION_BASED_REQUIRED_COLUMNS = {
        "data": ["chromosome", "position_start", "methylation_percentage*"]
    }
    @analysis_function
    def position_based(self, data: Optional[InputProcessor] = None , method="limma", features=None, test_factor="Group", processes=12, model="eBayes", min_std=0.1, fill_na:bool=True, rerun=False) -> pd.DataFrame:
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
        :type fill_na: bool, optional
        :param rerun: Rerun the analysis. If False, load previous output. Defaults to False.
        :type rerun: bool, optional
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
            print(self.saved_results)
            data = self.saved_results[self.merge_tables.__name__]
            
        parameters = self.__prepare_parameters(parameters, data=data)

        res = self.obj.position_based(**parameters)

        if self.pipeline:
            self.saved_results[self.position_based.__name__] = res

        return res

    MAP_WIN_2_POS_REQUIRED_COLUMNS = {
        "window_data": ["chromosome", "region_start", "region_end"],
        "position_data": ["chromosome", "position_start", "avg_case", "avg_ctr"]
    }
    @analysis_function
    def map_win_2_pos(self, window_data: Optional[InputProcessor] = None, position_data: Optional[InputProcessor] = None, processes=12, sub_window_size = 100, sub_window_step = 100, sub_window_min_diff=0, pipeline_window_result: str = "auto", rerun=False) -> pd.DataFrame:
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
        :param pipeline_window_result: The function results to use as input if DiffMethylTools is pipelined and no data is provided. Options are ``["auto", "filters", "generate_q_values", "window_based"]``, defaults to "auto"
        :type pipeline_window_result: str, optional
        :param rerun: Rerun the analysis. If False, load previous output. Defaults to False.
        :type rerun: bool, optional
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
                if self.saved_results.get(self.filters.__name__) is not None:
                    window_data = self.saved_results[self.window_based.__name__]
                elif self.saved_results.get(self.generate_q_values.__name__) is not None:
                    window_data = self.saved_results[self.generate_q_values.__name__]
                elif self.saved_results.get(self.window_based.__name__) is not None:
                    window_data = self.saved_results[self.window_based.__name__]
            elif pipeline_window_result == "filters":
                window_data = self.saved_results[self.filters.__name__]
            elif pipeline_window_result == "generate_q_values":
                window_data = self.saved_results[self.generate_q_values.__name__]
            elif pipeline_window_result == "window_based":
                window_data = self.saved_results[self.window_based.__name__]
        if position_data is not None:
            position_data = position_data.copy()
            position_data.process()
            position_data = position_data.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            position_data = self.saved_results[self.merge_tables.__name__]
        
        parameters.pop("pipeline_window_result")
        parameters = self.__prepare_parameters(parameters, window_data=window_data, position_data=position_data)

        res = self.obj.map_win_2_pos(**parameters)

        if self.pipeline:
            self.saved_results[self.map_win_2_pos.__name__] = res

        return res

    FILTERS_REQUIRED_COLUMNS = {
        "data": ["q-value", "diff"]
    }
    @analysis_function
    def filters(self, data: Optional[InputProcessor] = None, max_q_value=0.05, abs_min_diff=0.10, position_or_window: str = "auto", rerun=False) -> pd.DataFrame:
        """Filter data by q-value and minimum difference.

        .. note::
            Required columns for ``data``:
                - ``["q-value", "diff"]``

        :param data: Input data. Not necessary if the pipeline is in use, defaults to None
        :type data: InputProcessor, optional
        :param max_q_value: Maximum q-value filter, defaults to 0.05
        :type max_q_value: float, optional
        :param abs_min_diff: Absolute minimum difference filter, defaults to 0.25
        :type abs_min_diff: int, optional
        :param position_or_window: The position-based or window-based results to use as input if DiffMethylTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto"
        :type position_or_window: str, optional
        :param rerun: Rerun the analysis. If False, load previous output. Defaults to False.
        :type rerun: bool, optional
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
                if self.saved_results.get((self.generate_q_values.__name__, "window")) is not None:
                    data = self.saved_results[(self.generate_q_values.__name__, "window")]
                    position_or_window = "window"
                elif self.saved_results.get((self.generate_q_values.__name__, "position")) is not None:
                    data = self.saved_results[(self.generate_q_values.__name__, "position")]
                    position_or_window = "position"
                else:
                    # TODO invalid
                    pass
        
        
        parameters.pop("position_or_window")
        parameters = self.__prepare_parameters(parameters, data=data)

        res = self.obj.filters(**parameters)

        if self.pipeline:
            self.saved_results[(self.filters.__name__, position_or_window)] = res

        return res

    GENERATE_DMR_REQUIRED_COLUMNS = {
        "significant_position_data": ["chromosome", "position_start", "diff"],
        "position_data": ["chromosome", "position_start", "diff"]
    }
    @analysis_function
    def generate_DMR(self, significant_position_data: Optional[InputProcessor] = None, position_data: Optional[InputProcessor] = None, min_pos=3, neural_change_limit=7.5, neurl_perc=30, opposite_perc=10, significant_position_pipeline: str = "auto", rerun=False) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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
        :param neurl_perc: Neutral percentage, defaults to 30
        :type neurl_perc: int, optional
        :param opposite_perc: Opposite percentage, defaults to 10
        :type opposite_perc: int, optional
        :param significant_position_pipeline: The significant position-based or window-based results to use as input if DiffMethylTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto"
        :type significant_position_pipeline: str, optional
        :param rerun: Rerun the analysis. If False, load previous output. Defaults to False.
        :type rerun: bool, optional
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
                if self.saved_results.get((self.filters.__name__, "position")) is not None:
                    significant_position_data = self.saved_results[(self.filters.__name__, "position")]
                elif self.saved_results.get(self.map_win_2_pos.__name__) is not None:
                    significant_position_data = self.saved_results[self.map_win_2_pos.__name__]
                else:
                    # TODO invalid
                    pass
            elif significant_position_pipeline == "position":
                significant_position_data = self.saved_results[(self.filters.__name__, "position")]
            else:
                significant_position_data = self.saved_results[self.map_win_2_pos.__name__]

        if position_data is not None:
            position_data = position_data.copy()
            position_data.process()
            position_data = position_data.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            position_data = self.saved_results[self.merge_tables.__name__]
        
        
        parameters.pop("significant_position_pipeline")
        parameters = self.__prepare_parameters(parameters, significant_position_data=significant_position_data, position_data=position_data)

        res = self.obj.generate_DMR(**parameters)

        if self.pipeline:
            self.saved_results[(self.generate_DMR.__name__, "cluster_df")] = res[0]
            self.saved_results[(self.generate_DMR.__name__, "unclustered_dms_df")] = res[1]
            self.saved_results[(self.generate_DMR.__name__, "clustered_dms_df")] = res[2]
        return res
    
    MAP_POSITIONS_TO_GENES_REQUIRED_COLUMNS = {
        "positions": ["chromosome", "position_start", "diff"]
    }
    @analysis_function
    def map_positions_to_genes(self, positions: Optional[InputProcessor] = None, ref_folder:str = None, gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], min_pos_diff=0, gtf_file="gencode.chr_patch_hapl_scaff.annotation.gtf", bed_file="CpG_gencode_annotation.bed", pipeline_input_source = "auto", rerun=False) -> tuple[pd.DataFrame, pd.DataFrame]:
    #def map_positions_to_genes(self, positions: Optional[InputProcessor] = None, gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], min_pos_diff=0, bed_file="CpG_gencode_annotation.bed", gtf_file="outfile_w_hm450.bed", pipeline_input_source = "auto", rerun=False) -> tuple[pd.DataFrame, pd.DataFrame]:
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
        :param bed_file: BED annotation file with unflexible input format, defaults to "CpG_gencode_annotation.bed"
        :type bed_file: str, optional
        :param gtf_file: GTF annotation file with unflexible input format, defaults to "gencode.v42.chr_patch_hapl_scaff.annotation.gtf"
        :type gtf_file: str, optional
        :param pipeline_input_source: Pipeline input source for pipelining, options are ["auto", "map_win_2_pos", "generate_q_values", "filters"], defaults to "auto"
        :type pipeline_input_source: str, optional
        :param rerun: Rerun the analysis. If False, load previous output. Defaults to False.
        :type rerun: bool, optional
        :return: Mapped positions to genes in a list: ``[gene mapping dataframe, CCRE mapping dataframe]``
        :rtype: list[pd.DataFrame] 
        """
        assert (not self.pipeline and positions is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert pipeline_input_source in ["auto", "map_win_2_pos", "generate_q_values", "filters"], "Invalid parameter for pipeline_input_source. Options are: [""auto"", ""map_win_2_pos"", ""generate_q_values"", ""filters""]"

        parameters = locals().copy()

        if ref_folder!= None: gtf_file = Path(__file__).resolve().parent / ref_folder / gtf_file ; bed_file = Path(__file__).resolve().parent / ref_folder / bed_file

        if positions is not None:
            positions = positions.copy()
            positions.process()
            positions = positions.data_container
        else:
            # pipeline is in use here, due to assert and data is None
            if pipeline_input_source == "auto":
                if self.saved_results.get((self.filters.__name__, "position")) is not None:
                    positions = self.saved_results[(self.filters.__name__, "position")]
                elif self.saved_results.get((self.generate_q_values.__name__, "position")) is not None:
                    positions = self.saved_results[(self.generate_q_values.__name__, "position")]
                elif self.saved_results.get(self.map_win_2_pos.__name__) is not None:
                    positions = self.saved_results[self.map_win_2_pos.__name__]
            elif pipeline_input_source == "map_win_2_pos":
                positions = self.saved_results[self.map_win_2_pos.__name__]
            elif pipeline_input_source == "generate_q_values":
                positions = self.saved_results[(self.generate_q_values.__name__, "position")]
            else:
                positions = self.saved_results[(self.filters.__name__, "position")]
        
        
        parameters.pop("pipeline_input_source")
        parameters.pop("ref_folder")
        parameters = self.__prepare_parameters(parameters, positions=positions, gtf_file =gtf_file, bed_file= bed_file)

        res = self.obj.map_positions_to_genes(**parameters)

        if self.pipeline:
            self.saved_results[(self.map_positions_to_genes.__name__, "gene")] = res[0]
            self.saved_results[(self.map_positions_to_genes.__name__, "CCRE")] = res[1]

        return tuple(x.reset_index() for x in res)
    
    VOLCANO_PLOT_REQUIRED_COLUMNS = {
        "data": ["q-value", "diff"]
    }
    def volcano_plot(self, data: Optional[InputProcessor] = None, name : str = "volcano_plot.png", threshold : Optional[float] = 0.05, line: Optional[float] = None, x_range: tuple[int, int] = (-1, 1), y_max: int = None, title: str = None, x_label: str = None, y_label: str = None, position_or_window: str = "auto") -> None:
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
        :param line: Vertical line threshold for the `abs(line)` vertical line. Set to NOne to have no lines, defaults to None
        :type line: float, optional
        :param x_range: X-axis range, defaults to (-1, 1)
        :type x_range: tuple[int, int], optional
        :param y_max: Y-axis maximum, defaults to None
        :type y_max: int, optional
        :param title: Plot title, defaults to None for a generic title
        :type title: str, optional
        :param x_label: X-axis label, defaults to None for a generic label
        :type x_label: str, optional
        :param y_label: Y-axis label, defaults to None for a generic label
        :type y_label: str, optional
        :param position_or_window: The position-based or window-based results to use as input if DiffMethylTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto"
        :type position_or_window: str, optional
        """
        name = self.results_path + "/" + name
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
                if self.saved_results.get((self.generate_q_values.__name__, "window")) is not None:
                    data = self.saved_results[(self.generate_q_values.__name__, "window")]
                elif self.saved_results.get((self.generate_q_values.__name__, "position")) is not None:
                    data = self.saved_results[(self.generate_q_values.__name__, "position")]
                else:
                    # TODO invalid
                    pass
            elif position_or_window == "position":
                data = self.saved_results[(self.generate_q_values.__name__, "position")]
            else:
                data = self.saved_results[(self.generate_q_values.__name__, "window")]
            
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
        :param position_or_window: The position-based or window-based results to use as input if DiffMethylTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto"
        :type position_or_window: str, optional
        """        
        name = self.results_path + "/" + name
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
                if self.saved_results.get((self.generate_q_values.__name__, "window")) is not None:
                    data = self.saved_results[(self.generate_q_values.__name__, "window")]
                elif self.saved_results.get((self.generate_q_values.__name__, "position")) is not None:
                    data = self.saved_results[(self.generate_q_values.__name__, "position")]
                else:
                    # TODO invalid
                    pass
            elif position_or_window == "position":
                data = self.saved_results[(self.generate_q_values.__name__, "position")]
            else:
                data = self.saved_results[(self.generate_q_values.__name__, "window")]
                
        parameters.pop("position_or_window")
        parameters = self.__prepare_parameters(parameters, data=data)
        
        self.plots.manhattan_plot(**parameters)

    COVERAGE_PLOT_REQUIRED_COLUMNS = {
        "case_data": ["coverage", "positive_methylation_count", "negative_methylation_count"],
        "ctr_data": ["coverage", "positive_methylation_count", "negative_methylation_count"]
    }
    def coverage_plot(self, case_data: InputProcessor, ctr_data: InputProcessor, name : str = "coverage_plot.png", cov_min : int = 1, cov_max : int = -1, cov_max_percentile : float = 99.5, bins:int = 20, title: str = None, x_label: str = None, y_label: str = None) -> None:
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
        :param cov_min: Minimum coverage display, defaults to 1
        :type cov_min: int, optional
        :param cov_max: Maximum coverage display, defaults to -1
        :type cov_max: int, optional
        :param cov_max_percentile: Maximum coverage percentile display. Ranges from 0.0-100.0, defaults to 99.5. Overrides ``cov_max`` if set.
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
        name = self.results_path + "/" + name
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
    def graph_gene_regions(self, gene_data: Optional[InputProcessor] = None, ccre_data: Optional[InputProcessor] = None, name: str ="gene_regions.png", gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], intron_cutoff: int = -1, exon_cutoff: int = -1, upstream_cutoff: int = -1, CCRE_cutoff: int = -1, prom_cutoff:int = -1, title: str = None, x_label: str = None, intron_y_label: str = None, exon_y_label: str = None, upstream_y_label: str = None, CCRE_y_label: str = None, prom_y_label: str = None , position_or_window: str = "auto") -> None:
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
        :param position_or_window: The position-based or window-based results to use as input if DiffMethylTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto"
        :type position_or_window: str, optional
        """
        name = self.results_path + "/" + name
        assert (not self.pipeline and (gene_data is not None or ccre_data is not None)) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert position_or_window in ["auto", "position", "window"], "Invalid parameter for position_or_window. Options are: [""auto"", ""position"", ""window""]"

        parameters = locals().copy()

        if gene_data is not None:
            gene_data = gene_data.copy()
            gene_data.process()
            gene_data = gene_data.data_container
        else:
            if position_or_window == "auto":
                if self.saved_results.get((self.map_positions_to_genes.__name__, "gene")) is not None:
                    gene_data = self.saved_results[(self.map_positions_to_genes.__name__, "gene")]
                elif self.saved_results.get((self.map_windows_to_genes.__name__, "gene")) is not None:
                    gene_data = self.saved_results[(self.map_windows_to_genes.__name__, "gene")]
                else:
                    # TODO invalid
                    pass
            elif position_or_window == "position":
                gene_data = self.saved_results[(self.map_positions_to_genes.__name__, "gene")]
            else:
                gene_data = self.saved_results[(self.map_windows_to_genes.__name__, "gene")]

        if ccre_data is not None:
            ccre_data = ccre_data.copy()
            ccre_data.process()
            ccre_data = ccre_data.data_container
        else:
            if position_or_window == "auto":
                if self.saved_results.get((self.map_positions_to_genes.__name__, "CCRE")) is not None:
                    ccre_data = self.saved_results[(self.map_positions_to_genes.__name__, "CCRE")]
                elif self.saved_results.get((self.map_windows_to_genes.__name__, "CCRE")) is not None:
                    ccre_data = self.saved_results[(self.map_windows_to_genes.__name__, "CCRE")]
                else:
                    # TODO invalid
                    pass
            elif position_or_window == "position":
                ccre_data = self.saved_results[(self.map_positions_to_genes.__name__, "CCRE")]
            else:
                ccre_data = self.saved_results[(self.map_windows_to_genes.__name__, "CCRE")]
        
        parameters.pop("position_or_window")
        # print(gene_data)
        # print(ccre_data)
        parameters = self.__prepare_parameters(parameters, gene_data=gene_data, ccre_data=ccre_data)

        res = self.plots.graph_gene_regions(**parameters)

    GRAPH_UPSTREAM_GENE_METHYLATION_REQUIRED_COLUMNS = {
        "position_data": ["chromosome", "position_start", "diff"]
        # "gene_data": ["intron", "intron_diff", "exon", "exon_diff", "upstream", "upstream_diff"]
    }
    def graph_upstream_gene_methylation(self, position_data: Optional[InputProcessor] = None, ref_folder:str = None, region_data: Optional[InputProcessor] = None, name: str = "upstream_methylation.png", csv_name: str = "upstream_methylation.csv", csv: Optional[str] = None, left_distance: int = 1000, right_distance: int = 1000, window_size: int = 100, hypermethylated: bool = True, gene_hypermethylated_min: int = 20, window_hypermethylated_min: int = 5, min_hypermethylated_windows: int = 5, hypomethylated: bool = True, gene_hypomethylated_max: int = -20, window_hypomethylated_max: int = -5, min_hypomethylated_windows: int = 5, position_count: int = 5, clamp_positive: int = 50, clamp_negative: int = -50, title:str = None, gtf_file: str= "gencode.chr_patch_hapl_scaff.annotation.gtf", position_or_window: str = "auto", position_or_region:str = "region") -> None:
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
        :param name: Output PNG file name, defaults to "upstream_methylation.png"
        :type name: str, optional
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
        :param gene_hypermethylated_min: Minimum hypermethylation for genes, defaults to 0.20
        :type gene_hypermethylated_min: int, optional
        :param window_hypermethylated_min: Minimum hypermethylation for windows, defaults to 0.05
        :type window_hypermethylated_min: int, optional
        :param min_hypermethylated_windows: Minimum hypermethylated windows, defaults to 5
        :type min_hypermethylated_windows: int, optional
        :param hypomethylated: Include hypomethylated regions, defaults to True
        :type hypomethylated: bool, optional
        :param gene_hypomethylated_max: Maximum hypomethylation for genes, defaults to -0.20
        :type gene_hypomethylated_max: int, optional
        :param window_hypomethylated_max: Maximum hypomethylation for windows, defaults to -0.5
        :type window_hypomethylated_max: int, optional
        :param min_hypomethylated_windows: Minimum hypomethylated windows, defaults to 5
        :type min_hypomethylated_windows: int, optional
        :param position_count: Position count, defaults to 5
        :type position_count: int, optional
        :param clamp_positive: Positive clamp value, defaults to 0.50
        :type clamp_positive: int, optional
        :param clamp_negative: Negative clamp value, defaults to -0.50
        :type clamp_negative: int, optional
        :param title: Plot title, defaults to None for a generic title
        :type title: str, optional
        :param gtf_file: GTF file, defaults to "gencode.chr_patch_hapl_scaff.annotation.gtf"
        :type gtf_file: str, optional
        :param position_or_window: Position or window, options are ["auto", "position", "window"], defaults to "auto"
        :type position_or_window: str, optional
        """
        name = self.results_path + "/" + name

        if ref_folder!= None: gtf_file = Path(__file__).resolve().parent / ref_folder / gtf_file
        assert (not self.pipeline and position_data is not None) or (not self.pipeline and csv is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        assert position_or_window in ["auto", "position", "window"], "Invalid parameter for position_or_window. Options are: [""auto"", ""position"", ""window""]"
        assert position_or_region in ["position", "region"], "Invalid parameter for position_or_window. Options are: [""position"", ""region""]"

        parameters = locals().copy()


        #if gene_data is not None:
        #    gene_data = gene_data.copy()
        #    gene_data.process()
        #    gene_data = gene_data.data_container
        #    
        #elif csv is None:
        #    if position_or_window == "auto":
        #        if self.saved_results.get((self.map_positions_to_genes.__name__, "gene")) is not None:
        #            gene_data = self.saved_results[(self.map_positions_to_genes.__name__, "gene")]
        #        elif self.saved_results.get((self.map_windows_to_genes.__name__, "gene")) is not None:
        #            gene_data = self.saved_results[(self.map_windows_to_genes.__name__, "gene")]
        #        else:
        #            # TODO invalid
        #            pass
        #    elif position_or_window == "position":
        #        gene_data = self.saved_results[(self.map_positions_to_genes.__name__, "gene")]
        #    else:
        #        gene_data = self.saved_results[(self.map_windows_to_genes.__name__, "gene")]

        if position_data is not None:
            position_data = position_data.copy()
            position_data.process()
            position_data = position_data.data_container
        elif csv is None:
            position_data = self.saved_results[self.merge_tables.__name__]

        if region_data is not None:
            region_data = region_data.copy()
            region_data.process()
            region_data = region_data.data_container
        elif csv is None:
            region_data = self.saved_results[(self.generate_DMR.__name__, "cluster_df")]

        parameters.pop("position_or_window")
        parameters.pop("ref_folder")
        parameters = self.__prepare_parameters(parameters, position_data=position_data, region_data=region_data)

        self.plots.graph_upstream_gene_methylation(**parameters)


    PROCESS_REGIONS_REQUIRED_COLUMNS = {
        "region_file": ["chromosome","start","end"]
    }
    def process_regions(self, region_file: Optional[InputProcessor] = None, ref_folder:str = None, annotation_file:str = "CpG_gencode_annotation.bed", gene_bed_file:str = "gencode.v42.chr_patch_hapl_scaff.annotation.genes.bed", ccre_file:str ="encodeCcreCombined.bed") -> list:
        if ref_folder!= None : annotation_file = Path(__file__).resolve().parent /ref_folder / annotation_file ; gene_bed_file = Path(__file__).resolve().parent / ref_folder / gene_bed_file; ccre_file = Path(__file__).resolve().parent / ref_folder / ccre_file

        assert (not self.pipeline and region_file is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided."
        parameters = locals().copy()


        if region_file is not None:
            region_file = region_file.copy()
            region_file.process()
            region_file = region_file.data_container
        else:
            region_file = self.saved_results[(self.generate_DMR.__name__, "cluster_df")]

        parameters = self.__prepare_parameters(parameters, region_file = region_file)

        parameters.pop("ref_folder")

        res = self.obj.process_regions(**parameters)
        return res

    GRAPH_UPSTREAM_UCSC_REQUIRED_COLUMNS = {
        "position_data": ["chromosome", "position_start", "methylation_percentage*"]
    }
    def graph_upstream_UCSC(self, gene_name: str, position_data: Optional[InputProcessor] = None, ref_folder:str = None , name: str="UCSC_graph.bedGraph", before_tss: int = 5000, gtf_file: str = "gencode.chr_patch_hapl_scaff.annotation.gtf") -> None:
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
        :param gtf_file: GTF file, defaults to "gencode.chr_patch_hapl_scaff.annotation.gtf"
        :type gtf_file: str, optional
        """
        name = self.results_path + "/" + name


        if ref_folder != None: gtf_file = Path(__file__).resolve().parent / ref_folder / gtf_file ; bed_file = Path(__file__).resolve().parent / ref_folder / bed_file


        assert (not self.pipeline and position_data is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        parameters = locals().copy()

        if position_data is not None:
            position_data = position_data.copy()
            position_data.process()
            position_data = position_data.data_container
        else:
            position_data = self.saved_results[self.merge_tables.__name__]
        
        parameters.pop("ref_folder")
        parameters = self.__prepare_parameters(parameters, position_data=position_data)

        self.plots.graph_upstream_UCSC(**parameters)

    GRAPH_FULL_GENE_REQUIRED_COLUMNS = {
        "position_data": ["chromosome", "position_start", "methylation_percentage*"]
    }
    def graph_full_gene(self, gene_name:str, position_data: Optional[InputProcessor] = None, ref_folder:str = None, name="gene_methylation_graph.png", before_tss: int = 0, after_tss: Optional[int] = None, bin_size: int = 500, start_marker: bool = True, end_marker: bool = True, deviation_display: bool = True, aggregate_samples: bool=True, legend_size:int = 12, title: str = None, x_label:str = None, y_label:str=None, case_name: str = "Case", ctr_name: str = "Control",  gtf_file: str = "gencode.chr_patch_hapl_scaff.annotation.gtf") -> None:
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
        :param aggregate_samples: Aggregate samples and display mean, defaults to True
        :type aggregate_samples: bool, optional
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
        :param gtf_file: GTF file, defaults to "gencode.chr_patch_hapl_scaff.annotation.gtf"
        :type gtf_file: str, optional
        """
        name = self.results_path + "/" + name

        if ref_folder != None: gtf_file = Path(__file__).resolve().parent /ref_folder / gtf_file ; bed_file = Path(__file__).resolve().parent / ref_folder / bed_file

        assert (not self.pipeline and position_data is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided." 
        parameters = locals().copy()

        if position_data is not None:
            position_data = position_data.copy()
            position_data.process()
            position_data = position_data.data_container
        else:
            position_data = self.saved_results[self.merge_tables.__name__]
        
        parameters.pop("ref_folder")
        parameters = self.__prepare_parameters(parameters, position_data=position_data)

        self.plots.graph_full_gene(**parameters)

    PLOT_METHYLATION_CURVE_REQUIRED_COLUMNS ={
        "region_data":['chromosome', 'start', 'end'],
	"position_data": ['chrom', 'chromStart', 'blockSizes_case*', 'blockSizes_ctr*']
    }
    def plot_methylation_curve(self, region_data: Optional[InputProcessor] = None, ref_folder:str = None, position_data: Optional[InputProcessor] = None, name:str = ".", repeat_regions_df: str = "rmsk.txt", enhancer_promoter_df: str = "encodeCcreCombined.bed", repeat_regions_columns:list[int] = [5,6,7,11], enhancer_promoter_columns:list[int] = [0,1,2,12,13], window_size:int = 50, step_size:int = 25, chr_filter:str = None, start_filter:int = None, end_filter:int = None, sample_start_ind:int = 3) -> dict:
        
        assert (not self.pipeline and region_data is not None and position_data is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided."

        if ref_folder != None: repeat_regions_df = Path(__file__).resolve().parent / ref_folder / repeat_regions_df ; enhancer_promoter_df = Path(__file__).resolve().parent / ref_folder / enhancer_promoter_df

        parameters = locals().copy()
        if region_data is not None:
            region_data = region_data.copy()
            region_data.process()
            region_data = region_data.data_container
        else:
            region_data = self.saved_results[(self.generate_DMR.__name__, "cluster_df")]
        if position_data is not None:
            position_data = position_data.copy()
            position_data.process()
            position_data = position_data.data_container
        else:
            position_data = self.saved_results[(self.filters.__name__, "position")]	

        parameters.pop("ref_folder")
        parameters = self.__prepare_parameters(parameters, region_data=region_data, position_data=position_data)
        res = self.plots.plot_methylation_curve(**parameters)
        return pd.DataFrame(res)

    ALL_ANALYSIS_REQUIRED_COLUMNS = {
        "case_data": ["chromosome", "position_start", "coverage", "methylation_percentage", "positive_methylation_count", "negative_methylation_count", "strand"],
        "ctr_data": ["chromosome", "position_start", "coverage", "methylation_percentage", "positive_methylation_count", "negative_methylation_count", "strand"]
    }

    def all_analysis(self, case_data: InputProcessor, ctr_data: InputProcessor, ref_folder = None,window_based=False, min_cov_individual = 10, min_cov_group = 15, filter_samples_ratio=0.6, meth_group_threshold=0.2, cov_percentile = 100.0, min_samp_ctr = 2, min_samp_case = 2, max_q_value=0.05, abs_min_diff=0.0) -> pd.DataFrame:
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
        :param min_cov_individual: Minimum coverage filter (individual), defaults to 10
        :type min_cov_individual: int, optional
        :type min_cov: int, optional
        :param min_cov_group: Minimum coverage filter (group), defaults to 15
        :param filter_samples_ratio: Minimum sample ratio filter. Used with min_cov_group, defaults to 0.6
        :type filter_samples_ratio: float, optional
        :param meth_group_threshold: Methylation group threshold. Used with min_cov_group, defaults to 0.2
        :type meth_group_threshold: float, optional
        :type min_cov_group: int, optional
        :param cov_percentile: Maximum coverage filter (percentile of sample coverage). Ranges from 0.0-100.0, defaults to 100.0
        :type cov_percentile: float, optional
        :param min_samp_ctr: Minimum samples in control, defaults to 2
        :type min_samp_ctr: int, optional
        :param min_samp_case: Minimum samples in case, defaults to 2
        :type min_samp_case: int, optional
        :param max_q_value: Maximum q-value filter, defaults to 0.05
        :type max_q_value: float, optional
        :param abs_min_diff: Minimum absolute difference filter, defaults to 0.25
        :type abs_min_diff: float, optional
        :return: Final significant positions
        :rtype: pd.DataFrame
        """
        self.pipeline = False # True may take a lot of memory
        print("merging")
        min_cov_individual = int(min_cov_individual)
        min_cov_group = int(min_cov_group)
        filter_samples_ratio = float(filter_samples_ratio)
        meth_group_threshold = float(meth_group_threshold)
        cov_percentile = float(cov_percentile)
        min_samp_ctr = int(min_samp_ctr)
        min_samp_case = int(min_samp_case)
        merged = self.merge_tables(case_data, ctr_data, min_cov_individual = min_cov_individual, min_cov_group = min_cov_group, filter_samples_ratio=filter_samples_ratio, meth_group_threshold=meth_group_threshold, cov_percentile = cov_percentile, min_samp_ctr = min_samp_ctr, min_samp_case = min_samp_case)
        del case_data, ctr_data
        if window_based:
            print("windows")
            res = self.window_based(InputProcessor(merged))
            print("p-val")
            # res = self.generate_q_values(InputProcessor(res))
            print("filter")
            res = self.filters(InputProcessor(res), max_q_value=max_q_value, abs_min_diff=abs_min_diff)
            print("mapping")
            res = self.map_win_2_pos(InputProcessor(res), InputProcessor(merged))
        else:
            res = self.position_based(InputProcessor(merged))
            del merged
            # res = self.generate_q_values(InputProcessor(res))
            res_filter = self.filters(InputProcessor(res), max_q_value=max_q_value, abs_min_diff=abs_min_diff)
            ####################################
        DMR = self.generate_DMR(InputProcessor(res_filter), InputProcessor(res))
        pos_mapped = self.map_win_2_pos(InputProcessor(DMR[0]) , InputProcessor(res) )
        mapped = self.map_positions_to_genes(InputProcessor(pos_mapped), ref_folder= ref_folder)
        return mapped
    
    ALL_PLOTS_REQUIRED_COLUMNS = {
        "data": ["chromosome", "position_start", "region_start", "q-value", "diff", "methylation_percentage*", "coverage"],
        "gene_data": ["intron", "intron_diff", "exon", "exon_diff", "upstream", "upstream_diff"],
        "ccre_data": ["CCRE", "CCRE_diff"]
    }
    def all_plots(self, data: InputProcessor, ref_folder: str,  window_data: InputProcessor, gene: InputProcessor, ccre: InputProcessor) -> None:
        """Run all plot methods."""
        self.pipeline = False
        self.volcano_plot(data) #
        self.manhattan_plot(data) #
        #try:
        self.graph_upstream_gene_methylation(position_data=data, region_data=window_data, ref_folder = ref_folder,  position_count = 50, min_hypomethylated_windows= 20, min_hypermethylated_windows = 20, left_distance = 4000, right_distance = 100, clamp_negative = -100, clamp_positive = 100)
        #except:
        #    print("No graph_upstream_gene_methylation generated")
        # self.pie_chart(gene_data, ccre_data)
        self.graph_gene_regions(gene, ccre)
    MATCH_REGION_ANNOTATION_REQUIRED_COLUMNS ={
        "regions_df":['chrom', 'chromStart', 'chromEnd'],
    }
    def match_region_annotation(self, regions_df: Optional[InputProcessor] = None, ref_folder:str = None, bed_file:str = "CpG_gencode_annotation.bed", name:str="match_region_annotation", annotation_or_region: str = "region", show_counts: bool = False) -> list:

        assert (not self.pipeline and regions_df is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided."
        assert annotation_or_region in ["annotation", "region"], "Invalid parameter for annotation_or_region. Options are: [""annotation"", ""region""]"

        if ref_folder!= None : bed_file = Path(__file__).resolve().parent / ref_folder / bed_file
        parameters = locals().copy()
        if regions_df is not None:
            regions_df = regions_df.copy()
            regions_df.process()
            regions_df = regions_df.data_container
        else:
            regions_df = self.saved_results[(self.generate_DMR.__name__, "cluster_df")]

        parameters.pop("ref_folder")
        parameters = self.__prepare_parameters(parameters, regions_df=regions_df)
        res = self.plots.match_region_annotation(**parameters)
        return res    
    MATCH_POSITION_ANNOTATION_REQUIRED_COLUMNS ={
        "regions_df":['chrom', 'chromStart'],
    }
    def match_position_annotation(self, regions_df: Optional[InputProcessor] = None, ref_folder:str = None, bed_file:str = "CpG_gencode_annotation.bed", name:str="match_position_annotation") -> list:
        assert (not self.pipeline and regions_df is not None) or (self.pipeline), "If the pipeline isn't in use, data must be provided."
        if ref_folder!= None : bed_file = Path(__file__).resolve().parent / ref_folder / bed_file
        parameters = locals().copy()
        if regions_df is not None:
            regions_df = regions_df.copy()
            regions_df.process()
            regions_df = regions_df.data_container
        else:
            regions_df = self.saved_results[(self.generate_DMR.__name__, "cluster_df")]

        parameters.pop("ref_folder")
        parameters = self.__prepare_parameters(parameters, regions_df=regions_df)
        res = self.plots.match_position_annotation(**parameters)
        return res

def parse_arguments():
    parser = argparse.ArgumentParser(description="CLI for DiffMethylTools")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    obj = DiffMethylTools(pipeline=False)
    # Dynamically add subcommands for each method in DiffMethylTools
    for name, method in inspect.getmembers(obj, predicate=inspect.ismethod):
        if name.startswith("_"):  # Skip private methods
            continue

        names_argument_added = []

        method_parser = subparsers.add_parser(name, help=f"Run the {name} method")
        sig = inspect.signature(method)
        for param_name, param in sig.parameters.items():
            # print(param_name, param) #######################################
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

    # Initialize DiffMethylTools instance
    tool = DiffMethylTools(pipeline=False)

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
    elif type(result) == tuple:
        for i, l in enumerate(list(result)):
            f = open(f"{output}_{i}.txt", "w")
            f.write("\n".join(l))
            f.close()

if __name__ == "__main__":
    main()

