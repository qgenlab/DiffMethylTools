from typing import Optional
from pathlib import Path
import pandas as pd
from ..lib import InputProcessor
from .main_funcs import analysis_function
import warnings

class AnalysisMixin:
    MERGE_TABLES_REQUIRED_COLUMNS = {
        "case_data": ["chromosome", "position_start", "coverage", "methylation_percentage", "positive_methylation_count", "negative_methylation_count", "strand"],
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
        :param min_cov_group: Minimum coverage filter (group), defaults to 15
        :type min_cov_group: int, optional
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

        parameters = self._prepare_parameters(parameters, case_data=case_data.data_container, ctr_data=ctr_data.data_container)

        res = self.obj.merge_tables(**parameters)

        
        if self.pipeline:
            self.saved_results[self.merge_tables.__name__] = res

        return res
    POSITION_BASED_REQUIRED_COLUMNS = {
        "data": ["chromosome", "position_start", "methylation_percentage*"]
    }
    @analysis_function
    def position_based(self, data: Optional[InputProcessor] = None , method="limma", features=None, test_factor="Group", processes=12, model="eBayes", min_std=0.1, fill_na:bool=True, rerun=False) -> pd.DataFrame:
        """Perform position-based DML detection. Has options for using the gamma function, or the limma R package.

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
            
        parameters = self._prepare_parameters(parameters, data=data)

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
        parameters = self._prepare_parameters(parameters, window_data=window_data, position_data=position_data)

        res = self.obj.map_win_2_pos(**parameters)

        if self.pipeline:
            self.saved_results[self.map_win_2_pos.__name__] = res

        return res

    FILTERS_REQUIRED_COLUMNS = {
        "data": ["q-value", "diff"]
    }
    @analysis_function
    def filters(self, data: Optional[InputProcessor] = None, max_q_value=0.05, abs_min_diff=0.10, position_or_window: str = "auto", rerun=False) -> pd.DataFrame:
        """Filter positions by q-value and minimum difference.

        .. note::
            Required columns for ``data``:
                - ``["q-value", "diff"]``

        :param data: Input data. Not necessary if the pipeline is in use, defaults to None
        :type data: InputProcessor, optional
        :param max_q_value: Maximum q-value filter, defaults to 0.05
        :type max_q_value: float, optional
        :param abs_min_diff: Absolute minimum difference filter, defaults to 0.10
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
        parameters = self._prepare_parameters(parameters, data=data)

        res = self.obj.filters(**parameters)

        if self.pipeline:
            self.saved_results[(self.filters.__name__, position_or_window)] = res

        return res

    GENERATE_DMR_REQUIRED_COLUMNS = {
        "significant_position_data": ["chromosome", "position_start", "diff"],
        "position_data": ["chromosome", "position_start", "diff"]
    }
    @analysis_function
    def generate_DMR(self, significant_position_data: Optional[InputProcessor] = None, position_data: Optional[InputProcessor] = None, min_pos=3, neutral_change_limit=7.5, neutral_perc=30, opposite_perc=10, significant_position_pipeline: str = "auto", rerun=False) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        # TODO ask chris what to write here for documentation
        """Generate Differentially Methylated Regions (DMRs) by clustering DMLs.
        
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
        parameters = self._prepare_parameters(parameters, significant_position_data=significant_position_data, position_data=position_data)

        res = self.obj.generate_DMR(**parameters)

        if self.pipeline:
            self.saved_results[(self.generate_DMR.__name__, "cluster_df")] = res[0]
            self.saved_results[(self.generate_DMR.__name__, "unclustered_dms_df")] = res[1]
            self.saved_results[(self.generate_DMR.__name__, "clustered_dms_df")] = res[2]
        return res
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

        parameters = self._prepare_parameters(parameters, region_file = region_file)

        parameters.pop("ref_folder")

        res = self.obj.process_regions(**parameters)
        return res
    ALL_ANALYSIS_REQUIRED_COLUMNS = {
        "case_data": ["chromosome", "position_start", "coverage", "methylation_percentage", "positive_methylation_count", "negative_methylation_count", "strand"],
        "ctr_data": ["chromosome", "position_start", "coverage", "methylation_percentage", "positive_methylation_count", "negative_methylation_count", "strand"]
    }

    def all_analysis(self, case_data: InputProcessor, ctr_data: InputProcessor, ref_folder = None,window_based=False, min_cov_individual = 10, min_cov_group = 15, filter_samples_ratio=0.6, meth_group_threshold=0.2, cov_percentile = 100.0, min_samp_ctr = 2, min_samp_case = 2, max_q_value=0.05, abs_min_diff=0.0, features=None) -> pd.DataFrame:
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
        if ref_folder == None: warnings.warn("all_analysis requires a reference folder; if not provided, the script will stop after the DMR identification step.", category=UserWarning, stacklevel=2)
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
            res = self.position_based(InputProcessor(merged), features = features)
            del merged
            # res = self.generate_q_values(InputProcessor(res))
            res_filter = self.filters(InputProcessor(res), max_q_value=max_q_value, abs_min_diff=abs_min_diff)
            ####################################
        DMR = self.generate_DMR(InputProcessor(res_filter), InputProcessor(res))
        pos_mapped = self.map_win_2_pos(InputProcessor(DMR[0]) , InputProcessor(res) )
        mapped = self.map_positions_to_genes(InputProcessor(pos_mapped), ref_folder= ref_folder)
        return mapped
