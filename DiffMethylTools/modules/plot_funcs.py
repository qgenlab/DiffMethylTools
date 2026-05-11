from typing import Optional
from pathlib import Path
import pandas as pd
from ..lib import InputProcessor
from .main_funcs import analysis_function

class PlotsMixin:
    VOLCANO_PLOT_REQUIRED_COLUMNS = {
        "data": ["q-value", "diff"]
    }
    def volcano_plot(self, data: Optional[InputProcessor] = None, name : str = "plots/volcano_plot.png", threshold : Optional[float] = 0.05, line: Optional[float] = None, x_range: tuple[int, int] = (-1, 1), y_max: int = None, title: str = None, x_label: str = None, y_label: str = None, position_or_window: str = "auto") -> None:
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
        :param line: Vertical line threshold for the `abs(line)` vertical line. Set to None to have no lines, defaults to None
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
        parameters = self._prepare_parameters(parameters, data=data)

        self.plots.volcano_plot(**parameters)

    MANHATTAN_PLOT_REQUIRED_COLUMNS = {
        "data": ["chromosome", "q-value", "position_start", "region_start"]
    }
    def manhattan_plot(self, data: Optional[InputProcessor] = None, name: str = "plots/manhattan_plot.png", threshold: Optional[float] = 0.05, title: str = None, x_label: str = None, y_label: str = None, position_or_window: str = "auto") -> None:
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
        parameters = self._prepare_parameters(parameters, data=data)
        
        self.plots.manhattan_plot(**parameters)

    COVERAGE_PLOT_REQUIRED_COLUMNS = {
        "case_data": ["coverage", "positive_methylation_count", "negative_methylation_count"],
        "ctr_data": ["coverage", "positive_methylation_count", "negative_methylation_count"]
    }
    def coverage_plot(self, case_data: InputProcessor, ctr_data: InputProcessor, name : str = "plots/coverage_plot.png", cov_min : int = 1, cov_max : int = -1, cov_max_percentile : float = 99.5, bins:int = 20, title: str = None, x_label: str = None, y_label: str = None) -> None:
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

        parameters = self._prepare_parameters(parameters, case_data=case_data.data_container, ctr_data=ctr_data.data_container)

        self.plots.coverage_plot(**parameters)

    GRAPH_GENE_REGIONS_REQUIRED_COLUMNS = {
        "gene_data": ["intron", "intron_diff", "exon", "exon_diff", "upstream", "upstream_diff"],
        "ccre_data": ["CCRE", "CCRE_diff"]
    }
    def graph_gene_regions(self, gene_data: Optional[InputProcessor] = None, ccre_data: Optional[InputProcessor] = None, name: str ="plots/gene_regions.png", gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], intron_cutoff: int = -1, exon_cutoff: int = -1, upstream_cutoff: int = -1, CCRE_cutoff: int = -1, prom_cutoff:int = -1, title: str = None, x_label: str = None, intron_y_label: str = None, exon_y_label: str = None, upstream_y_label: str = None, CCRE_y_label: str = None, prom_y_label: str = None , position_or_window: str = "auto") -> None:
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
        parameters = self._prepare_parameters(parameters, gene_data=gene_data, ccre_data=ccre_data)

        res = self.plots.graph_gene_regions(**parameters)

    GRAPH_UPSTREAM_GENE_METHYLATION_REQUIRED_COLUMNS = {
        "position_data": ["chromosome", "position_start", "diff"]
        # "gene_data": ["intron", "intron_diff", "exon", "exon_diff", "upstream", "upstream_diff"]
    }
    def graph_upstream_gene_methylation(self, position_data: Optional[InputProcessor] = None, ref_folder:str = None, region_data: Optional[InputProcessor] = None, name: str = "plots/upstream_methylation.png", csv_name: str = "upstream_methylation.csv", csv: Optional[str] = None, left_distance: int = 1000, right_distance: int = 1000, window_size: int = 100, hypermethylated: bool = True, gene_hypermethylated_min: int = 20, window_hypermethylated_min: int = 5, min_hypermethylated_windows: int = 5, hypomethylated: bool = True, gene_hypomethylated_max: int = -20, window_hypomethylated_max: int = -5, min_hypomethylated_windows: int = 5, position_count: int = 5, clamp_positive: int = 50, clamp_negative: int = -50, title:str = None, gtf_file: str= "gencode.chr_patch_hapl_scaff.annotation.gtf", position_or_window: str = "auto", position_or_region:str = "region") -> None:
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
        parameters = self._prepare_parameters(parameters, position_data=position_data, region_data=region_data)

        self.plots.graph_upstream_gene_methylation(**parameters)
    
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
        parameters = self._prepare_parameters(parameters, position_data=position_data)

        self.plots.graph_upstream_UCSC(**parameters)

    GRAPH_FULL_GENE_REQUIRED_COLUMNS = {
        "position_data": ["chromosome", "position_start", "methylation_percentage*"]
    }
    def graph_full_gene(self, gene_name:str, position_data: Optional[InputProcessor] = None, ref_folder:str = None, name="plots/gene_methylation_graph.png", before_tss: int = 0, after_tss: Optional[int] = None, bin_size: int = 500, start_marker: bool = True, end_marker: bool = True, deviation_display: bool = True, aggregate_samples: bool=True, legend_size:int = 12, title: str = None, x_label:str = None, y_label:str=None, case_name: str = "Case", ctr_name: str = "Control",  gtf_file: str = "gencode.chr_patch_hapl_scaff.annotation.gtf") -> None:
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
        parameters = self._prepare_parameters(parameters, position_data=position_data)

        self.plots.graph_full_gene(**parameters)

    PLOT_METHYLATION_CURVE_REQUIRED_COLUMNS ={
        "region_data":['chromosome', 'start', 'end'],
	"position_data": ['chrom', 'chromStart', 'blockSizes_case*', 'blockSizes_ctr*']
    }
    def plot_methylation_curve(self, region_data: Optional[InputProcessor] = None, ref_folder:str = None, position_data: Optional[InputProcessor] = None, name:str = "plots/.", repeat_regions_df: str = "rmsk.txt", enhancer_promoter_df: str = "encodeCcreCombined.bed", repeat_regions_columns:list[int] = [5,6,7,11], enhancer_promoter_columns:list[int] = [0,1,2,12,13], window_size:int = 50, step_size:int = 25, chr_filter:str = None, start_filter:int = None, end_filter:int = None, sample_start_ind:int = 3) -> dict:
        """Generate plots showing methylation curves across specific genomic regions, including annotations for repeats and enhancers/promoters.

        .. note::
            If ``pipeline`` mode is disabled, both ``region_data`` and ``position_data`` must be provided manually.
            If ``pipeline`` mode is enabled, the function automatically retrieves results from ``generate_DMR`` and ``filters``.

        .. note::
            This method utilizes a sliding window approach (defined by ``window_size`` and ``step_size``) to smooth methylation values across regions.

        :param region_data: DMR or cluster data. If None and pipelined, uses 'cluster_df' from generate_DMR.
        :type region_data: InputProcessor, optional
        :param ref_folder: Path to the reference genome folder containing annotation files.
        :type ref_folder: str, optional
        :param position_data: All position-level data. If None and pipelined, uses output from filters.
        :type position_data: InputProcessor, optional
        :param name: Directory path or prefix where the resulting plots will be saved, defaults to "."
        :type name: str, optional
        :param repeat_regions_df: Filename for repeat regions annotation (e.g., RepeatMasker), defaults to "rmsk.txt"
        :type repeat_regions_df: str, optional
        :param enhancer_promoter_df: Filename for enhancer/promoter annotation, defaults to "encodeCcreCombined.bed"
        :type enhancer_promoter_df: str, optional
        :param repeat_regions_columns: List of column indices to extract from the repeat regions file, defaults to [5,6,7,11]
        :type repeat_regions_columns: list[int], optional
        :param enhancer_promoter_columns: List of column indices to extract from the enhancer/promoter file, defaults to [0,1,2,12,13]
        :type enhancer_promoter_columns: list[int], optional
        :param window_size: Size of the sliding window in base pairs for smoothing, defaults to 50
        :type window_size: int, optional
        :param step_size: Step size for the sliding window in base pairs, defaults to 25
        :type step_size: int, optional
        :param chr_filter: Optional chromosome name to restrict plotting to a specific chromosome.
        :type chr_filter: str, optional
        :param start_filter: Optional genomic start coordinate to filter regions.
        :type start_filter: int, optional
        :param end_filter: Optional genomic end coordinate to filter regions.
        :type end_filter: int, optional
        :param sample_start_ind: Column index where individual sample methylation data begins in the input matrix, defaults to 3
        :type sample_start_ind: int, optional
        :return: A dictionary (converted to DataFrame) containing the plotted region metrics.
        :rtype: dict
        """
        
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
        parameters = self._prepare_parameters(parameters, region_data=region_data, position_data=position_data)
        res = self.plots.plot_methylation_curve(**parameters)
        return pd.DataFrame(res)
    ALL_PLOTS_REQUIRED_COLUMNS = {
        "data": ["chromosome", "position_start", "region_start", "q-value", "diff", "methylation_percentage*", "coverage"],
        "gene_data": ["intron", "intron_diff", "exon", "exon_diff", "upstream", "upstream_diff"],
        "ccre_data": ["CCRE", "CCRE_diff"]
    }
    def all_plots(self, data: InputProcessor, ref_folder: str,  window_data: InputProcessor, gene: InputProcessor, ccre: InputProcessor) -> None:
        """
        Execute a suite of visualization methods including Volcano plots, Manhattan plots, 
        and gene-centric methylation graphs.

        .. note::
            This method forces ``self.pipeline = False`` to ensure that all plots are generated 
            using the specific data objects provided as arguments rather than cached state.

        :param data: Position-level methylation data (e.g., DMLs).
        :type data: InputProcessor
        :param ref_folder: Path to the reference genome folder containing necessary genomic annotations.
        :type ref_folder: str
        :param window_data: Window-based methylation data (e.g., DMRs).
        :type window_data: InputProcessor
        :param gene: Gene annotation data.
        :type gene: InputProcessor
        :param ccre: Candidate Cis-Regulatory Elements (cCRE) data.
        :type ccre: InputProcessor
        :return: None
        """
        self.pipeline = False
        self.volcano_plot(data) #
        self.manhattan_plot(data) #
        #try:
        self.graph_upstream_gene_methylation(position_data=data, region_data=window_data, ref_folder = ref_folder,  position_count = 50, min_hypomethylated_windows= 20, min_hypermethylated_windows = 20, left_distance = 4000, right_distance = 100, clamp_negative = -100, clamp_positive = 100)
        #except:
        #    print("No graph_upstream_gene_methylation generated")
        # self.pie_chart(gene_data, ccre_data)
        # self.match_position_annotation(window_data, gene_data, ccre_data)
        self.graph_gene_regions(gene, ccre)
