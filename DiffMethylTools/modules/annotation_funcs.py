from typing import Optional
from pathlib import Path
import pandas as pd
from ..lib import InputProcessor
from .main_funcs import analysis_function

class AnnotationMixin:
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
        parameters = self._prepare_parameters(parameters, positions=positions, gtf_file =gtf_file, bed_file= bed_file)

        res = self.obj.map_positions_to_genes(**parameters)

        if self.pipeline:
            self.saved_results[(self.map_positions_to_genes.__name__, "gene")] = res[0]
            self.saved_results[(self.map_positions_to_genes.__name__, "CCRE")] = res[1]

        return tuple(x.reset_index() for x in res)
    MATCH_REGION_ANNOTATION_REQUIRED_COLUMNS ={
        "regions_df":['chrom', 'chromStart', 'chromEnd'],
    }
    def match_region_annotation(self, regions_df: Optional[InputProcessor] = None, ref_folder:str = None, bed_file:str = "CpG_gencode_annotation.bed", name:str="match_region_annotation", annotation_or_region: str = "region", show_counts: bool = False) -> list:
        """
        Intersect identified regions (DMRs) with genomic annotations (e.g., GENCODE) to determine 
        genomic context (promoters, exons, introns, etc.).

        .. note::
            Required columns in ``regions_df``: ``['chrom', 'chromStart', 'chromEnd']``.
            If ``pipeline`` is enabled and ``regions_df`` is None, it defaults to the 'cluster_df' 
            from the ``generate_DMR`` step.

        :param regions_df: Input dataframe containing genomic regions to annotate.
        :type regions_df: InputProcessor, optional
        :param ref_folder: Path to the reference genome folder.
        :type ref_folder: str, optional
        :param bed_file: Filename of the annotation BED file located in the ref_folder, defaults to "CpG_gencode_annotation.bed"
        :type bed_file: str, optional
        :param name: Prefix for output files, defaults to "match_region_annotation"
        :type name: str, optional
        :param annotation_or_region: Determines the perspective of the output overlap ("annotation" or "region"), defaults to "region"
        :type annotation_or_region: str, optional
        :param show_counts: If True, returns count statistics of the annotations, defaults to False
        :type show_counts: bool, optional
        :return: A list of DataFrames containing annotated regions.
        :rtype: list
        """
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
        parameters = self._prepare_parameters(parameters, regions_df=regions_df)
        res = self.plots.match_region_annotation(**parameters)
        return res    
    MATCH_POSITION_ANNOTATION_REQUIRED_COLUMNS ={
        "regions_df":['chrom', 'chromStart'],
    }
    def match_position_annotation(self, regions_df: Optional[InputProcessor] = None, ref_folder:str = None, bed_file:str = "CpG_gencode_annotation.bed", name:str="match_position_annotation") -> list:
        """
        Intersect single genomic positions (DMLs) with genomic annotations to determine local context.

        .. note::
            Required columns in ``regions_df``: ``['chrom', 'chromStart']``.
            Similar to ``match_region_annotation``, but optimized for single-base coordinate matching
            rather than range-based intersection.

        :param regions_df: Input dataframe containing genomic positions to annotate.
        :type regions_df: InputProcessor, optional
        :param ref_folder: Path to the reference genome folder.
        :type ref_folder: str, optional
        :param bed_file: Filename of the annotation BED file, defaults to "CpG_gencode_annotation.bed"
        :type bed_file: str, optional
        :param name: Prefix for output files, defaults to "match_position_annotation"
        :type name: str, optional
        :return: A list of DataFrames containing annotated positions.
        :rtype: list
        """
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
        parameters = self._prepare_parameters(parameters, regions_df=regions_df)
        res = self.plots.match_position_annotation(**parameters)
        return res
