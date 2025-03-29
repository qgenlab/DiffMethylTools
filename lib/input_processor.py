import pandas as pd
import polars as pl
from typing import Optional
from collections import deque
import re

# TODO maybe enforce types for columns?

ACCEPTABLE_COLUMNS = {"chromosome", "position_start", "position_end", "region_start", "region_end", "strand", 
                      "methylation_percentage", "coverage", 
                      "positive_methylation_count", "negative_methylation_count", "avg_case", "avg_ctr", 
                      "diff", "p-val", "nbr", "q-value", "gene", "intron", "intron_diff", "exon", "exon_diff", 
                      "upstream", "upstream_diff", "CCRE", "CCRE_diff", "tag", "nearest_genes"}
ACCEPTABLE_STYLES = {"DiffMethTools_input", "Bismark", "DiffMethTools_header"}

# Format pre-selected styles
DIFFMETHTOOLS_INPUT_STYLE = {"chromosome":0, "position_start":1, "position_end":2, "strand":5, "coverage":9, "methylation_percentage":10}
BISMARK_STYLE = {"chromosome":0, "position_start":1, "strand":2, "positive_methylation_count":3, "negative_methylation_count":4}

DIFFMETHTOOLS_HEADER_STYLE = {"chromosome":0, "position_start":1, "position_end":2, "strand":3}

class FormatDefinition():
    """Example Usage:
    
    >>> format = FormatDefinition(style="DiffMethTools_input", column_mapping={"coverage":8}, sep=",")
    >>> format.mapping
    {'chromosome': 0, 'start': 1, 'end': 2, 'strand': 5, 'coverage': 8, 'methylation_percentage': 10}
    >>> format.separator
    ','
    """
    def __init__(self, style:Optional[str]=None, column_mapping:dict[str, int] = {}, sep:Optional[str]="\t"):
        assert style is None or style in ACCEPTABLE_STYLES, f"Acceptable format styles are: {str(ACCEPTABLE_STYLES)}."

        final_mapping = {}
        if style == "DiffMethTools_input":
            final_mapping = DIFFMETHTOOLS_INPUT_STYLE
        elif style == "Bismark":
            final_mapping = BISMARK_STYLE
        elif style == "DiffMethTools_header":
            final_mapping = DIFFMETHTOOLS_HEADER_STYLE
        
        # Merge style and column_mapping. Overwrite any style values with what's in column_mapping.
        final_mapping.update(column_mapping)

        # force to be list
        if "methylation_percentage" in final_mapping and not isinstance(final_mapping["methylation_percentage"], list):
            final_mapping["methylation_percentage"] = [final_mapping["methylation_percentage"]]
        if "coverage" in final_mapping and not isinstance(final_mapping["coverage"], list):
            final_mapping["coverage"] = [final_mapping["coverage"]]
            final_mapping["coverage_KEY"] = final_mapping["coverage"]
            del final_mapping["coverage"]

        self.final_mapping = final_mapping
        self.sep = sep

    @property 
    def mapping(self):
        return self.final_mapping
    
    @property
    def separator(self):
        return self.sep
    
    def __str__(self):
        return str(self.final_mapping)

# {"chromosome", "start", "end", "strand", "methylation_percentage", "coverage", "positive_methylation_count", "negative_methylation_count"}
INPUT_COLUMN_RENAME = {
    "chromosome":"chrom", "position_start":"chromStart", "position_end":"chromEnd",
    "region_start":"start", "region_end":"end",
    "methylation_percentage.*":"blockSizes_{ctr_case}_{name}", "coverage_KEY.*":"coverage_{ctr_case}_{name}", 
    "positive_methylation_count":"positive_{name}", "negative_methylation_count":"negative_{name}", "q_value":"q-value", "p_val":"p-val"
    }

class InputProcessor():
    """Example Usage:

    >>> obj = input_processor([data_1,data_2,...,data_n], [format_1,format_2,...,format_n])
    >>> obj.process()
    >>> obj.data_container
    [processed_1,processed_2,...,processed_n]

    Format is only required if data is not in the correct format. Otherwise, it can be none.
    """
    def __init__(self, data:list[str]|list[pd.DataFrame]|str|pd.DataFrame, format:Optional[list[FormatDefinition]|FormatDefinition] = None, has_header:list[bool]|bool = True, names:Optional[list[str]|str]= None, ctr_case:Optional[list[str]|str]=None, zero_indexed_positions:bool = True):
        ## assert relationships
        # assert an n-to-n relationship if both are list.
        if format is not None:
            if isinstance(format, list) and isinstance(data, list):
                assert len(format) == len(data), "Invalid input. data and format should have one of the following relationships, regarding the length of lists: [n:n, n:1, 1:1]."
            # assert a 1-to-many or 1-1 relationship
            else:
                assert isinstance(format, FormatDefinition) and isinstance(data, (list, pd.DataFrame, str)), "Invalid input. data and format should have one of the following relationships, regarding the length of lists: [n:n, n:1, 1:1]."

        self.original_params = locals().copy()
        self.original_params.pop("self")

        if isinstance(names, str):
            names = [names]
        if isinstance(ctr_case, str):
            ctr_case = [ctr_case]

        coverage_names = names.copy() if names is not None else None
        methylation_names = names.copy() if names is not None else None
        positive_names = names.copy() if names is not None else None
        negative_names = names.copy() if names is not None else None
        
        coverage_ctr_case = ctr_case.copy() if ctr_case is not None else None
        methylation_ctr_case = ctr_case.copy() if ctr_case is not None else None

        if not isinstance(has_header, list):
            if isinstance(data, list):
                has_header = [has_header] * len(data)
            else: 
                has_header = [has_header]

        if not isinstance(data, list):
            data = [data]
        if not isinstance(methylation_names, list):
            methylation_names = [methylation_names]
        if not isinstance(coverage_names, list):
            coverage_names = [coverage_names]
        if not isinstance(positive_names, list):
            positive_names = [positive_names]
        if not isinstance(negative_names, list):
            negative_names = [negative_names]
        if not isinstance(methylation_ctr_case, list):
            methylation_ctr_case = [methylation_ctr_case]
        if not isinstance(coverage_ctr_case, list):
            coverage_ctr_case = [coverage_ctr_case]
        

        self.raw_data = data
        self.raw_format = format  
        self.methylation_names = deque(methylation_names) if methylation_names is not None else None
        self.methylation_names_copy = self.methylation_names.copy()
        self.no_methylation_names = 1
        self.coverage_names = deque(coverage_names) if coverage_names is not None else None
        self.coverage_names_copy = self.coverage_names.copy()
        self.no_coverage_names = 1
        self.positive_names = deque(positive_names) if positive_names is not None else None
        self.positive_names_copy = self.positive_names.copy()
        self.no_positive_names = 1
        self.negative_names = deque(negative_names) if negative_names is not None else None
        self.negative_names_copy = self.negative_names.copy
        self.no_negative_names = 1
        self.methylation_ctr_case = deque(methylation_ctr_case) if methylation_ctr_case is not None else None
        self.methylation_ctr_case_copy = self.methylation_ctr_case.copy()
        self.coverage_ctr_case = deque(coverage_ctr_case) if coverage_ctr_case is not None else None
        self.coverage_ctr_case_copy = self.coverage_ctr_case.copy()
        self.zero_indexed_positions = zero_indexed_positions
        self.has_header = deque(has_header) 

    def copy(self):
        return InputProcessor(**self.original_params)


    @property
    def data_container(self):
        return self.data
    
    def process(self):
        data = []
        
        # Helper function to read data from different sources
        def read_data(source, sep):
            if isinstance(source, pd.DataFrame):
                return pl.from_pandas(source)
            if isinstance(source, str):
                has_header = self.has_header.popleft() if len(self.has_header) != 0 else True
                return pl.read_csv(source, separator=sep, has_header=has_header)
            raise ValueError("Unsupported data source type")

        def rename_columns(df):
            # rename based on regex patterns in INPUT_COLUMN_RENAME.
            renamed_columns = {}
            
            row_to_attribute_name = {"methylation_percentage":"methylation", "coverage_KEY":"coverage", "positive_methylation_count":"positive", "negative_methylation_count":"negative"}
            for col in df.columns:

                new_col = col
                for pattern, replacement in INPUT_COLUMN_RENAME.items():
                    if re.match(pattern, col):
                        name = None
                        for x in ["methylation_percentage", "coverage_KEY", "positive_methylation_count", "negative_methylation_count"]:
                            if x in col:
                                names_list = getattr(self, f"{row_to_attribute_name[x]}_names")
                                if names_list is not None:
                                    if len(names_list) == 0:
                                        setattr(self, f"{row_to_attribute_name[x]}_names", getattr(self, f"{row_to_attribute_name[x]}_names_copy").copy())
                                    name = names_list.popleft()
                                else:
                                    name = f"file_{getattr(self, f'no_{row_to_attribute_name[x]}_names')}"
                                    setattr(self, f"no_{row_to_attribute_name[x]}_names", getattr(self, f"no_{row_to_attribute_name[x]}_names") + 1)
                                break
                        if name is None:
                            name = "N/A"

                        if "methylation_percentage" in col:
                            if self.methylation_ctr_case is not None:
                                if len(self.methylation_ctr_case) == 0:
                                    self.methylation_ctr_case = self.methylation_ctr_case_copy.copy()

                                ctr_case = self.methylation_ctr_case.popleft()
                            else:
                                ctr_case = "N/A"
                        elif "coverage_KEY" in col:
                            if self.coverage_ctr_case is not None:
                                if len(self.coverage_ctr_case) == 0:
                                    self.coverage_ctr_case = self.coverage_ctr_case_copy.copy()

                                ctr_case = self.coverage_ctr_case.popleft()
                            else:
                                ctr_case = "N/A"
                        else:
                            ctr_case = "N/A"


                        new_col = re.sub(f"^{pattern}$", replacement.format(ctr_case=ctr_case, name=name), col)

                        break
                renamed_columns[col] = new_col

            return df.rename(renamed_columns)

        for i, file in enumerate(self.raw_data):
            if isinstance(self.raw_format, list):
                format_def = self.raw_format[i]
            else:
                format_def = self.raw_format
            
            df = read_data(file, format_def.separator if format_def is not None else ",")

            # map names according to the ones from the FormatDefinition. 
            if format_def is not None and format_def.mapping != {}:
                
                selected_columns = []
                old_column_names = []

                for col, idx in format_def.mapping.items():
                    if col in ["methylation_percentage", "coverage_KEY"] and isinstance(idx, list):
                        for sub_idx in idx:
                            old_column_names.append(df.columns[sub_idx])
                            selected_columns.append(pl.col(df.columns[sub_idx]).alias(f"{col}_{sub_idx}"))
                    else:
                        old_column_names.append(df.columns[idx])
                        selected_columns.append(pl.col(df.columns[idx]).alias(col))

                # include columns not in the mapping
                for col in reversed(df.columns):
                    if col not in old_column_names:
                        selected_columns.insert(0, pl.col(col))
                df = df.select(selected_columns) 


            # remap any of those names to the conventions used in the program
            renamed_df = rename_columns(df)

            if not self.zero_indexed_positions:
                if "chromStart" in renamed_df.columns:
                    renamed_df = renamed_df.with_columns(pl.col("chromStart") + 1)
                if "chromEnd" in renamed_df.columns:
                    renamed_df = renamed_df.with_columns(pl.col("chromEnd") + 1)

            renamed_df = renamed_df.to_pandas()

            if "gene" in renamed_df.columns:
                renamed_df.set_index("gene", inplace=True, drop=True)

            data.append(renamed_df)

        self.data = data if len(data) > 1 else data[0]
