import argparse
import pandas as pd
from typing import Callable, Optional, Type, get_args, get_origin, Union
from types import UnionType
from .lib import InputProcessor, Plots, FormatDefinition, Analysis
import inspect
import os
import subprocess
import sys
import warnings

from .modules.main_funcs import MainMixin
from .modules.analysis_funcs import AnalysisMixin
from .modules.annotation_funcs import AnnotationMixin
from .modules.plot_funcs import PlotsMixin
from .modules import cpg_annotation

class DiffMethylTools(MainMixin, AnalysisMixin, AnnotationMixin, PlotsMixin):
    """
    Main pipeline handler for DiffMethylTools.
    Core functions are inherited from the respective modules for maintainability.
    """
    pass

# def run_setup(genome):
#     """Runs the setup scripts directly inside the package folder."""
#     package_home = os.path.dirname(os.path.abspath(__file__))
#     script_path = os.path.join(package_home, "bin", f"get_files_{genome}.sh")
#     if not os.path.exists(script_path):
#         print(f"Error: Could not find setup script at {script_path}")
#         sys.exit(1)
#     print(f"Downloading {genome} reference data into {package_home}...")
#     try:
#         subprocess.run(
#             ["bash", script_path], 
#             cwd=package_home, 
#             check=True
#         )
#         print(f"Success! {genome} data is ready.")
#     except subprocess.CalledProcessError as e:
#         print(f"Error during setup: {e}")
#         sys.exit(1)

def run_setup(genome):
    """Runs the setup scripts or the annotation pipeline based on the input."""
    if str(genome).lower().endswith(('.yml', '.yaml')):
        print(f"Running annotation pipeline using config: {genome}")
        try:
            annotator = cpg_annotation.run_setup(genome)
            # annotator.generate("Final_CpG_Annotated.bed")
            print("Success! Annotation pipeline completed.")
        except Exception as e:
            print(f"\n[ERROR] Pipeline Halted: {str(e)}")
            sys.exit(1)
        return 
    package_home = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(package_home, "bin", f"get_files_{genome}.sh")
    if not os.path.exists(script_path):
        print(f"Error: Could not find setup script at {script_path}")
        sys.exit(1)
    print(f"Downloading {genome} reference data into {package_home}...")
    try:
        subprocess.run(
            ["bash", script_path], 
            cwd=package_home, 
            check=True
        )
        print(f"Success! {genome} data is ready.")
    except subprocess.CalledProcessError as e:
        print(f"Error during setup: {e}")
        sys.exit(1)




def parse_arguments():
    parser = argparse.ArgumentParser(description="DiffMethylTools")
    
    parser.add_argument(
    "--version",
    action="version",
    version="%(prog)s 1.1.0"
    )

    parser.add_argument("--setup", help="Download required reference files for the specified genome.")
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    obj = DiffMethylTools(pipeline=False)
    # Dynamically add subcommands for each method in DiffMethylTools
    for name, method in inspect.getmembers(obj, predicate=inspect.ismethod):
        if name.startswith("_"):  # Skip private methods
            continue

        names_argument_added = []

        method_parser = subparsers.add_parser(name, help=f"Run the {name} method")
        sig = inspect.signature(method)
        if name == "all_analysis" or name == "merge_tables":
            method_parser.add_argument(f"--input_format", type=str, default=None, help=f"(Input the input file format (CR for Bismark cytosine report or BED for BED methylation file) default: {None})")
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

    if args.setup:
        run_setup(args.setup)
        sys.exit(0) 
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
            input_format = getattr(args, "input_format", None)

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
            if input_format != None:
                format = FormatDefinition(input_format)

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

    # print(method_args)
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

