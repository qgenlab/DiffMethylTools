import os
import yaml
import pandas as pd
from functools import wraps
from ..lib import InputProcessor, Plots, FormatDefinition, Analysis

def analysis_function(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if os.path.exists(self.results_path + "/data/state.yaml"):
            with open(self.results_path + "/data/state.yaml", "r") as file:
                state = yaml.safe_load(file)

            func_state = state.get(func.__name__, False)
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
                        result = tuple(res)
                    return result
                else:
                    state[func.__name__] = False
                    with open(self.results_path + "/data/state.yaml", "w") as file:
                        yaml.dump(state, file)

        if self.results_path is not None:  
            os.makedirs(self.results_path, exist_ok=True)
            os.makedirs(self.results_path + "/plots", exist_ok=True)
            os.makedirs(self.results_path + "/data", exist_ok=True)

        png_name = kwargs.get("name", None)
        csv_name = kwargs.get("csv_name", None)

        if png_name is not None and "/" not in png_name:
            kwargs["name"] = self.results_path + "/plots/" + png_name

        if csv_name is not None and "/" not in csv_name:
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
                        if func.__name__ in ["map_positions_to_genes", "map_windows_to_genes"]:
                            if i == 0:
                                dataframe.to_csv(self.results_path + "/data/" + func.__name__ + "_genes.csv", index=False)
                            else:
                                dataframe.to_csv(self.results_path + "/data/" + func.__name__ + "_CCRE.csv", index=False)
                        else:
                            dataframe.to_csv(self.results_path + "/data/" + func.__name__ + "_" + str(i) + ".csv", index=False)
                    else:
                        if func.__name__ in ["map_positions_to_genes", "map_windows_to_genes"]:
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

class MainMixin:
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

    def _prepare_parameters(self, parameters, **kwargs):
        parameters.pop("self", None)
        if "rerun" in parameters.keys():
            parameters.pop("rerun")
        for key in parameters:
            if key in kwargs.keys():
                parameters[key] = kwargs[key]
        return parameters