import inspect
import re
import os
from DiffMethylTools import DiffMethylTools

# Define the groups and their target filenames
GROUPS = {
    "Analysis": {
        "file": "analysis.md",
        "desc": "Core commands for processing and analyzing methylation data.",
        "commands": ["all_analysis", "merge_tables", "position_based", "process_regions", "generate_DMR", "filters", "generate_q_values", "map_win_2_pos"]
    },
    "Annotation": {
        "file": "annotation.md",
        "desc": "Commands for mapping data to genomic features.",
        "commands": ["match_region_annotation", "match_position_annotation", "map_positions_to_genes", "graph_upstream_UCSC"]
    },
    "Plotting": {
        "file": "plotting.md",
        "desc": "Commands for generating visual plots and graphs.",
        "commands": ["all_plots", "volcano_plot", "manhattan_plot", "coverage_plot", "plot_methylation_curve", "graph_gene_regions", "graph_full_gene", "graph_upstream_gene_methylation"]
    }
}

def get_sphinx_data(func):
    doc = inspect.getdoc(func)
    if not doc: return {}, "No description available."
    main_desc = doc.split(':param')[0].strip()
    params = re.findall(r':param\s+(.*?):\s+(.*?)(?=\n\s*:(?:type|param|return|rtype)|$)', doc, re.DOTALL)
    param_dict = {name.strip(): desc.strip().replace('\n', ' ') for name, desc in params}
    return param_dict, main_desc

def build_cli_pages():
    os.makedirs("docs", exist_ok=True)
    all_methods = dict(inspect.getmembers(DiffMethylTools, predicate=inspect.isfunction))

    for group_name, info in GROUPS.items():
        # Start a new markdown page for this specific group
        md = [f"# {group_name} Commands\n", f"{info['desc']}\n", "---\n"]
        
        for func_name in info["commands"]:
            if func_name not in all_methods:
                continue
                
            func = all_methods[func_name]
            sphinx_params, description = get_sphinx_data(func)
            sig = inspect.signature(func)
            
            md.append(f"## `{func_name}`")
            md.append(f"{description}\n")
            md.append("| Argument | Default | Description |")
            md.append("| :--- | :--- | :--- |")
            
            for p_name, p_val in sig.parameters.items():
                if p_name == 'self': continue
                d = sphinx_params.get(p_name, "---")
                default = p_val.default if p_val.default != inspect.Parameter.empty else "*Required*"
                md.append(f"| `--{p_name}` | `{default}` | {d} |")
            md.append("\n---\n")

        # Save this specific group to its own file (e.g., docs/analysis.md)
        file_path = os.path.join("docs", info["file"])
        with open(file_path, "w") as f:
            f.write("\n".join(md))
        print(f"✅ Created {file_path}")

if __name__ == "__main__":
    build_cli_pages()
