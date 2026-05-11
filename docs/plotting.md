# Plotting Commands

Commands for generating visual plots and graphs.

---

## `all_plots`
Execute a suite of visualization methods including Volcano plots, Manhattan plots, 
and gene-centric methylation graphs.

.. note::
    This method forces ``self.pipeline = False`` to ensure that all plots are generated 
    using the specific data objects provided as arguments rather than cached state.

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--data` | `*Required*` | Position-level methylation data (e.g., DMLs). |
| `--ref_folder` | `*Required*` | Path to the reference genome folder containing necessary genomic annotations. |
| `--window_data` | `*Required*` | Window-based methylation data (e.g., DMRs). |
| `--gene` | `*Required*` | Gene annotation data. |
| `--ccre` | `*Required*` | Candidate Cis-Regulatory Elements (cCRE) data. |

---

## `volcano_plot`
Generate a volcano plot.

.. note::
    Required columns for ``data``:
        - ``["q-value", "diff"]``

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--data` | `None` | Input data. Not necessary if the pipeline is in use, defaults to None |
| `--name` | `plots/volcano_plot.png` | Output file name, defaults to "volcano_plot.png" |
| `--threshold` | `0.05` | Q-value threshold for horizontal line. Set to None to have no line, defaults to 0.05 |
| `--line` | `None` | Vertical line threshold for the `abs(line)` vertical line. Set to None to have no lines, defaults to None |
| `--x_range` | `(-1, 1)` | X-axis range, defaults to (-1, 1) |
| `--y_max` | `None` | Y-axis maximum, defaults to None |
| `--title` | `None` | Plot title, defaults to None for a generic title |
| `--x_label` | `None` | X-axis label, defaults to None for a generic label |
| `--y_label` | `None` | Y-axis label, defaults to None for a generic label |
| `--position_or_window` | `auto` | The position-based or window-based results to use as input if DiffMethylTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto" |

---

## `manhattan_plot`
Generate a manhattan plot.

.. note::
    Required columns for ``data`` are in one of the following column formats:
        - ``["chromosome", "q-value", "position_start"]``
        - ``["chromosome", "q-value", "region_start"]``

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--data` | `None` | Input data. Not necessary if the pipeline is in use, defaults to None |
| `--name` | `plots/manhattan_plot.png` | Output file name, defaults to "manhattan_plot.png" |
| `--threshold` | `0.05` | Q-value threshold for horizontal line. Set to None to have no line, defaults to 0.05 |
| `--title` | `None` | Plot title, defaults to None for a generic title |
| `--x_label` | `None` | X-axis label, defaults to None for a generic label |
| `--y_label` | `None` | Y-axis label, defaults to None for a generic label |
| `--position_or_window` | `auto` | The position-based or window-based results to use as input if DiffMethylTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto" |

---

## `coverage_plot`
Generate a coverage plot.

.. note::
    ``case_data`` and ``ctr_data`` must be from a list of samples, where each contains one of the following column formats:
        - ``["coverage"]``
        - ``["positive_methylation_count", "negative_methylation_count"]``

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--case_data` | `*Required*` | The case data to plot. |
| `--ctr_data` | `*Required*` | The control data to plot. |
| `--name` | `plots/coverage_plot.png` | Output file name, defaults to "coverage_plot.png" |
| `--cov_min` | `1` | Minimum coverage display, defaults to 1 |
| `--cov_max` | `-1` | Maximum coverage display, defaults to -1 |
| `--cov_max_percentile` | `99.5` | Maximum coverage percentile display. Ranges from 0.0-100.0, defaults to 99.5. Overrides ``cov_max`` if set. |
| `--bins` | `20` | Number of bins, defaults to 20 |
| `--title` | `None` | Plot title, defaults to None for a generic title |
| `--x_label` | `None` | X-axis label, defaults to None for a generic label |
| `--y_label` | `None` | Y-axis label, defaults to None for a generic label |

---

## `plot_methylation_curve`
Generate plots showing methylation curves across specific genomic regions, including annotations for repeats and enhancers/promoters.

.. note::
    If ``pipeline`` mode is disabled, both ``region_data`` and ``position_data`` must be provided manually.
    If ``pipeline`` mode is enabled, the function automatically retrieves results from ``generate_DMR`` and ``filters``.

.. note::
    This method utilizes a sliding window approach (defined by ``window_size`` and ``step_size``) to smooth methylation values across regions.

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--region_data` | `None` | DMR or cluster data. If None and pipelined, uses 'cluster_df' from generate_DMR. |
| `--ref_folder` | `None` | Path to the reference genome folder containing annotation files. |
| `--position_data` | `None` | All position-level data. If None and pipelined, uses output from filters. |
| `--name` | `plots/.` | Directory path or prefix where the resulting plots will be saved, defaults to "." |
| `--repeat_regions_df` | `rmsk.txt` | Filename for repeat regions annotation (e.g., RepeatMasker), defaults to "rmsk.txt" |
| `--enhancer_promoter_df` | `encodeCcreCombined.bed` | Filename for enhancer/promoter annotation, defaults to "encodeCcreCombined.bed" |
| `--repeat_regions_columns` | `[5, 6, 7, 11]` | List of column indices to extract from the repeat regions file, defaults to [5,6,7,11] |
| `--enhancer_promoter_columns` | `[0, 1, 2, 12, 13]` | List of column indices to extract from the enhancer/promoter file, defaults to [0,1,2,12,13] |
| `--window_size` | `50` | Size of the sliding window in base pairs for smoothing, defaults to 50 |
| `--step_size` | `25` | Step size for the sliding window in base pairs, defaults to 25 |
| `--chr_filter` | `None` | Optional chromosome name to restrict plotting to a specific chromosome. |
| `--start_filter` | `None` | Optional genomic start coordinate to filter regions. |
| `--end_filter` | `None` | Optional genomic end coordinate to filter regions. |
| `--sample_start_ind` | `3` | Column index where individual sample methylation data begins in the input matrix, defaults to 3 |

---

## `graph_gene_regions`
Generate a graph of gene regions.

.. note::
    Required columns for ``gene_data``:
        - If ``"intron"`` is in ``gene_regions``: ``["intron", "intron_diff"]``
        - If ``"exon"`` is in ``gene_regions``: ``["exon", "exon_diff"]``
        - If ``"upstream"`` is in ``gene_regions``: ``["upstream", "upstream_diff"]``
    Required columns for ``ccre_data``:
        - If ``"CCRE"`` is in ``gene_regions``: ``["CCRE", "CCRE_diff"]``

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--gene_data` | `None` | Gene data. Not necessary if the pipeline is in use, defaults to None |
| `--ccre_data` | `None` | CCRE data. Not necessary if the pipeline is in use, defaults to None |
| `--name` | `plots/gene_regions.png` | Output file name, defaults to "gene_regions.png" |
| `--gene_regions` | `['intron', 'exon', 'upstream', 'CCRE']` | Gene regions to map to. Options are any combination of ``["intron", "exon", "upstream", "CCRE"]``, defaults to ``["intron", "exon", "upstream", "CCRE"]`` |
| `--intron_cutoff` | `-1` | Intron count vertical display cutoff, defaults to -1 (for no cutoff) |
| `--exon_cutoff` | `-1` | Exon count vertical display cutoff, defaults to -1 (for no cutoff) |
| `--upstream_cutoff` | `-1` | Upstream count vertical display cutoff, defaults to -1 (for no cutoff) |
| `--CCRE_cutoff` | `-1` | CCRE count vertical display cutoff, defaults to -1 (for no cutoff) |
| `--prom_cutoff` | `-1` | --- |
| `--title` | `None` | Plot title, defaults to None for a generic title |
| `--x_label` | `None` | X-axis label, defaults to None for a generic label |
| `--intron_y_label` | `None` | Intron Y-axis label, defaults to None for a generic label |
| `--exon_y_label` | `None` | Exon Y-axis label, defaults to None for a generic label |
| `--upstream_y_label` | `None` | Upstream Y-axis label, defaults to None for a generic label |
| `--CCRE_y_label` | `None` | CCRE Y-axis label, defaults to None for a generic label |
| `--prom_y_label` | `None` | --- |
| `--position_or_window` | `auto` | The position-based or window-based results to use as input if DiffMethylTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto" |

---

## `graph_full_gene`
Generate a graph of full gene methylation.

.. note::
    Required columns for ``position_data``:
        - ``["chromosome", "position_start", "methylation_percentage*"]``
        - There must be a ``methylation_percentage`` column for each sample to plot.

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--gene_name` | `*Required*` | The name of the gene to graph. |
| `--position_data` | `None` | Position data. Not necessary if the pipeline is in use, defaults to None |
| `--ref_folder` | `None` | --- |
| `--name` | `plots/gene_methylation_graph.png` | Output PNG file name, defaults to "gene_methylation_graph.png" |
| `--before_tss` | `0` | Distance before transcription start site, defaults to 0 |
| `--after_tss` | `None` | Distance after transcription start site, defaults to None |
| `--bin_size` | `500` | Bin size, defaults to 500 |
| `--start_marker` | `True` | Include start marker, defaults to True |
| `--end_marker` | `True` | Include end marker, defaults to True |
| `--deviation_display` | `True` | Display standard deviation regions, defaults to True |
| `--aggregate_samples` | `True` | Aggregate samples and display mean, defaults to True |
| `--legend_size` | `12` | Legend size, defaults to 12 |
| `--title` | `None` | Plot title, defaults to None for a generic title |
| `--x_label` | `None` | X-axis label, defaults to None for a generic label |
| `--y_label` | `None` | Y-axis label, defaults to None for a generic label |
| `--case_name` | `Case` | Case name in the legend, defaults to "Case" |
| `--ctr_name` | `Control` | Control name in the legend, defaults to "Control" |
| `--gtf_file` | `gencode.chr_patch_hapl_scaff.annotation.gtf` | GTF file, defaults to "gencode.chr_patch_hapl_scaff.annotation.gtf" |

---

## `graph_upstream_gene_methylation`
Generate a graph of upstream gene methylation.

.. note::
    Required columns for ``position_data``:
        - ``["chromosome", "position_start", "diff"]``
    Required columns for ``gene_data``:
        - If ``"intron"`` is in ``gene_regions``: ``["intron", "intron_diff"]``
        - If ``"exon"`` is in ``gene_regions``: ``["exon", "exon_diff"]``
        - If ``"upstream"`` is in ``gene_regions``: ``["upstream", "upstream_diff"]``

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--position_data` | `None` | Position data. Not necessary if the pipeline is in use, defaults to None |
| `--ref_folder` | `None` | --- |
| `--region_data` | `None` | --- |
| `--name` | `plots/upstream_methylation.png` | Output PNG file name, defaults to "upstream_methylation.png" |
| `--csv_name` | `upstream_methylation.csv` | Output CSV file name, defaults to "upstream_methylation.csv" |
| `--csv` | `None` | Input CSV file, defaults to None |
| `--left_distance` | `1000` | Left distance, defaults to 5000 |
| `--right_distance` | `1000` | Right distance, defaults to 5000 |
| `--window_size` | `100` | Window size, defaults to 100 |
| `--hypermethylated` | `True` | Include hypermethylated regions, defaults to True |
| `--gene_hypermethylated_min` | `20` | Minimum hypermethylation for genes, defaults to 0.20 |
| `--window_hypermethylated_min` | `5` | Minimum hypermethylation for windows, defaults to 0.05 |
| `--min_hypermethylated_windows` | `5` | Minimum hypermethylated windows, defaults to 5 |
| `--hypomethylated` | `True` | Include hypomethylated regions, defaults to True |
| `--gene_hypomethylated_max` | `-20` | Maximum hypomethylation for genes, defaults to -0.20 |
| `--window_hypomethylated_max` | `-5` | Maximum hypomethylation for windows, defaults to -0.5 |
| `--min_hypomethylated_windows` | `5` | Minimum hypomethylated windows, defaults to 5 |
| `--position_count` | `5` | Position count, defaults to 5 |
| `--clamp_positive` | `50` | Positive clamp value, defaults to 0.50 |
| `--clamp_negative` | `-50` | Negative clamp value, defaults to -0.50 |
| `--title` | `None` | Plot title, defaults to None for a generic title |
| `--gtf_file` | `gencode.chr_patch_hapl_scaff.annotation.gtf` | GTF file, defaults to "gencode.chr_patch_hapl_scaff.annotation.gtf" |
| `--position_or_window` | `auto` | Position or window, options are ["auto", "position", "window"], defaults to "auto" |
| `--position_or_region` | `region` | --- |

---
