# Analysis Commands

Core commands for processing and analyzing methylation data.

---

## `all_analysis`
Run all analysis methods.

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

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--case_data` | `*Required*` | Case data (list of files) |
| `--ctr_data` | `*Required*` | Control data (list of files) |
| `--ref_folder` | `None` | --- |
| `--window_based` | `False` | Window-based analysis, defaults to True |
| `--min_cov_individual` | `10` | Minimum coverage filter (individual), defaults to 10 |
| `--min_cov_group` | `15` | Minimum coverage filter (group), defaults to 15 |
| `--filter_samples_ratio` | `0.6` | Minimum sample ratio filter. Used with min_cov_group, defaults to 0.6 |
| `--meth_group_threshold` | `0.2` | Methylation group threshold. Used with min_cov_group, defaults to 0.2 |
| `--cov_percentile` | `100.0` | Maximum coverage filter (percentile of sample coverage). Ranges from 0.0-100.0, defaults to 100.0 |
| `--min_samp_ctr` | `2` | Minimum samples in control, defaults to 2 |
| `--min_samp_case` | `2` | Minimum samples in case, defaults to 2 |
| `--max_q_value` | `0.05` | Maximum q-value filter, defaults to 0.05 |
| `--abs_min_diff` | `0.0` | Minimum absolute difference filter, defaults to 0.25 |
| `--features` | `None` | --- |

---

## `merge_tables`
Merge case and control data tables.

.. note::
    ``case_data`` and ``ctr_data`` must be from a list of samples, where each contains one of the following column formats:
        - ``["chromosome", "position_start", "coverage", "methylation_percentage"]``
        - ``["chromosome", "position_start", "positive_methylation_count", "negative_methylation_count"]``
        - Additionally, ``"position_end"`` and ``"strand"`` columns can be used. These will be considered during the join process.

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--case_data` | `*Required*` | The case data to be merged. |
| `--ctr_data` | `*Required*` | The control data to be merged. |
| `--min_cov_individual` | `10` | Minimum coverage filter (individual), defaults to 10 |
| `--min_cov_group` | `15` | Minimum coverage filter (group), defaults to 15 |
| `--filter_samples_ratio` | `0.6` | Minimum sample ratio filter. Used with min_cov_group, defaults to 0.6 |
| `--meth_group_threshold` | `0.2` | Methylation group threshold. Used with min_cov_group, defaults to 0.2 |
| `--cov_percentile` | `100.0` | Maximum coverage filter (percentile of sample coverage). Ranges from 0.0-100.0, defaults to 100.0 |
| `--min_samp_ctr` | `2` | Minimum samples in control, defaults to 2 |
| `--min_samp_case` | `2` | Minimum samples in case, defaults to 2 |
| `--rerun` | `False` | Rerun the analysis. If False, load previous output. Defaults to False. |
| `--small_mean` | `1` | --- |

---

## `position_based`
Perform position-based DML detection. Has options for using the gamma function, or the limma R package.

.. note::
    ``data`` must contain the following column format:
        - ``["chromosome", "position_start", "methylation_percentage*"]``
        - There must be a ``methylation_percentage`` column for each sample to test.

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--data` | `None` | Input data. Not necessary if the pipeline is in use, defaults to None |
| `--method` | `limma` | Position-based method to use, options are ``["gamma", "limma"]``, defaults to "limma" |
| `--features` | `None` | Features .csv file used with the "limma" method, has indicies with the names of samples in ``data`` and columns for the features, defaults to None |
| `--test_factor` | `Group` | Test factor used with the "limma" method, defaults to "Group" |
| `--processes` | `12` | Number of CPU processes, defaults to 12 |
| `--model` | `eBayes` | Limma model to use, options are ``["eBayes", "treat"]`` defaults to "eBayes" |
| `--min_std` | `0.1` | Minimum standard deviation filter, defaults to 0.1 |
| `--fill_na` | `True` | Fill NA values with row group average, defaults to True |
| `--rerun` | `False` | Rerun the analysis. If False, load previous output. Defaults to False. |

---

## `process_regions`
No description available.

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--region_file` | `None` | --- |
| `--ref_folder` | `None` | --- |
| `--annotation_file` | `CpG_gencode_annotation.bed` | --- |
| `--gene_bed_file` | `gencode.v42.chr_patch_hapl_scaff.annotation.genes.bed` | --- |
| `--ccre_file` | `encodeCcreCombined.bed` | --- |

---

## `generate_DMR`
Generate Differentially Methylated Regions (DMRs) by clustering DMLs.

.. note::
    Required columns for ``significant_position_data``:
        - ``["chromosome", "position_start", "diff"]``
    Required columns for ``position_data``:
        - ``["chromosome", "position_start", "diff"]``

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--significant_position_data` | `None` | Significant position data. Not necessary if the pipeline is in use, defaults to None |
| `--position_data` | `None` | All position data. Not necessary if the pipeline is in use, defaults to None |
| `--min_pos` | `3` | Minimum positions, defaults to 3 |
| `--neutral_change_limit` | `7.5` | Neutral change limit, defaults to 7.5 |
| `--neutral_perc` | `30` | Neutral percentage, defaults to 30 |
| `--opposite_perc` | `10` | Opposite percentage, defaults to 10 |
| `--significant_position_pipeline` | `auto` | The significant position-based or window-based results to use as input if DiffMethylTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto" |
| `--rerun` | `False` | Rerun the analysis. If False, load previous output. Defaults to False. |

---

## `filters`
Filter positions by q-value and minimum difference.

.. note::
    Required columns for ``data``:
        - ``["q-value", "diff"]``

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--data` | `None` | Input data. Not necessary if the pipeline is in use, defaults to None |
| `--max_q_value` | `0.05` | Maximum q-value filter, defaults to 0.05 |
| `--abs_min_diff` | `0.1` | Absolute minimum difference filter, defaults to 0.10 |
| `--position_or_window` | `auto` | The position-based or window-based results to use as input if DiffMethylTools is pipelined and no data is provided. Options are ``["auto", "position", "window"]``, defaults to "auto" |
| `--rerun` | `False` | Rerun the analysis. If False, load previous output. Defaults to False. |

---

## `map_win_2_pos`
Map windows to positions.

.. note::
    Required columns for ``window_data``:
        - ``["chromosome", "region_start", "region_end"]``

    Required columns for ``position_data``:
        - ``["chromosome", "position_start", "avg_case", "avg_ctr"]``

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--window_data` | `None` | Window data to map with. Not necessary if the pipeline is in use, defaults to None |
| `--position_data` | `None` | Position data to map the windows to. Not necessary if the pipeline is in use, defaults to None |
| `--processes` | `12` | Number of CPU processes, defaults to 12 |
| `--sub_window_size` | `100` | Sub-window size for a deeper difference filtering, defaults to 100 |
| `--sub_window_step` | `100` | Sub-window step size, defaults to 100 |
| `--sub_window_min_diff` | `0` | Sub-window minimum difference, defaults to 0 |
| `--pipeline_window_result` | `auto` | The function results to use as input if DiffMethylTools is pipelined and no data is provided. Options are ``["auto", "filters", "generate_q_values", "window_based"]``, defaults to "auto" |
| `--rerun` | `False` | Rerun the analysis. If False, load previous output. Defaults to False. |

---
