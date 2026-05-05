# Annotation Commands

Commands for mapping data to genomic features.

---

## `match_region_annotation`
Intersect identified regions (DMRs) with genomic annotations (e.g., GENCODE) to determine 
genomic context (promoters, exons, introns, etc.).

.. note::
    Required columns in ``regions_df``: ``['chrom', 'chromStart', 'chromEnd']``.
    If ``pipeline`` is enabled and ``regions_df`` is None, it defaults to the 'cluster_df' 
    from the ``generate_DMR`` step.

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--regions_df` | `None` | Input dataframe containing genomic regions to annotate. |
| `--ref_folder` | `None` | Path to the reference genome folder. |
| `--bed_file` | `CpG_gencode_annotation.bed` | Filename of the annotation BED file located in the ref_folder, defaults to "CpG_gencode_annotation.bed" |
| `--name` | `match_region_annotation` | Prefix for output files, defaults to "match_region_annotation" |
| `--annotation_or_region` | `region` | Determines the perspective of the output overlap ("annotation" or "region"), defaults to "region" |
| `--show_counts` | `False` | If True, returns count statistics of the annotations, defaults to False |

---

## `match_position_annotation`
Intersect single genomic positions (DMLs) with genomic annotations to determine local context.

.. note::
    Required columns in ``regions_df``: ``['chrom', 'chromStart']``.
    Similar to ``match_region_annotation``, but optimized for single-base coordinate matching
    rather than range-based intersection.

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--regions_df` | `None` | Input dataframe containing genomic positions to annotate. |
| `--ref_folder` | `None` | Path to the reference genome folder. |
| `--bed_file` | `CpG_gencode_annotation.bed` | Filename of the annotation BED file, defaults to "CpG_gencode_annotation.bed" |
| `--name` | `match_position_annotation` | Prefix for output files, defaults to "match_position_annotation" |

---

## `map_positions_to_genes`
Map positions to genes.

.. note::
    Required columns for ``positions``:
        - ``["chromosome", "position_start", "diff"]``

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--positions` | `None` | Position data. Not necessary if the pipeline is in use, defaults to None |
| `--ref_folder` | `None` | --- |
| `--gene_regions` | `['intron', 'exon', 'upstream', 'CCRE']` | Gene regions to map to. Options are any combination of ``["intron", "exon", "upstream", "CCRE"]``, defaults to ``["intron", "exon", "upstream", "CCRE"]`` |
| `--min_pos_diff` | `0` | Minimum position difference for mapping, defaults to 0 |
| `--gtf_file` | `gencode.chr_patch_hapl_scaff.annotation.gtf` | GTF annotation file with unflexible input format, defaults to "gencode.v42.chr_patch_hapl_scaff.annotation.gtf" |
| `--bed_file` | `CpG_gencode_annotation.bed` | BED annotation file with unflexible input format, defaults to "CpG_gencode_annotation.bed" |
| `--pipeline_input_source` | `auto` | Pipeline input source for pipelining, options are ["auto", "map_win_2_pos", "generate_q_values", "filters"], defaults to "auto" |
| `--rerun` | `False` | Rerun the analysis. If False, load previous output. Defaults to False. |

---

## `graph_upstream_UCSC`
Generate a UCSC graph of upstream gene methylation.

.. note::
    Required columns for ``position_data``:
        - ``["chromosome", "position_start", "methylation_percentage*"]``
        - There must be a ``methylation_percentage`` column for each sample to plot.

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--gene_name` | `*Required*` | Gene name |
| `--position_data` | `None` | Position data |
| `--ref_folder` | `None` | --- |
| `--name` | `UCSC_graph.bedGraph` | Output BEDGraph file name, defaults to "UCSC_graph.bedGraph" |
| `--before_tss` | `5000` | Distance before transcrption start site (TSS), defaults to 5000 |
| `--gtf_file` | `gencode.chr_patch_hapl_scaff.annotation.gtf` | GTF file, defaults to "gencode.chr_patch_hapl_scaff.annotation.gtf" |

---
