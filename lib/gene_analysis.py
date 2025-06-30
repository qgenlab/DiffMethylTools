from lib import InputProcessor
import polars as pl
import pandas as pd
import numpy as np
import multiprocessing
from joblib import Parallel, delayed
from pathlib import Path

def assert_required_columns(df, required_columns):
        from lib import Analysis
        Analysis().assert_required_columns(df, required_columns)

def process_chromosome_group(chrom, group_filtered_table, df_chrom_lookup):
    result_rows = []
    import time
    start_chrom = time.time()

    if chrom not in df_chrom_lookup:
        return result_rows
    
    df = df_chrom_lookup[chrom]
    starts = df["start"].to_numpy()
    ends = df["end"].to_numpy()
    diffs = df["diff"].to_numpy()

    chrom_starts = group_filtered_table["chromStart"].to_numpy()
    genes = group_filtered_table["gene"].to_numpy()

    for start, end, diff in zip(starts, ends, diffs):
        mask = (chrom_starts >= start) & (chrom_starts <= end)
        overlapping_genes = genes[mask]
        unique_genes = np.unique(overlapping_genes)
        result_rows.extend({"chrom": chrom, "diff": diff, "gene": gene} for gene in unique_genes)

    return result_rows

def find_nearest_gene(ccre_name, max_distance_tag, bed_gene_strings, bed_chrom_starts, gene_starts, gene_ends, gene_genes):
    query_indices = np.where(bed_gene_strings == ccre_name)[0]

    if len(query_indices) == 0:
        return None

    query_chrom_starts = bed_chrom_starts[query_indices]

    # get all difference between chrom_starts and gene_starts/ends. There are multiple chrom_starts of different dimensions, so simple subtraction will not work
    start_diffs = np.abs(query_chrom_starts[:, None] - gene_starts).min(axis=1)
    end_diffs = np.abs(query_chrom_starts[:, None] - gene_ends).min(axis=1)

    min_start_index = np.argmin(start_diffs)
    min_end_index = np.argmin(end_diffs)

    min_start_diff = start_diffs[min_start_index]
    min_end_diff = end_diffs[min_end_index]

    if min(min_start_diff, min_end_diff) > max_distance_tag:
        return None
    elif min_start_diff <= min_end_diff:
        return gene_genes[min_start_index]
    else:
        return gene_genes[min_end_index]
    
def process_chunk(chunk, bed_gene_strings, bed_chrom_starts, gene_starts, gene_ends, gene_genes):
    results = []
    for ccre_name, max_distance_tag in chunk:
        result = find_nearest_gene(
            ccre_name, max_distance_tag.astype(int), bed_gene_strings, bed_chrom_starts, gene_starts, gene_ends, gene_genes
        )
        results.append(result)
    return results

def window_based_gene(positions: InputProcessor.data_container, gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], min_pos_diff=0, bed_file="", gtf_file="", enhd_thr = 500000, enhp_thr = 50000, prom_thr = 2000, processes=12):
    """ 
    
    Required Columns: ["chrom", "start", "end", "diff"]
    
    """

    assert_required_columns(positions, ["chrom", "start", "end", "diff"])
    
    df = pl.from_pandas(positions).sort(["chrom","start", "end"])

    bed = pl.read_csv(bed_file, separator="\t", columns=[0,1,4], new_columns=["chrom","chromStart","gene_string"],has_header=False)


    gene_locations = pl.read_csv(gtf_file, 
                        skip_rows=5, 
                        separator="\t", 
                        has_header=False,
                        columns= [0,2,3,4,6,8],
                        new_columns=["chrom", "gene_type", "start", "end", "strand", "gene_info"]
                    ) \
                    .filter((pl.col("gene_type") == "gene") & (pl.col("chrom").str.contains("chr"))) \
                    .select(pl.exclude(["gene_type"]))

    
    gene_locations = gene_locations.with_columns(pl.col("gene_info").str.split(";"))
    

    gene_locations = gene_locations.with_columns(gene_type=pl.col("gene_info").list.get(1)) \
            .with_columns(pl.col("gene_type").str.slice(12)) \
            .with_columns(pl.col("gene_type").str.head(-1))
    

    gene_locations = gene_locations.with_columns(gene=pl.col("gene_info").list.get(2)) \
            .with_columns(pl.col("gene").str.slice(12)) \
            .with_columns(pl.col("gene").str.head(-1))
    
    gene_locations = gene_locations.with_columns(pl.col("start")-1)

    result = pd.DataFrame() 
    ccre_result = pd.DataFrame()
    # change to extend/append dataframe by row
    for pos in gene_regions:
        print(pos)
        filtered_table = bed.filter(
            pl.col("gene_string").str.to_lowercase().str.contains(pos.lower())
        ).with_columns(
            # Split the rows by '/' to get value:key pairs
            pl.col("gene_string").str.split("/").alias("pairs")
        ).explode("pairs")\
        .filter(
            pl.col("pairs").str.to_lowercase().str.contains(pos.lower())
        ).with_columns(
            # Split each pair by ':' to separate values and keys
            pl.col("pairs").str.split(":").list.get(0 if pos != "CCRE" else 1).alias("gene")
        ).select(pl.col(["chrom","chromStart","gene"])).sort(["chrom","chromStart"])

        import time

        t = time.time()

        result_rows = []

        grouped = filtered_table.group_by("chrom")
        
        df_chrom_lookup = {i[0]:x for i,x in df.group_by("chrom")}

        results = Parallel(n_jobs=processes)(delayed(process_chromosome_group)(chrom[0], group_filtered_table, df_chrom_lookup) for chrom, group_filtered_table in grouped)
        
        result_rows = [item for sublist in results for item in sublist]
        loop_result = pd.DataFrame(result_rows)

        loop_result = loop_result[loop_result["diff"].abs() > min_pos_diff]
        genes = loop_result["gene"].value_counts().rename(pos)
        diffs = (loop_result.groupby(loop_result["gene"])["diff"].sum().div(genes)).rename(pos + "_diff")
        if pos == "CCRE":
            ccre_result = pd.concat([ccre_result, genes.to_frame(), diffs.to_frame()]).fillna(0)
        else:
            result = pd.concat([result, genes.to_frame(), diffs.to_frame()]).fillna(0)

    non_ccre_position = [x for x in gene_regions if x != "CCRE"]
    result = result.groupby(result.index).sum()
    result[non_ccre_position] = result[non_ccre_position].astype(int)

    # TODO test this in window-based. 
    if "CCRE" in gene_regions: 
        ccre_result = ccre_result.groupby(ccre_result.index).sum()
        
        new_tags = []
        for tag in ccre_result.index.str.split("@").str[1].tolist():
            new_tag = tag
            if "_" in tag:
                new_tag = tag.split("_")[0]

            if "enhD" in new_tag:
                new_tags.extend(["enhD"])
            elif "enhP" in new_tag:
                new_tags.extend(["enhP"])
            elif "prom" in new_tag:
                new_tags.extend(["prom"])
            else:
                new_tags.extend([None])

        ccre_result["tag"] = new_tags

        nearest_genes = []

        bed = bed.with_columns(
                pl.col("gene_string")
                .str.extract(r".*CCRE:(.+?)@.*") 
                .alias("gene_string")
        ).drop_nulls()

        # Extract ccre_name and tag as NumPy arrays
        ccre_names = ccre_result.index.str.split("@").str[0].to_numpy()
        tags = ccre_result["tag"].to_numpy()

        # Map tags to max_distance_tag using NumPy
        tag_to_threshold = np.vectorize(lambda tag: enhd_thr if tag == "enhD" else enhp_thr if tag == "enhP" else prom_thr if tag == "prom" else -1)
        max_distance_tags = tag_to_threshold(tags)

        # Extract gene_string from bed as a NumPy array
        bed_gene_strings = bed.select("gene_string").to_numpy().flatten()
        bed_chrom_starts = bed.select("chromStart").to_numpy().flatten()

        # Extract gene_locations columns as NumPy arrays
        gene_chroms = gene_locations.select("chrom").to_numpy().flatten()
        gene_starts = gene_locations.select("start").to_numpy().flatten()
        gene_ends = gene_locations.select("end").to_numpy().flatten()
        gene_genes = gene_locations.select("gene").to_numpy().flatten()

        # Prepare nearest_genes array
        nearest_genes = []

        chunks = np.array_split(list(zip(ccre_names, max_distance_tags)), processes)

        # Find the nearest gene for each ccre name, with the filter in mind using NumPy
        print("starting parallel")
        import time
        t = time.time()
        nearest_genes = Parallel(n_jobs=processes, backend='multiprocessing')(
            delayed(process_chunk)(chunk, bed_gene_strings, bed_chrom_starts, gene_starts, gene_ends, gene_genes)
            for chunk in chunks
        )

                # Flatten the results from all chunks
        nearest_genes = [gene for sublist in nearest_genes for gene in sublist]

        print("Time taken", time.time() - t)
                

        # for ccre_name, (_, _, tag) in ccre_result.iterrows():
        #     ccre_name = ccre_name[:ccre_name.find("@")]

        #     query_df = bed.filter(pl.col("gene_string") == ccre_name)

        #     result_df = query_df.join_asof(
        #         gene_locations,
        #         left_on="chromStart",
        #         right_on="start",
        #         by="chrom",
        #         strategy="nearest"
        #     )

        #     result_df = result_df.with_columns(start_diff=(pl.col("chromStart") - pl.col("start")).abs(), end_diff=(pl.col("chromStart") - pl.col("end")).abs())
        #     min_start = result_df.sort("start_diff").head(1)
        #     min_end = result_df.sort("end_diff").head(1)

        #     max_distance_tag = enhd_thr if tag == "enhD" else enhp_thr if tag == "enhP" else prom_thr if tag == "prom" else None

        #     if max_distance_tag is not None and min(min_start["start_diff"][0], min_end["end_diff"][0]) > max_distance_tag:
        #         closest_gene = None
        #     elif min_start["start_diff"][0] <= min_end["end_diff"][0]:
        #         closest_gene = min_start["gene"][0]
        #     else: 
        #         closest_gene = min_end["gene"][0]

        #     nearest_genes.extend([closest_gene])

        print("Time taken", time.time() - t)

        ccre_result["nearest_genes"] = nearest_genes
        ccre_result["CCRE"] = ccre_result["CCRE"].astype(int)

    print(ccre_result)
    return result, ccre_result
    
def position_based_gene(positions: InputProcessor.data_container, gene_regions: list[str]|str = ["intron", "exon", "upstream", "CCRE"], min_pos_diff=0, bed_file="", gtf_file=""):
    """
    
    Required Columns: ["chrom", "chromStart", "diff"]
    
    """
    def max_abs_windowed_diff(loop_result, name, window_size=100, pos_col="chromStart"):
        loop_result = loop_result.copy()
        loop_result["bin"] = (loop_result[pos_col] // window_size)
        binned_diff = (
            loop_result.groupby(["gene", "bin"])["diff"]
            .mean()
            .reset_index()
        )
        idx = binned_diff.groupby("gene")["diff"].apply(lambda x: x.abs().idxmax())
        max_abs_per_gene = binned_diff.loc[idx].set_index("gene")["diff"].rename(name)
        return max_abs_per_gene

    if gtf_file == "gencode.v41.chr_patch_hapl_scaff.annotation.gtf": gtf_file = Path(__file__).resolve().parent.parent / gtf_file
    if bed_file == "outfile_w_hm450.bed": bed_file = Path(__file__).resolve().parent.parent / bed_file

    assert_required_columns(positions, ["chrom", "chromStart", "diff"])
    

    df = positions
    
    if "strand" in df.columns:
        df.loc[df["strand"] == "-", "chromStart"] -= 1
        if "chromEnd" in df.columns:
            df.loc[df["strand"] == "-", "chromEnd"] -= 1
        df["strand"] = "+"

    df = pl.from_pandas(df)

    bed = pl.read_csv(bed_file, separator="\t", columns=[0,1,4], new_columns=["chrom","chromStart","gene_string"],has_header=False)

    gene_locations = pl.read_csv(gtf_file, 
                        skip_rows=5, 
                        separator="\t", 
                        has_header=False,
                        columns= [0,2,3,4,6,8],
                        new_columns=["chrom", "gene_type", "start", "end", "strand", "gene_info"]
                    ) \
                    .filter((pl.col("gene_type") == "gene") & (pl.col("chrom").str.contains("chr"))) \
                    .select(pl.exclude(["gene_type"]))

    
    gene_locations = gene_locations.with_columns(pl.col("gene_info").str.split(";"))
    

    gene_locations = gene_locations.with_columns(gene_type=pl.col("gene_info").list.get(1)) \
            .with_columns(pl.col("gene_type").str.slice(12)) \
            .with_columns(pl.col("gene_type").str.head(-1))
    

    gene_locations = gene_locations.with_columns(gene=pl.col("gene_info").list.get(2)) \
            .with_columns(pl.col("gene").str.slice(12)) \
            .with_columns(pl.col("gene").str.head(-1))

    # gene_strand_map = gene_locations.select(["gene", "strand"]).to_pandas().drop_duplicates() ###################
    
    gene_locations = gene_locations.with_columns(pl.col("start")-1)
    
    print(bed)
    # df_ = pa.Table.from_pandas(df)
    result = pd.DataFrame() 
    ccre_result = pd.DataFrame()
    # change to extend/append dataframe by row
    for pos in gene_regions:
        filtered_table = bed.filter(
            pl.col("gene_string").str.to_lowercase().str.contains(pos.lower()) ###############################
        ).with_columns(
        #  Split the rows by '/' to get value:key pairs
            pl.col("gene_string").str.split("/").alias("pairs")
        ).explode("pairs")\
        .filter(
            pl.col("pairs").str.to_lowercase().str.contains(pos.lower())
        ).with_columns(
        #    # Split each pair by ':' to separate values and keys
            pl.col("pairs").str.split(":").list.get(0 if pos != "CCRE" else 1).alias("gene")
        ).select(pl.col(["chrom","chromStart","gene"])) #######################################
        # print(pos)
        # print(df)
        # print(pos, filtered_table)
        loop_result = df.join(filtered_table, how="inner", on=["chrom","chromStart"])
        # print(loop_result)
        loop_result = loop_result.to_pandas()
        # print(pos)
        # print(loop_result)
        loop_result = loop_result[loop_result["diff"].abs() > min_pos_diff]
        # print(loop_result)
        genes = loop_result["gene"].value_counts().rename(pos)
        # diffs = (loop_result.groupby(loop_result["gene"])["diff"].sum().div(genes)).rename(pos + "_diff")
        diffs = max_abs_windowed_diff(loop_result, pos + "_diff")
        if pos == "CCRE":
            ccre_result = pd.concat([ccre_result, genes.to_frame(), diffs.to_frame()]).fillna(0)
        else:
            result = pd.concat([result, genes.to_frame(), diffs.to_frame()]).fillna(0)

    non_ccre_position = [x for x in gene_regions if x != "CCRE"]
    result = result.groupby(result.index).sum()
    result[non_ccre_position] = result[non_ccre_position].astype(int)
    
    if "CCRE" in gene_regions: 
        ccre_result = ccre_result.groupby(ccre_result.index).sum()
        
        new_tags = []
        for tag in ccre_result.index.str.split("@").str[1].tolist():
            new_tag = tag
            if "_" in tag:
                new_tag = tag.split("_")[0]

            if "enh" in new_tag:
                new_tags.extend(["enhancer"])
            elif "prom" in new_tag:
                new_tags.extend(["promoter"])
            else:
                new_tags.extend([new_tag])

        ccre_result["tag"] = new_tags

        nearest_genes = []

        bed = bed.with_columns(
                pl.col("gene_string")
                .str.extract(r".*CCRE:(.+?)@.*") 
                .alias("gene_string")
        ).drop_nulls()

        # TODO multithread this

        for ccre_name, _ in ccre_result.iterrows():
            ccre_name = ccre_name[:ccre_name.find("@")]

            query_df = bed.filter(pl.col("gene_string") == ccre_name)

            result_df = query_df.join_asof(
                gene_locations,
                left_on="chromStart",
                right_on="start",
                by="chrom",
                strategy="nearest"
            )

            result_df = result_df.with_columns(start_diff=(pl.col("chromStart") - pl.col("start")).abs(), end_diff=(pl.col("chromStart") - pl.col("end")).abs())
            min_start = result_df.sort("start_diff").head(1)
            min_end = result_df.sort("end_diff").head(1)

            if min_start["start_diff"][0] <= min_end["end_diff"][0]:
                closest_gene = min_start["gene"][0]
            else: 
                closest_gene = min_end["gene"][0]

            nearest_genes.extend([closest_gene])

        ccre_result["nearest_genes"] = nearest_genes
        ccre_result["CCRE"] = ccre_result["CCRE"].astype(int)

    # result = result.merge(gene_strand_map, left_index=True, right_on="gene", how="left") ###############################
    # if "CCRE" in gene_regions: ##############################
    #    ccre_result = ccre_result.merge(gene_strand_map, left_on="nearest_genes", right_on="gene", how="left") ###############################

    return result, ccre_result
def match_region_annotation(regions_df: InputProcessor.data_container, bed_file:str = "CpG_gencodev42ccrenb_repeat_epic1v2hm450.bed"):
    """
    Reads region file and annotation file, finds matches, and counts occurrences.
    """
    if bed_file == "CpG_gencodev42ccrenb_repeat_epic1v2hm450.bed": bed_file = Path(__file__).resolve().parent.parent / bed_file
    # regions_df = pd.read_csv(regions, sep=',', header=0, names=['chr', 'start', 'end'], usecols=[0,1,2])
    annotation_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'id', 'annotation'])
    total_counts_1 = defaultdict(int)
    total_counts_2 = defaultdict(int)
    total_counts_3 = defaultdict(int)
    total_counts_4 = defaultdict(int)
    total_counts_p4 = defaultdict(int)
    total_gene = {}
    for _, region in regions_df.iterrows():
        matched_annotations = annotation_df[
            (annotation_df['chr'] == region['chromosome']) &
            (annotation_df['start'] <= region['end']) &
            (annotation_df['end'] >= region['start'])
        ]
        this_counts_1 = defaultdict(int)
        this_counts_2 = defaultdict(int)
        en_list = defaultdict(int)
        this_counts_3 = defaultdict(int)
        this_counts_4 = defaultdict(int)
        gend_counts = defaultdict(int)
        gend_counts_diff = defaultdict(int)
        for _, annotation in matched_annotations.iterrows():
            categories, encode_types, repeat, cpg_epic, gene_list = parse_annotation(annotation['annotation'])
            for _ge in gene_list:
               for _de in gene_list[ _ge ]:
                  if (_ge + ':'+_de) not in gend_counts: gend_counts[  _ge + ':'+_de ] = 1; gend_counts_diff[  _ge + ':'+_de + "_diff"] = [annotation["diff"]]
                  else: gend_counts[  _ge + ':'+_de ] += 1; gend_counts_diff[  _ge + ':'+_de + "_diff"].append(annotation["diff"])
            for cat in categories:
                this_counts_1[cat] += 1
                this_counts_2['Gene_Body' if 'Gene_' in cat else cat] += 1
            if len(repeat)==0: this_counts_3['No_repeat'] += 1;
            else:
                for rep_c in repeat:
                    this_counts_3[rep_c] += 1
            if len(cpg_epic)==0:
                this_counts_4['non-EPIC'] += 1;
            else:
                this_counts_4['EPIC'] += 1;
            for encode_type in encode_types:
                this_counts_1[encode_type] += 1
                this_counts_2[encode_type] += 1
                en_list[encode_type] += 1
        de_gede = []
        for _gede in gend_counts:
           if _gede[-3:]=='_nb' and _gede[:-3] in gend_counts:
              gend_counts[ _gede[:-3] ] += gend_counts[ _gede ]
              gend_counts_diff[_gede[:-3]+"_diff"].extend(gend_counts_diff[ _gede ])
              de_gede.append( _gede)
        for _gede in de_gede:
           del gend_counts[ _gede ]
           del gend_counts_diff[ _gede ]
        if 'Gene_Intron' in this_counts_1 or 'Gene_Exon' in this_counts_1 or len(en_list)>0:
           if 'IG' in this_counts_1: del this_counts_1['IG']
           if 'IG' in this_counts_2: del this_counts_2['IG']
        if 'enhD' in en_list or 'enhP' in en_list or 'prom' in en_list or 'K4m3' in en_list or 'enh' in en_list:
           if 'Gene_Intron' in this_counts_1: del this_counts_1['Gene_Intron']
        if 'No_repeat' in this_counts_3 and this_counts_3['No_repeat']< matched_annotations.shape[0] - this_counts_3['No_repeat']:
           del this_counts_3['No_repeat']
        if 'Gene_Exon' in this_counts_1:
           print("Exon", region)
        en_dis = []
        en_cod_type = list(this_counts_1.keys())
        for _c in en_cod_type:
           if '_' in _c and _c.split('_')[0] in ['enhD', 'enhP', 'prom', 'K4m3', 'enh'] and _c.split('_')[0] in en_cod_type:
              en_dis.append( _c )
        for _c in en_dis:
           if _c in this_counts_1: del this_counts_1[_c]
           if _c in this_counts_2: del this_counts_2[_c]
        for _c in this_counts_1:
           total_counts_1[ _c ] += 1
        for _c in this_counts_2:
           total_counts_2[ _c ] += 1
        for _c in this_counts_3:
           total_counts_3[ _c ] += 1
        for _c in this_counts_4:
           total_counts_4[ _c ] += 1
        for _c in this_counts_4:
           total_counts_p4[ _c ] += this_counts_4[ _c ]
    cols = ['CDS', 'CTCF', 'CTCF_nb', 'K4m3', 'K4m3_nb', 'enhD', 'enhD_nb', 'enhP', 'enhP_nb', 'exon', 'intron', 'prom', 'prom_nb']
    d = {}
    for e, val in gend_counts.items():
        k = e.split(":")
        diff_l = gend_counts_diff[e+"_diff"]
        if len(k) != 2:
            continue
        if k[0] in cols:
            col, gene = k[0], k[1]
        else:
            gene, col = k[0], k[1]
        if col not in cols:
            continue
        if gene not in d:
            d[gene] = {}
        d[gene][col] = d[gene].get(col, 0) + val
        d[gene][col+"_diff"] = sum(diff_l)/len(diff_l)
    return total_counts_1, total_counts_2, total_counts_3, total_counts_4, total_counts_p4, d
