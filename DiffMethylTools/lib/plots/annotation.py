import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import Optional, Union
from pathlib import Path
from collections import defaultdict
from ..analysis import Analysis
from collections import defaultdict, Counter


from ..input_processor import InputProcessor


class AnnotationMixin():
    def match_region_annotation(self, regions_df: InputProcessor.data_container, bed_file: InputProcessor.data_container, name: str = "match_region_annotation", annotation_or_region: str = "region", show_counts: bool = False):
        """
        Reads region file and annotation file, finds matches, and counts occurrences.
        """
        print(bed_file)
        if bed_file == "CpG_gencodev42ccrenb_repeat_epic1v2hm450.bed": bed_file = Path(__file__).resolve().parent.parent / bed_file
        annotation_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'id', 'annotation'])
        total_counts_1 = defaultdict(int)
        total_counts_2 = defaultdict(int)
        total_counts_3 = defaultdict(int)
        total_counts_4 = defaultdict(int)
        total_counts_p4 = defaultdict(int)
        total_gene = {}
        special_types = ['enhD', 'enhP', 'prom', 'K4m3', 'enh']
        add_to_name = ""
        # print(regions_df)
        chroms = regions_df['chrom'].unique()
        # print(chroms)
        for chrom in list(chroms):
            print(chrom)
            regions_df_chr = regions_df[regions_df["chrom"] == chrom]
            annotation_df_chr = annotation_df[annotation_df['chr'] == chrom]
            for _, region in regions_df_chr.iterrows():
                matched_annotations = annotation_df_chr[
                (annotation_df_chr['start'] <= region['end']) &
                (annotation_df_chr['end'] >= region['start'])
                ]
                this_counts_1 = defaultdict(int)
                this_counts_2 = defaultdict(int)
                en_list = defaultdict(int)
                this_counts_3 = defaultdict(int)
                this_counts_4 = defaultdict(int)
                gend_counts = defaultdict(int)
                if annotation_or_region == "annotation":
                    add_to_name = "annotation_based"
                    all_categories = []
                    all_encode_types = []
                    all_repeat = []
                    all_cpg_epic = []
                    all_gene_list = {}
                    for _, annotation in matched_annotations.iterrows():
                        try:
                            categories, encode_types, repeat, cpg_epic, gene_list = Analysis.parse_annotation(annotation['annotation'])
                        except AttributeError:
                            continue
                        all_categories.extend(categories)
                        all_encode_types.extend(encode_types)
                        all_repeat.extend(repeat)
                        # print(all_repeat)
                        all_cpg_epic.extend(cpg_epic)
                        for g in gene_list:
                            if g not in all_gene_list:
                                all_gene_list[g] = set(gene_list[g])
                            else:
                            #for k in gene_list[g]:
                                all_gene_list[g] |= set(gene_list[g])
                                # all_gene_list[g][k] = all_gene_list[g].get(k, 0) + gene_list[g][k]
                    for cat in set(all_categories):
                        if cat in special_types:
                            this_counts_1[cat] += 1
                            key = 'Gene_Body' if 'Gene_' in cat else cat
                            this_counts_2[key] += 1
                        else:
                            count_cat = all_categories.count(cat)
                            this_counts_1[cat] += 1 ############################################# count_cat
                            key = 'Gene_Body' if 'Gene_' in cat else cat
                            this_counts_2[key] += 1 ############################################### count_cat
                    if len(all_repeat) == 0:
                        this_counts_3['No_repeat'] += 1
                    else:
                        for rep_c in set(all_repeat):
                            this_counts_3[rep_c] += 1 ################################################all_repeat.count(rep_c)
                    if len(all_cpg_epic) == 0:
                        this_counts_4['non-EPIC'] += 1
                    else:
                        this_counts_4['EPIC'] += 1
                    for encode_type in set(all_encode_types):
                        this_counts_1[encode_type] += 1
                        this_counts_2[encode_type] += 1
                        en_list[encode_type] += 1
                    for _ge in all_gene_list:
                        for _de in all_gene_list[_ge]:
                            key = _ge + ':' + _de
                            gend_counts[key] += 1
                        #if key not in gend_counts:
                        #gend_counts[key] = all_gene_list[_ge][_de]
                        # else:
                        #    gend_counts[key] += all_gene_list[_ge][_de]
                else:
                    add_to_name = "region_based"
                    for _, annotation in matched_annotations.iterrows():
                        try:
                            categories, encode_types, repeat, cpg_epic, gene_list = Analysis.parse_annotation(annotation['annotation'])
                        except AttributeError:
                            continue
                        for _ge in gene_list:
                            for _de in gene_list[_ge]:
                                key = _ge + ':' + _de
                                gend_counts[key] = gend_counts.get(key, 0) + 1
                        for cat in categories:
                            this_counts_1[cat] += 1
                            this_counts_2['Gene_Body' if 'Gene_' in cat else cat] += 1
                            # print("s1 - 1: ", this_counts_1)
                            # print("s1 - 2: ", this_counts_2)
                        if len(repeat) == 0:
                            this_counts_3['No_repeat'] += 1
                        else:
                            for rep_c in repeat:
                                this_counts_3[rep_c] += 1
                        if len(cpg_epic) == 0:
                            this_counts_4['non-EPIC'] += 1
                        else:
                            this_counts_4['EPIC'] += 1
                        for encode_type in encode_types:
                            this_counts_1[encode_type] += 1
                            this_counts_2[encode_type] += 1
                            en_list[encode_type] += 1
                # print("s2 - 1: ", this_counts_1)
                # print("s2 - 2: ", this_counts_2)
                de_gede = []
                for _gede in gend_counts:
                    if _gede.endswith('_nb') and _gede[:-3] in gend_counts:
                        gend_counts[_gede[:-3]] += gend_counts[_gede]
                        de_gede.append(_gede)
                for _gede in de_gede:
                    del gend_counts[_gede]
                total_gene = Counter(total_gene) + Counter(dict(gend_counts))
                if 'Gene_Intron' in this_counts_1 or 'Gene_Exon' in this_counts_1 or len(en_list) > 0:
                    this_counts_1.pop('IG', None)
                    this_counts_2.pop('IG', None)
            # print("s3 - 1: ", this_counts_1)
            # print("s3 - 2: ", this_counts_2)
                if any(et in en_list for et in ['enhD', 'enhP', 'prom', 'K4m3', 'enh']):
                    # this_counts_1.pop('Gene_Intron', None)
                    this_counts_2['Gene_Body'] -= this_counts_1['Gene_Intron']
                    this_counts_1.pop('Gene_Intron', None)
            # print("s4 - 1: ", this_counts_1)
            # print("s4 - 2: ", this_counts_2)
                if annotation_or_region == "region" and 'No_repeat' in this_counts_3 and this_counts_3['No_repeat'] < matched_annotations.shape[0] - this_counts_3['No_repeat']:
                    this_counts_3.pop('No_repeat', None)
                en_dis = []
                en_cod_type = list(this_counts_1.keys())
                for _c in en_cod_type:
                    if '_' in _c and _c.split('_')[0] in special_types and _c.split('_')[0] in en_cod_type:
                        en_dis.append(_c)
                for _c in en_dis:
                    this_counts_1.pop(_c, None)
                    this_counts_2.pop(_c, None)
            # print("s5 - 1: ", this_counts_1)
            # print("s5 - 2: ", this_counts_2)
                for _c in this_counts_1:
                    total_counts_1[_c] += this_counts_1[_c]
                for _c in this_counts_2:
                    total_counts_2[_c] += this_counts_2[_c]
                for _c in this_counts_3:
                    total_counts_3[_c] += this_counts_3[_c]
                for _c in this_counts_4:
                    total_counts_4[_c] += this_counts_4[_c]
                for _c in this_counts_4:
                    total_counts_p4[_c] += this_counts_4[_c]
        self.__annotation_chart(total_counts_1, "Fraction of Occurrences: Gene Intron, Gene Exon, IG, and ENCODE Types", add_to_name + "_" + name+'_inex2.png')
        self.__annotation_chart(total_counts_2, "Fraction of Occurrences: Gene Body, IG, and ENCODE Types", add_to_name + "_" + name+'_body2.png')
        self.__annotation_chart(total_counts_3, "", add_to_name + "_" + name+'_repeat2.png')
        self.__annotation_chart(total_counts_4, "", add_to_name + "_" + name+'_epic2.png')
        self.__annotation_chart(total_counts_p4, "", add_to_name + "_" + name+'_epicP2.png')
        return [pd.DataFrame.from_dict(list(dict(total_counts_1).items())), pd.DataFrame.from_dict(list(dict(total_counts_2).items())), pd.DataFrame.from_dict(list(dict(total_counts_3).items())), pd.DataFrame.from_dict(list(dict(total_counts_4).items())), pd.DataFrame.from_dict(list(dict(total_counts_p4).items())), pd.DataFrame.from_dict(list(dict(total_gene).items()))]
    def match_position_annotation(self, regions_df: InputProcessor.data_container, bed_file: InputProcessor.data_container, name:str="match_position_annotation", show_counts = False):
        """
        Reads position file and annotation file, finds matches, and counts occurrences.
        """
        annotation_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'id', 'annotation'])
        total_counts_1 = defaultdict(int)
        total_counts_2 = defaultdict(int)
        total_counts_3 = defaultdict(int)
        total_counts_4 = defaultdict(int)
        total_counts_p4 = defaultdict(int)
        total_gene = {}
        # regions_df = pd.read_csv(regions_df)
        merged_df = pd.merge(annotation_df, regions_df, left_on=["chr", "start"], right_on=["chrom", "chromStart"], how="inner")
        this_counts_1 = defaultdict(int)
        this_counts_2 = defaultdict(int)
        en_list = defaultdict(int)
        this_counts_3 = defaultdict(int)
        this_counts_4 = defaultdict(int)
        gend_counts = defaultdict(int)
        merged_df = merged_df.dropna(subset=["annotation"])
        for _, annotation in merged_df.iterrows():
            if annotation['annotation'] == "nan" : continue
            # print(annotation['annotation'])
            categories, encode_types, repeat, cpg_epic, gene_list = Analysis.parse_annotation(annotation['annotation'])
            for _ge in gene_list:
               for _de in gene_list[ _ge ]:
                  if (_ge + ':'+_de) not in gend_counts: gend_counts[  _ge + ':'+_de ] = 1;
                  else: gend_counts[  _ge + ':'+_de ] += 1
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
              de_gede.append( _gede)
        for _gede in de_gede:
           del gend_counts[ _gede ]
        total_gene = Counter(total_gene) + Counter(dict(gend_counts))
        en_dis = []
        print(this_counts_1, this_counts_2, this_counts_3, this_counts_4) #, gend_counts)
        en_cod_type = list(this_counts_1.keys())
        for _c in en_cod_type:
           if '_' in _c and _c.split('_')[0] in ['enhD', 'enhP', 'prom', 'K4m3', 'enh'] and _c.split('_')[0] in en_cod_type:
              en_dis.append( _c )
        for _c in en_dis:
           if _c in this_counts_1: del this_counts_1[_c]
           if _c in this_counts_2: del this_counts_2[_c]
        for _c in this_counts_1:
           total_counts_1[ _c ] = 1 + total_counts_1[ _c ]
        for _c in this_counts_2:
           total_counts_2[ _c ] = 1 + total_counts_2[ _c ]
        for _c in this_counts_3:
           total_counts_3[ _c ] = 1 + total_counts_3[ _c ]
        for _c in this_counts_4:
           total_counts_4[ _c ] = 1 + total_counts_4[ _c ]
        for _c in this_counts_4:
           total_counts_p4[ _c ] = this_counts_4[ _c ] + total_counts_p4[ _c ]
        self.__annotation_chart(total_counts_1, "Fraction of Occurrences: Gene Intron, Gene Exon, IG, and ENCODE Types", name+'_inex2.png', nb = show_counts)
        self.__annotation_chart(total_counts_2, "Fraction of Occurrences: Gene Body, IG, and ENCODE Types", name+'_body2.png', nb = show_counts)
        self.__annotation_chart(total_counts_3, "", name+'_repeat2.png', nb = show_counts)
        self.__annotation_chart(total_counts_4, "", name+'_epic2.png', nb = show_counts)
        self.__annotation_chart(total_counts_p4, "", name+'_epicP2.png', nb = show_counts)
        return [pd.DataFrame.from_dict(list(dict(total_counts_1).items())), pd.DataFrame.from_dict(list(dict(total_counts_2).items())), pd.DataFrame.from_dict(list(dict(total_counts_3).items())), pd.DataFrame.from_dict(list(dict(total_counts_4).items())), pd.DataFrame.from_dict(list(dict(total_counts_p4).items())), pd.DataFrame.from_dict(list(dict(total_gene).items()))]
    def __annotation_chart(self, data, title, fig_name, keep_small=False, nb = False):
        """
        Plot and save a pie chart.
        """
        print(data)
        labels = sorted(list(data.keys()))
        sizes = [data[_lb_] for _lb_ in labels]
        total = sum(sizes)
        filtered_labels = []
        filtered_sizes = []
        var = 1 if not keep_small else 0
        for label, size in zip(labels, sizes):
            percent = 100 * size / total
            if percent >= var:
                filtered_labels.append(label)
                filtered_sizes.append(size)
        def make_autopct(nb):
            def autopct_func(pct):
                absolute = int(round(pct * total / 100.0))
                if nb:
                    return f'{pct:.1f}%\n({absolute})' if pct >= 1 else ''
                else:
                    return f'{pct:.1f}%' if pct >= 1 else ''
            return autopct_func
        colors = plt.get_cmap('Blues')(np.linspace(0.2, 0.7, len(filtered_sizes)))
        plt.figure(figsize=(7, 7))
        plt.pie(
            filtered_sizes,
            labels=filtered_labels,
            autopct=make_autopct(nb),
            startangle=140,
            colors=colors,
            textprops={'fontsize': 18},
            labeldistance=1.1
        )
        plt.title(title, fontsize=16)
        plt.savefig(fig_name, dpi=600, bbox_inches='tight')
        plt.close('all')
    def __match_to_gene(self, positions_df, bed_file, name: str = "match_region_annotation"):
        """
            Reads region file and annotation file, finds matches, and returns gene annotation matrix with average diff.
        """
        annotation_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'id', 'annotation'])
        # regions_df = pd.read_csv(regions_df)
        merged_df = pd.merge(annotation_df, positions_df, left_on=["chr", "start"], right_on=["chrom", "chromStart"], how="inner").dropna(subset=["annotation"])
        print(merged_df)
        gend_counts = defaultdict(int)
        gend_counts_diff = defaultdict(list)
        for _, annotation in merged_df.iterrows():
            categories, encode_types, repeat, cpg_epic, gene_list = Analysis.parse_annotation(annotation['annotation'])
            for _ge in gene_list:
                for _de in gene_list[_ge]:
                    key = f"{_ge}:{_de}"
                    gend_counts[key] += 1
                    gend_counts_diff[key + "_diff"].append(annotation["diff"])
        cols = {'CDS', 'CTCF', 'CTCF_nb', 'K4m3', 'K4m3_nb', 'enhD', 'enhD_nb', 'enhP', 'enhP_nb', 'exon', 'intron', 'prom', 'prom_nb'}
        d = {}
        for key, val in gend_counts.items():
            k = key.split(":")
            if len(k) != 2: continue
            col, gene = (k[1], k[0]) if k[0] not in cols else (k[0], k[1])
            if col not in cols: continue
            if gene not in d:
                d[gene] = {}
            d[gene][col] = d[gene].get(col, 0) + val
            d[gene][col + "_diff"] = sum(gend_counts_diff[key + "_diff"]) / len(gend_counts_diff[key + "_diff"])
        return pd.DataFrame.from_dict(d, orient='index').fillna(0)
