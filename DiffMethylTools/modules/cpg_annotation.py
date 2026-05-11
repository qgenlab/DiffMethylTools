import os,sys,yaml,shutil,subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm



class CpGAnnotator:
    def __init__(self,config_path):
        self.config=self._load_yaml(config_path)
        self.ref_name=self.config.get("ref_name")
        if not self.ref_name: raise ValueError("Process Rejected: Missing 'ref_name'")
        self.out_dir=Path(f"modules/{self.ref_name}")
        self.out_dir.mkdir(parents=True,exist_ok=True)
        self.gtf_info=self.config.get("gtf",{})
        self.cpg_info=self.config.get("CpG",{})
        if not self.gtf_info.get("path") or not os.path.exists(self.gtf_info.get("path")): raise FileNotFoundError("Missing GTF")
        if not self.cpg_info.get("path") or not os.path.exists(self.cpg_info.get("path")): raise FileNotFoundError("Missing CpG")
        self.ccre_info=self.config.get("ccre",{})
        self.rmsk_info=self.config.get("rmsk",{})
        self.epic_v2_info=self.config.get("epic_v2",{})
        self.epic_v1_info=self.config.get("epic_v1",{})
        self.hm450_info=self.config.get("hm450",{})
        self.updown_size=2000
        self.gtf_dict={}
    def _load_yaml(self,path):
        with open(path,'r') as file: return yaml.safe_load(file) or {}
    def _prepare_output_files(self):
        gtf_dest=self.out_dir/"gencode.chr_patch_hapl_scaff.annotation.gtf"
        if self.gtf_info.get("path") and os.path.exists(self.gtf_info.get("path")):
            if not gtf_dest.exists(): shutil.copy(self.gtf_info.get("path"),gtf_dest)
        else: open(gtf_dest,'a').close()
        ccre_dest=self.out_dir/"encodeCcreCombined.bed"
        if self.ccre_info.get("path") and os.path.exists(self.ccre_info.get("path")):
            if not ccre_dest.exists(): shutil.copy(self.ccre_info.get("path"),ccre_dest)
        else: pd.DataFrame(columns=range(15)).to_csv(ccre_dest,sep='\t',index=False,header=False)
        rmsk_dest=self.out_dir/"rmsk.txt"
        if self.rmsk_info.get("path") and os.path.exists(self.rmsk_info.get("path")):
            if not rmsk_dest.exists(): shutil.copy(self.rmsk_info.get("path"),rmsk_dest)
        else: pd.DataFrame(columns=['bin','swScore','milliDiv','milliDel','milliIns','chr','start','end','genoLeft','strand','repName','repeat_cat','repFamily','repStart','repEnd','repLeft','id']).to_csv(rmsk_dest,sep='\t',index=False,header=False)
    def _get_gene_exon_info(self,header_rown=5):
        gtf_df=pd.read_csv(self.gtf_info.get("path"),sep=self.gtf_info.get("sep","\t"),header=None,skiprows=header_rown)
        this_gene_info=[[],[]]
        col2_types={}
        for _,t_row in tqdm(gtf_df.iterrows(),total=len(gtf_df)):
            if t_row[2] not in col2_types: col2_types[t_row[2]]=0
            col2_types[t_row[2]]+=1
            if t_row[2]=='gene':
                if not len(this_gene_info[1])==0:
                    if this_gene_info[0][0] not in self.gtf_dict: self.gtf_dict[this_gene_info[0][0]]=[]
                    self.gtf_dict[this_gene_info[0][0]].append([max(0,this_gene_info[0][1]-self.updown_size),this_gene_info[0][1],('Upstream' if this_gene_info[0][3][0]=='+' else 'Downstream')+":"+this_gene_info[0][3][1:]])
                    self.gtf_dict[this_gene_info[0][0]].append([this_gene_info[0][1],this_gene_info[0][2],'gene:'+this_gene_info[0][3][1:]])
                    this_gene_info[1]=sorted(this_gene_info[1])
                    this_region=this_gene_info[1][0]
                    this_region[2]=[this_region[2]]
                    if len(this_gene_info[1])>1:
                        for next_region in this_gene_info[1][1:]:
                            if next_region[0]<=this_region[1]:
                                if next_region[1]>this_region[1]: this_region[1]=next_region[1]
                                if next_region[2] not in this_region[2]: this_region[2].append(next_region[2])
                            else:
                                this_region[2]='+'.join(this_region[2])+":"+this_gene_info[0][3][1:]
                                self.gtf_dict[this_gene_info[0][0]].append(this_region)
                                this_region=next_region
                                this_region[2]=[this_region[2]]
                    this_region[2]='+'.join(this_region[2])+":"+this_gene_info[0][3][1:]
                    self.gtf_dict[this_gene_info[0][0]].append(this_region)
                    self.gtf_dict[this_gene_info[0][0]].append([this_gene_info[0][2],this_gene_info[0][2]+self.updown_size,('downstream' if this_gene_info[0][3][0]=='+' else 'upstream')+":"+this_gene_info[0][3][1:]])
                gene_name=""
                for fieldi in t_row[8].split(';'):
                    f_lsp=fieldi.split()
                    if f_lsp[0]=='gene_name':
                        gene_name=f_lsp[1][1:-1]
                        break
                this_gene_info=[[t_row[0],int(t_row[3])-1,int(t_row[4]),t_row[6]+gene_name],[]]
            elif t_row[2]!='transcript': this_gene_info[1].append([int(t_row[3])-1,int(t_row[4]),t_row[2]])
    def _annotate_base_genes(self,output_file,cpgtype='site'):
        c_chr,c_gene_ind='',0
        mw=open(output_file,'w')
        with open(self.cpg_info.get("path"),'r') as mr:
            for line in tqdm(mr):
                line=line.strip()
                if len(line)==0: continue
                if cpgtype=='site':
                    lsp=line.split()
                    site_info=[lsp[0],int(lsp[1])-1,int(lsp[1])]+lsp[2:]
                else: site_info=line.split()
                if site_info[0]!=c_chr: c_gene_ind,c_chr=0,site_info[0]
                site_info.append([])
                if c_chr not in self.gtf_dict:
                    site_info[1],site_info[2]=str(site_info[1]),str(site_info[2])
                    site_info[-1]='/'.join(site_info[-1])
                    mw.write("{}\n".format("\t".join(site_info)))
                    continue
                if site_info[1]<self.gtf_dict[c_chr][0][0]:
                    site_info[-1].append("IG:{}^{}".format(c_chr,self.gtf_dict[c_chr][0][2].split(':')[1]))
                    site_info[1],site_info[2]=str(site_info[1]),str(site_info[2])
                    site_info[-1]='/'.join(site_info[-1])
                    mw.write("{}\n".format("\t".join(site_info)))
                    continue
                gene_ind_list,c_element_ind,before_next_gene,is_save=[],c_gene_ind,False,False
                while c_element_ind<len(self.gtf_dict[c_chr]):
                    t_type=self.gtf_dict[c_chr][c_element_ind][2].split(":")[0]
                    if site_info[1]<self.gtf_dict[c_chr][c_element_ind][0]:
                        if t_type in ['Upstream','Downstream']:
                            gene_ind_list.append(c_element_ind)
                            if site_info[1]>=self.gtf_dict[c_chr][c_element_ind-1][1]: site_info[-1].append("IG:{}@{}^{}@{}".format(self.gtf_dict[c_chr][c_element_ind-1][2].split(':')[1],self.gtf_dict[c_chr][c_element_ind-1][1]-site_info[1],self.gtf_dict[c_chr][c_element_ind][2].split(':')[1],self.gtf_dict[c_chr][c_element_ind][0]-site_info[1]))
                            site_info[1],site_info[2]=str(site_info[1]),str(site_info[2])
                            site_info[-1]='/'.join(site_info[-1])
                            mw.write("{}\n".format("\t".join(site_info)))
                            before_next_gene,is_save=True,True
                        elif t_type=='gene': pass
                        elif site_info[1]>=self.gtf_dict[c_chr][c_element_ind-2 if self.gtf_dict[c_chr][c_element_ind-1][2].split(":")[0]=='gene' else c_element_ind-1][1]: site_info[-1].append(self.gtf_dict[c_chr][c_element_ind][2].split(':')[1]+':intron')
                    elif site_info[1]>=self.gtf_dict[c_chr][c_element_ind][0] and site_info[1]<self.gtf_dict[c_chr][c_element_ind][1]:
                        if t_type in ['Upstream','Downstream']:
                            site_info[-1].append(self.gtf_dict[c_chr][c_element_ind][2].split(':')[1]+':'+t_type)
                            gene_ind_list.append(c_element_ind)
                        elif t_type=='gene': pass
                        else: site_info[-1].append(self.gtf_dict[c_chr][c_element_ind][2].split(':')[1]+':'+t_type)
                    elif site_info[1]>=self.gtf_dict[c_chr][c_element_ind][1]:
                        if t_type in ['Upstream','Downstream']: gene_ind_list.append(c_element_ind)
                    c_element_ind+=1
                    if before_next_gene: break
                if c_element_ind>=len(self.gtf_dict[c_chr]) and len(site_info)==0: site_info[-1].append("IG:{}^{}".format(self.gtf_dict[c_chr][-1][2].split(':')[1],c_chr))
                if not is_save:
                    site_info[1],site_info[2]=str(site_info[1]),str(site_info[2])
                    site_info[-1]='/'.join(site_info[-1])
                    mw.write("{}\n".format("\t".join(site_info)))
                for c_gene_ind in gene_ind_list:
                    if int(site_info[1])<self.gtf_dict[c_chr][c_gene_ind+1][1]+self.updown_size: break
        mw.close()
    def _annotate_ccre(self,input_file,output_file):
        ccre_df=pd.read_csv(self.out_dir/"encodeCcreCombined.bed",sep='\t',header=None)
        ccre_dict={name:group[[1,2,12,13]].values for name,group in ccre_df.groupby(0)} if not ccre_df.empty else {}
        c_chr,c_ind='',0
        with open(input_file,'r') as mr,open(output_file,'w') as mw:
            for line in tqdm(mr):
                site_info=line.strip().split()
                if not site_info: continue
                if len(site_info)<5: site_info.append([])
                else: site_info[-1]=site_info[-1].split('/')
                site_pos,chrom=int(site_info[1]),site_info[0]
                if chrom!=c_chr: c_chr,c_ind=chrom,0
                if c_chr in ccre_dict:
                    df_chr=ccre_dict[c_chr]
                    while c_ind<len(df_chr):
                        start,end,col12,col13=df_chr[c_ind]
                        if site_pos<start: break
                        elif site_pos<end:
                            site_info[-1].append(f"ENCODE:{col13}@{col12}")
                            break
                        c_ind+=1
                site_info[-1]='/'.join(site_info[-1])
                mw.write("\t".join(map(str,site_info))+"\n")
    def _annotate_repeats(self,input_file,output_file):
        rep_df=pd.read_csv(self.out_dir/"rmsk.txt",sep='\t',header=None,names=['bin','swScore','milliDiv','milliDel','milliIns','chr','start','end','genoLeft','strand','repName','repeat_cat','repFamily','repStart','repEnd','repLeft','id'])
        rep_dict={name:group[['start','end','repeat_cat']].values for name,group in rep_df.groupby('chr')} if not rep_df.empty else {}
        c_chr,c_ind='',0
        with open(input_file,'r') as mr,open(output_file,'w') as mw:
            for line in tqdm(mr):
                site_info=line.strip().split()
                if not site_info: continue
                pos,chrom=int(site_info[1]),site_info[0]
                site_info[-1]=site_info[-1].split('/')
                if chrom!=c_chr: c_chr,c_ind=chrom,0
                if chrom in rep_dict:
                    df_chr=rep_dict[chrom]
                    while c_ind<len(df_chr) and df_chr[c_ind,1]<=pos: c_ind+=1
                    temp_ind,rep_list=c_ind,set()
                    while temp_ind<len(df_chr) and df_chr[temp_ind,0]<=pos:
                        if df_chr[temp_ind,0]<=pos<df_chr[temp_ind,1]: rep_list.add(df_chr[temp_ind,2])
                        temp_ind+=1
                    if rep_list: site_info[-1].append(f"REPEAT:{'+'.join(sorted(list(rep_list)))}")
                site_info[-1]='/'.join(site_info[-1])
                mw.write("\t".join(map(str,site_info))+"\n")
    def _annotate_epic(self,input_file,output_file,platform_info,label,names,skip,pos_offset=0):
        df=pd.read_csv(platform_info.get("path"),sep=platform_info.get("sep",","),header=None,usecols=platform_info.get("columns"),names=names,skiprows=skip).dropna()
        df['start']=df['start'].astype(int)-(pos_offset if pos_offset>0 else 0)
        df['end']=df['start']+1 if 'end' not in df.columns else df['end'].astype(int)
        p_dict={n:g[['start','end','cpg']].sort_values('start').values for n,g in df.groupby('chr')}
        c_chr,c_ind='',0
        with open(input_file,'r') as mr,open(output_file,'w') as mw:
            for line in tqdm(mr):
                site_info=line.strip().split()
                if not site_info: continue
                pos,chrom=int(site_info[1]),site_info[0]
                site_info[-1]=site_info[-1].split('/')
                if chrom!=c_chr: c_chr,c_ind=chrom,0
                if chrom in p_dict:
                    df_c=p_dict[chrom]
                    while c_ind<len(df_c):
                        start,end,cpg=df_c[c_ind]
                        if pos<start: break
                        elif pos<end:
                            site_info[-1].append(f"{label}:{cpg}")
                            break
                        c_ind+=1
                site_info[-1]='/'.join(filter(None,site_info[-1]))
                mw.write("\t".join(map(str,site_info))+"\n")
    def generate(self,final_output_name="CpG_FullyAnnotated.bed"):
        self._prepare_output_files()
        current_input=self.cpg_info.get("path")
        tmp_files=[]
        tmp1=self.out_dir/"tmp_genes.bed"
        self._get_gene_exon_info()
        self._annotate_base_genes(tmp1)
        current_input=tmp1
        tmp_files.append(tmp1)
        if self.ccre_info and self.ccre_info.get("path") and os.path.exists(self.ccre_info.get("path")):
            tmp2=self.out_dir/"tmp_ccre.bed"
            self._annotate_ccre(current_input,tmp2)
            current_input=tmp2
            tmp_files.append(tmp2)
        if self.rmsk_info and self.rmsk_info.get("path") and os.path.exists(self.rmsk_info.get("path")):
            tmp3=self.out_dir/"tmp_repeats.bed"
            self._annotate_repeats(current_input,tmp3)
            current_input=tmp3
            tmp_files.append(tmp3)
        if self.epic_v2_info and self.epic_v2_info.get("path") and os.path.exists(self.epic_v2_info.get("path")):
            tmp4=self.out_dir/"tmp_epicv2.bed"
            self._annotate_epic(current_input,tmp4,self.epic_v2_info,'EPICv2',['cpg','chr','start'],8,pos_offset=1)
            current_input=tmp4
            tmp_files.append(tmp4)
        if self.epic_v1_info and self.epic_v1_info.get("path") and os.path.exists(self.epic_v1_info.get("path")):
            tmp5=self.out_dir/"tmp_epicv1.bed"
            self._annotate_epic(current_input,tmp5,self.epic_v1_info,'EPIC',['cpg','chr','start','end'],8,pos_offset=0)
            current_input=tmp5
            tmp_files.append(tmp5)
        if self.hm450_info and self.hm450_info.get("path") and os.path.exists(self.hm450_info.get("path")):
            tmp6=self.out_dir/"tmp_hm450.bed"
            self._annotate_epic(current_input,tmp6,self.hm450_info,'HM450',['chr','start','end','cpg'],1,pos_offset=0)
            current_input=tmp6
            tmp_files.append(tmp6)
        print(f"moving {current_input} to {self.out_dir/final_output_name}")
        shutil.move(current_input,self.out_dir/final_output_name)
        for f in tmp_files:
            if f.exists(): f.unlink()



def run_setup(genome):
    if str(genome).lower().endswith(('.yml','.yaml')):
        try: CpGAnnotator(genome).generate("CpG_gencode_annotation.bed")
        except Exception: sys.exit(1)
        return
    package_home=os.path.dirname(os.path.abspath(__file__))
    script_path=os.path.join(package_home,"bin",f"get_files_{genome}.sh")
    if not os.path.exists(script_path): sys.exit(1)
    try: subprocess.run(["bash",script_path],cwd=package_home,check=True)
    except subprocess.CalledProcessError: sys.exit(1)


