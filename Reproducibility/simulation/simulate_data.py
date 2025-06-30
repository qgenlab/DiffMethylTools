import pandas as pd
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
import queue
import random
from scipy.stats import truncnorm
from statsmodels.nonparametric.smoothers_lowess import lowess 
import sys




def cluster_cpgs(bed_file, max_gap=100):
    df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'strand'])
    df = df[df["strand"] == "+"]
    df = df.sort_values(['chr', 'start'])
    regions = []
    current_chr = None
    region_start = region_end = None
    for idx, row in df.iterrows():
        chr_, start, end = row['chr'], row['start'], row['end']
        if current_chr != chr_ or start - region_end > max_gap:
            if current_chr is not None:
                regions.append((current_chr, region_start, region_end))
            current_chr = chr_
            region_start = start
            region_end = end
        else:
            region_end = max(region_end, end)
    regions.append((current_chr, region_start, region_end))
    out_df = pd.DataFrame(regions, columns=['chr', 'start', 'end'])
    return out_df



path = "."
CpG_file = sys.argv[1]
Methylation_file = sys.argv[2]

output_df = cluster_cpgs(CpG_file, max_gap=100)
output_df[(output_df["end"] - output_df["start"]) > 20].to_csv("./CpG_regions_long_hg38.bed", sep='\t', header=False, index=False)

methyl_data = pd.read_csv(Methylation_file, sep="\t",  names=["chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes"])

output_df = pd.read_csv("./CpG_regions_long_hg38.bed", sep='\t', names=["chr", "start", "end"])


valid_region_data = queue.Queue()
num = 0


chrs = output_df['chr'].unique()

for chr in chrs:
	print(chr)
	output_df_1 = output_df[output_df["chr"] == chr]
	methyl_data_1 = methyl_data[methyl_data["chrom"] == chr]
	for _, region in output_df_1.iterrows():
		chrom, r_start, r_end = region["chr"], region["start"], region["end"]   
		overlapping = methyl_data_1[(methyl_data_1["chromStart"] >= r_start) & (methyl_data_1["chromEnd"] <= r_end)]    
		if len(overlapping) > 20:
			valid_region_data.put(overlapping)
			num += 1
			if num % 100 == 0: print(num)







def generate_case_from_control(control_df, min_diff=5, direction='increase', smoothing_sigma=5):
    df = control_df.copy()
    meth = df['y'].values
    n = len(meth)
    offset = np.random.uniform(min_diff, min_diff, size=n)
    if direction == 'both':
        sign = np.random.choice([-1, 1], size=n)
    elif direction == 'decrease':
        sign = -1
    else:  
        sign = 1
    signed_offset = offset * sign
    smooth_offset = gaussian_filter1d(signed_offset, sigma=smoothing_sigma)
    case_meth = meth + smooth_offset
    case_meth = np.clip(case_meth, 0, 100)
    df['y'] = case_meth
    return df






def generate_samples(smoothed_df, n, std_case):
    positions = smoothed_df['x'].values
    smoothed_values = smoothed_df['y'].values
    samples = truncnorm.rvs(0, 100, loc=smoothed_values[:, None], scale=std_case, size=(len(smoothed_values), n))
    sample_df = pd.DataFrame(samples, columns=[f"meth_{i+1}" for i in range(n)])
    sample_df = sample_df.clip(lower=0, upper=100)
    sample_df.insert(0, "x", positions) 	
    return sample_df


smoothie = []

valid_region_data = list(valid_region_data.queue) 
for i, region in enumerate(valid_region_data):
	y = region["blockSizes"]
	x = region["chromStart"]
	cov = region["blockCount"]
	chr = region["chrom"]
	smoothed = lowess(y, x, frac=0.3)
	x_smooth, y_smooth = smoothed[:, 0], smoothed[:, 1]
	d = {"chr":chr, "x":x_smooth, "coverage":cov, "y":y_smooth}
	smoothie.append(pd.DataFrame(d))
	if i%1000 == 0: print(i, pd.DataFrame(d))


smoothie_case = []

diff = [0, 5, 10, 15, 20, 30, 40, 50]

weights = [0.4] *1 + [0.3 / 1] * 1 + [0.2 / 2] * 2  + [0.05 / 4] * 4

n_samples = 3

direction = ['decrease', 'increase']
samples_case = []
samples_ctr = []
for i, region in enumerate(smoothie):
	random_min_diff = random.choices(diff, weights=weights, k=1)[0]
	random_direction = random.choice(direction)
	if region["y"].mean() > 75: random_direction = 'decrease'; 
	elif region["y"].mean() < 25: random_direction = 'increase'; 
	data = generate_case_from_control(region, min_diff=random_min_diff, direction=random_direction)
	data["direction"] = [random_direction] * data.shape[0]
	data["min_diff"] = [random_min_diff] * data.shape[0]
	smoothie_case.append(data)
	stdv_case = random.randint(2, 15)
	stdv_ctr = random.randint(2, 15)
	data["stdv_case"] = [stdv_case] * region.shape[0]
	samples_case.append(generate_samples(data, n_samples, stdv_case))
	region["stdv_ctr"] = [stdv_ctr] * region.shape[0]
	samples_ctr.append(generate_samples(region, n_samples, stdv_ctr))
	if i%1000 == 0: print(i, samples_ctr[-1], samples_case[-1])



import numpy as np

def simulate_coverage_from_methylation(real_coverage, real_meth, sim_meth, base_stdv=0.1, noise=2):
    real_coverage = np.asarray(real_coverage)
    real_meth = np.asarray(real_meth)
    sim_meth = np.asarray(sim_meth)
    deviation = np.abs(sim_meth - real_meth) / 100  
    scaling = np.exp(-deviation / base_stdv) 
    simulated_cov = real_coverage * scaling
    simulated_cov += np.random.normal(loc=0, scale=noise, size=simulated_cov.shape)
    simulated_cov = np.clip(simulated_cov, 1, None).astype(int)
    return simulated_cov




new_samples_case = []
new_samples_ctr = []
for i, cpg in enumerate(smoothie):
	coverage = cpg["coverage"]
	base_std = cpg["stdv_ctr"].iloc[0]
	real_meth = cpg["y"]
	sim_meth_ = samples_ctr[i].filter(like="meth_")
	new_df_ctr = samples_ctr[i]
	for j in range(1, n_samples + 1):
		new_df_ctr[f'cov_{j}'] = simulate_coverage_from_methylation(coverage, real_meth, sim_meth_[f'meth_{j}'], base_std)
	new_samples_ctr.append(new_df_ctr)
	coverage = smoothie_case[i]["coverage"]
	base_std = smoothie_case[i]["stdv_case"].iloc[0]
	real_meth = smoothie_case[i]["y"]
	sim_meth_ = samples_case[i].filter(like="meth_")
	new_df_case = samples_case[i]
	for j in range(1, n_samples + 1):
		new_df_case[f'cov_{j}'] = simulate_coverage_from_methylation(coverage, real_meth, sim_meth_[f'meth_{j}'], base_std)
	new_samples_case.append(new_df_case)
	if i%1000 == 0: print(i, samples_case[i], new_df_case, samples_ctr[i],new_df_ctr)





case = pd.concat(new_samples_case, axis = 0)
ctr = pd.concat(new_samples_ctr, axis = 0)

case.reset_index(inplace = True)
ctr.reset_index(inplace = True)



og_case = pd.concat(smoothie_case, axis = 0)
og_ctr = pd.concat(smoothie, axis = 0)


og_case.reset_index(inplace = True)
og_ctr.reset_index(inplace = True)


case.drop('index', axis=1, inplace=True)
ctr.drop('index', axis=1, inplace=True)


og_case.drop('index', axis=1, inplace=True)
og_ctr.drop('index', axis=1, inplace=True)

all_case = pd.concat([case.drop('x', axis=1), og_case], axis = 1)

all_ctr = pd.concat([ctr.drop('x', axis=1), og_ctr], axis = 1)




mean1 = (all_ctr.filter(like="meth_") / 100).mean(axis=1)
mean2 = (all_case.filter(like="meth_") / 100).mean(axis=1)

std1 = (all_ctr.filter(like="meth_") / 100).std(axis=1)
std2 = (all_case.filter(like="meth_") / 100).std(axis=1)

pooled_std = np.sqrt(((std1 ** 2).apply(lambda x: max(abs(x), 0.075)) + (std2 ** 2).apply(lambda x: max(abs(x), 0.075)) / 2))

cohen_d = (mean1 - mean2) / pooled_std
cohen_d = cohen_d.replace([np.inf, -np.inf], 0).fillna(0)




for k in range(1, n_samples+1):
	df_x = all_case[["chr", "x", f"cov_{k}",f"meth_{k}"]].rename(columns={'x': 'start', f"cov_{k}": 'coverage', f"meth_{k}": 'methylation'})
	df_x = df_x.loc[:, ~df_x.columns.duplicated()]
	df_x.to_csv(f"{path}/simulation_case_{k}.csv", index=False, sep="\t")
	df_x = all_ctr[["chr", "x", f"cov_{k}",f"meth_{k}"]].rename(columns={'x': 'start', f"cov_{k}": 'coverage', f"meth_{k}": 'methylation'})
	df_x = df_x.loc[:, ~df_x.columns.duplicated()]
	df_x.to_csv(f"{path}/simulation_ctr_{k}.csv", index=False, sep="\t")

cases = []
ctrs = []
n_samples = 3

for k in range(1, n_samples+1):
	cases.append(pd.read_csv(f"{path}/simulation_case_{k}.csv",  sep="\t"))
	ctrs.append(pd.read_csv(f"{path}/simulation_ctr_{k}.csv",  sep="\t"))



x1 = all_ctr[["chr", "x", "stdv_ctr", "y"]].rename(columns={'y': 'y_ctr'})
x1 = x1.loc[:, ~x1.columns.duplicated()]
x2 = all_case[["direction", "min_diff", "stdv_case", "y"]].rename(columns={'y': 'y_case'})

pd.concat([x1,x2, cohen_d], axis=1).to_csv("{path}/res_sim.csv", sep="\t", index=False)




std = pd.concat([(all_ctr.filter(like="meth_") / 100), (all_case.filter(like="meth_") / 100)], axis = 1).std(axis=1)

pos = x1[["chr", "x"]]


pd.concat([pos,std, cohen_d], axis=1).rename(columns={0: 'std_all', 1:'cohen_d'}).to_csv("{path}/sim.pos.cohen.csv", sep="\t", index=False)


