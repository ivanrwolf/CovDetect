#! /usr/local/bin/miniconda2/bin/python
# Coverage detect version 1.0
# This software was written by Ivan R. Wolf
# Python modules required: Numpy, Pandas
#
import numpy as np
import pandas as pd
import collections
import argparse

def create_blocks(df, scfld, basepairs):
# Function to create blocks of genomic regions detected
	blocks = []
	x = 0
	s = []
	#basepairs = 100
	if isinstance(df.loc[scfld]['base'], np.float64):
		return False
	else:
		for i in df.loc[scfld]['base'].tolist():
			if (i - x + 1) <= basepairs: # Determine range of blocks
				s.append(i)
				x = i
			elif len(s) == 0:
					s.append(i)
					x = i
			else:
				if min(s) != max(s):
					blocks.append((scfld, min(s),max(s)))
					s = []
					s.append(i)
					x = i
				else:
					s = []
					s.append(i)
					x = i
	return blocks

def save_log(gsize, sigmaCi_0B, CovMed_0B, cutoff_0B, sigmaCi_1B, CovMed_1B, cutoff_1B, selected_scfld, AvgCovNorm1B, AvgCovNorm0B, meanRatio_0B, stdv_0B, output_log):
	# Save log

	f=open(output_log,"w")
	print >>f, "Genome size:", gsize
	print >>f, "0B Coverage sum:", sigmaCi_0B
	print >>f, "0B Avg Coverage:", CovMed_0B
	print >>f, "0B Cutoff Value:", cutoff_0B

	print >>f, "1B Coverage sum:", sigmaCi_1B
	print >>f, "1B Avg Coverage:", CovMed_1B
	print >>f, "1B Cutoff Value:", cutoff_1B

	print >>f, "Scafold selected for normalization:", selected_scfld

	print >>f, "Average of normalized coverage for 1B:", AvgCovNorm1B
	print >>f, "Average of normalized coverage for 0B:", AvgCovNorm0B

	print >>f, "Mean Ratio:", meanRatio_0B
	print >>f, "Mean Ratio STDV:", stdv_0B
	f.close()

	return

def main(args):
	output_table = "Seleted_bases_" + str(args.input) + "_STDV" + str(args.standev) + "_BP" + str(args.basepairs) + ".txt"
	output_blocks = "Seleted_bases_" + str(args.input) + "_STDV" + str(args.standev) + "_BP" + str(args.basepairs) + ".blocks.txt"
	output_log = "Seleted_bases_" + str(args.input) + "_STDV" + str(args.standev) + "_BP" + str(args.basepairs) + ".log.txt"

	df = pd.read_csv(args.input,index_col=0,header=None,engine='c',sep='\t') # Load Data
	df.columns = ["base","cov_0B","cov_1B"] # Rename Columns

	gsize = df.shape[0]

	sigmaCi_0B = df['cov_0B'].sum()
	CovMed_0B = sigmaCi_0B/gsize
	cutoff_0B = CovMed_0B / 2

	sigmaCi_1B = df['cov_1B'].sum()
	CovMed_1B = sigmaCi_1B/gsize
	cutoff_1B = CovMed_1B / 2

	df['size'] = pd.DataFrame.from_dict(collections.Counter(df.index.values),orient='index') # Add a size colum
	df['norm_0B'] = (df['cov_0B'] / df['size']) / CovMed_0B # Add a column with normalized value for 0B
	df['norm_1B'] = (df['cov_1B'] / df['size']) / CovMed_1B # Add a column with normalized value for 1B

	# Get scafold info
	scfld_info_0B = pd.concat([df.groupby(level=0).var()['cov_0B'], 
	pd.DataFrame.from_dict(collections.Counter(df.index.values),orient='index')], axis=1)
	scfld_info_0B.columns = ["cov_var","size"]
	scfld_0B_list = scfld_info_0B.sort_values(['cov_var', 'size'], ascending=[True, False]).index.values.tolist()


	scfld_info_1B = pd.concat([df.groupby(level=0).var()['cov_1B'], 
	pd.DataFrame.from_dict(collections.Counter(df.index.values),orient='index')], axis=1) # Join information
	scfld_info_1B.columns = ["cov_var","size"]
	scfld_1B_list = scfld_info_1B.sort_values(['cov_var', 'size'], ascending=[True, False]).index.values.tolist()

	# Select scafold and calculate mean ratio
	selected_scfld = [i for i, j in zip(scfld_0B_list, scfld_1B_list) if i == j][0]

	AvgCovNorm1B = df.loc[selected_scfld]['norm_1B'].mean()
	AvgCovNorm0B = df.loc[selected_scfld]['norm_0B'].mean()

	meanRatio_0B = AvgCovNorm1B / AvgCovNorm0B

	stdv_0B = df.loc[selected_scfld]['norm_0B'].std()

	# Filter table based on cutoffs
	df = df[df['cov_1B'] >= cutoff_1B]
	df = df[df['cov_0B'] >= cutoff_0B]

	# Add the ratio for 1B
	df['ratio_1B'] = df['norm_1B'] / df['norm_0B']

	# Save bases with additional info
	#(df[df['ratio_1B'] >= (meanRatio_0B + stdv_0B)]).to_csv("Selected_bases.txt",sep='\t')

	df = df[df['ratio_1B'] >= (meanRatio_0B + (stdv_0B * args.standev))] # Determine number of stdv

	save_log(gsize, sigmaCi_0B, CovMed_0B, cutoff_0B, sigmaCi_1B, CovMed_1B, cutoff_1B, selected_scfld, AvgCovNorm1B, AvgCovNorm0B, meanRatio_0B, stdv_0B, output_log)
	df.to_csv(output_table,sep='\t')

	block_df = pd.DataFrame()
	for scfld in set(df.index.values):
		block_list = create_blocks(df,scfld, args.basepairs)
		if block_list == False:
			continue
		else:
			block_df = pd.concat([block_df, pd.DataFrame(block_list)])

	block_df.columns = ["Scafold","start","end"]
	block_df.to_csv(output_blocks,sep='\t',index=False)

def init():
	parser = argparse.ArgumentParser(description="Retrieve regions with coverage above threshould between two different experiments")
	parser.add_argument("input", help="Input File")
	parser.add_argument("-bp", "--basepairs", type=int, default=100, help="Determine base pairs for sliding bases")
	parser.add_argument("-stdv", "--standev", type=int, default=2, help="Number of standard deviations used to select the bases")
	args = parser.parse_args()
	return args

if __name__ == "__main__":
	args = init()
	main(args)
