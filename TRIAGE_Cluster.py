import pandas as pd 
import numpy as np 
from scipy import stats
from sklearn.cluster import DBSCAN

np.random.seed(0)

# PART 1: Assign TRIAGE gene and its RTS score to each gene, which will be used later to calcuate repressive tendency density, used in the contour clustering. 
def getv1gene(exprdf, cellmeta, fdir = './', cellcolname = 'cell_barcode'):
	# load in TRIAGE version 1 rts
	rtslist = pd.read_table(fdir + 'TRIAGE1/repressive_hg19_cont_prop.txt', header = None, names = ['geneName', 'RTS'])
	# order rtslist by rts and only keep genes in expression matrix
	v1_new = rtslist[rtslist['geneName'].isin(exprdf.index)].sort_values('RTS', ascending = False)
	v1_new.index = range(0, len(v1_new)) # reset indices in v1_new
	# order expression matrix genes by rts
	exprdford = exprdf = exprdf.loc[v1_new['geneName'].values]
	# asisgn highest RTS score gene to each cell by getting id of first non-zero row in each col
	genelabels = exprdford.ne(0).idxmax()
	# put into dataframe and merge with existing metadata
	genelabelsdf = pd.DataFrame(list(zip(genelabels.index, genelabels.values)), columns = [cellcolname, 'v1gene'])
	outmeta = genelabelsdf.merge(v1_new.rename(columns = {'geneName':'v1gene'}), how = 'left', on = 'v1gene')
	outmeta = outmeta.merge(cellmeta, how = 'outer', on = cellcolname)
	return(outmeta)

################################
# PART 2: Take output from getv1gene and a given density calculation bandwidth, calculate repressive tendency density across umap space and find 'peak' regions of high rt density, i.e., highly defined cell types

# function called in getPeaks
# adds clusters that contain other clusters into the blacklist
# i.e. go through all the peaks in the given list, and blacklist each string without the last cluster id. so if the peak is 0_1_2_3_4, blacklist 0, 0_1, 0_1_2, 0_1_2_3 because they all contain the 0_1_2_3_4 peak (cluster 4 on the 5th level)
def getBlacklist(blacklistthese):
	out = []
	if isinstance(blacklistthese, list):
		for peak in blacklistthese:
			levels = peak.split('_')
			out = out + ['_'.join(levels[0:n]) for n in range(1, len(levels))]
	else: # if there's only one peak in blacklistthese
		levels = blacklistthese.split('_')
		out = out + ['_'.join(levels[0:n]) for n in range(1,len(levels))]
	return list(set(out))


# goes through all cluster strings (for each cell, the cluster it is in on every contour clust level) and gets only those that represent peaks (are a region of cells that don't contain smaller regions within them)
def getPeaks(clusterstrings):
	usc = list(set(clusterstrings)) # get unique clusterstrings
	peaks = [] # what we want to keep
	blacklist = [] # the clusters that have smaller clusters inside them
	# loop through clusterstrings and get peaks
	while len(usc) > 0 and '_' in usc[0]: 
		newpeaks = [u for u in usc if u.endswith('nan') == False] # we only want to keep those that don't end with none/NA
		print(newpeaks)
		peaks = list(set(peaks + newpeaks))
		# blacklist all the clusters that the peaks on this level are in 
		newblacklist = getBlacklist(peaks)
		print('Newly blacklisted', len(newblacklist), 'clusters')
		# print(set(newblacklist))
		blacklist = list(set(blacklist + newblacklist))
		# remove any peaks that are in the updated blacklist
		usc = [peak for peak in ['_'.join(u.split('_')[:-1]) for u in usc] if peak not in blacklist]
		if len(usc) == 1 and usc[0] == 'nan':
			usc = []
		print('Remaning clusters: ', len(usc))
	peaks = peaks + [u for u in usc if u != 'nan']
	# add nan_ * (total removed levels) to end of each peak so that their total length matches the original cluster string labels
	peaknames = ['_'.join(peak.split('_') + ['nan'] * (len(clusterstrings[0].split('_')) - len(peak.split('_')))) for peak in peaks]
	return(peaknames)


# the function that runs the contour clustering itself, have to choose bandwith prior to running it, usually do this by making a few kde pots of different bandwith (bw) values 
# df is pandas DataFrame with cell barcodes, UMAP coordinate columns, any other metadata, and the RTS values for the TRIAGE gene chosen for each cell
def runcclust(bandwidth, df):
	print('calculating density landscape...')
	kde = stats.gaussian_kde(df.loc[:, ["UMAP1", "UMAP2"]].T, weights=df["RTS"])
	kde.set_bandwidth(kde.factor * bandwidth)
	df["density"] = kde(df.loc[:, ["UMAP1", "UMAP2"]].T)
	# cut density into 10 contours
	dens_sorted = np.sort(df['density'])
	p = 1. * np.arange(len(df[['density']])) / (len(df[['density']]) - 1)
	thresholds = [dens_sorted[np.where(p == min(p, key=lambda list_value:abs(list_value - x)))] for x in np.arange(0, 1, 0.1)]
	print('finding clusters on each contour level...')
	df_new = df
	for threshold in thresholds:
		ctypes = df[df['density']>threshold[0]]
		dbscans = DBSCAN(eps = 0.2).fit(ctypes.loc[:,['UMAP1', 'UMAP2']])
		# print(np.asarray(np.unique(dbscans.labels_, return_counts=True)).T)
		#sns.relplot(data=ctypes, x='UMAP1', y='UMAP2', hue=dbscans.labels_.astype(str), s=1, edgecolor=None, legend='full', facet_kws={'legend_out':True}); plt.suptitle(threshold[0]); plt.show()
		ctypes['cont_'+str(threshold[0])] = dbscans.labels_
        # "cell" is the cell names column in metadata
		df_new = pd.merge(df_new, ctypes[['cell', 'cont_'+str(threshold[0])]], how='left', on='cell')
	# getPeaks function is from misc/contourclust_getpeakfns.py 
	print('getting peaks..')
	cccolidx = [col for col in df_new.columns if col.startswith('cont_')] 
	cccols = df_new[cccolidx[1:]] 
	strcc = cccols.apply(lambda row: '_'.join(row.values.astype('str')), axis = 1)
	peaks = getPeaks(strcc)
	df_new['peaks'] = [s if s in peaks else np.nan for s in strcc]
	df_new['clust'] = [clust if clust != -1 else np.nan for clust in pd.factorize(df_new['peaks'])[0]]
	firstcellsidx = [list(df_new['clust']).index(clust) for clust in df_new['clust'].unique()[1:]]
	firstcells = [df_new['clust'][i] if i in firstcellsidx else np.nan for i in range(0, len(df_new))]
	df_new['firstcells'] = firstcells
	# return metadata table
	return(df_new)



from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns 
import numpy as np
import matplotlib
from scipy import stats
from sklearn.cluster import DBSCAN
from scipy.io import mmread
import sys
from pathlib import Path


path = Path("/home/uqysun19/90days/manuscript_TRIAGEclustering/bottom2top_RTS/Gottgens_data")


df = pd.read_csv(path / "Gottgens_originalExpr_TRIAGEclustering/Gottgens_originalExpr_metadata_with_v1gene_RTS_density.csv", index_col=0)

# make initial contour plot based on RTS values, choose bandwidth and save #############
c = np.linspace(0,1,10)
colors = plt.get_cmap('Spectral', 10)(c)
cmap = matplotlib.colors.ListedColormap(colors[::-1])
a4_dims = (11.7, 8.27)

bw = (float(sys.argv[1]) / 10)
fig, ax = plt.subplots(figsize=a4_dims)
sns.kdeplot(data=df, x="UMAP1", y="UMAP2", weights="RTS", bw_adjust=bw, levels=10, cmap=cmap, fill = True, cbar=True, ax= ax)
sns.kdeplot(data=df, x="UMAP1", y="UMAP2", weights="RTS", bw_adjust=bw, levels=10, color='black', ax=ax, **{"linestyles":"solid", "linewidths":0.1})
ax.set_title('bw:'+ str(bw))
# fig.show()

fig.savefig(path / 'Gottgens_originalExpr_TRIAGEclustering/multi_bw/Gottgens_originalExpr_contour_bw{}.pdf'.format(str(bw)), dpi = 300)

# run contour clustering ##############################
# takes bandwidth we chose and runs clustering on each of the 10 levels and retrieves the regions where local RTS score is at a peak
df_new = runcclust(bw, df)

# show plot of cclustering output and save ##############################
colours = ['#e5e5e5'] + sns.color_palette('husl', (len(set(df_new['peaks'])))-1).as_hex()

# umap of peaks without label
a = sns.relplot(data = df_new, x = 'UMAP1', y = 'UMAP2', hue = df_new['clust'].astype('str'), s = 1, 
                edgecolor = None, palette = colours,legend = 'full', 
                facet_kws = {'legend_out':True}, rasterized=True)
a.savefig(path / 'Gottgens_originalExpr_TRIAGEclustering/multi_bw/Gottgens_originalExpr_bw{}_nolabelUMAP.pdf'.format(str(bw)), dpi = 300)

# umap of peaks with label
sns.relplot(data = df_new, x = 'UMAP1', y = 'UMAP2', hue = df_new['clust'].astype('str'), 
	s = 1, edgecolor = None, palette = colours,legend = 'full', facet_kws = {'legend_out':True}, rasterized=True)
for i in range(0, df_new.shape[0]):
	if df_new['firstcells'][i].astype('str') != 'nan':
		plt.text(x = df_new.UMAP1[i]+0.3, y = df_new.UMAP2[i]+0.3, s = df_new.firstcells[i].astype('int'), size = 5)
# plt.show()

plt.rcParams['pdf.fonttype'] = 42
plt.savefig(path / 'Gottgens_originalExpr_TRIAGEclustering/multi_bw/Gottgens_originalExpr_bw{}_labelledUMAP.pdf'.format(str(bw)), dpi = 300)

# save ouptut
df_new.to_csv(path / 'Gottgens_originalExpr_TRIAGEclustering/multi_bw/Gottgens_originalExpr_bw{}_metadata.csv'.format(str(bw)))

# #get metadata and matrix with only peak cells
# df_sub = df_new[df_new['clust'].notna()]
# df_sub.to_csv(path / 'Gottgens_originalExpr_TRIAGEclustering/multi_bw/Gottgens_originalExpr_bw{}_cleanmeta.csv'.format(str(bw)))

# mtx = pd.read_csv(path / "Gottgens_original_logged_matrix.csv", index_col=0)
# mtx_sub = mtx[mtx.columns.intersection(df_sub['cell'])]
# mtx_sub.to_csv(path / 'Gottgens_originalExpr_TRIAGEclustering/multi_bw/Gottgens_originalExpr_bw{}_submtx.csv'.format(str(bw)))




































