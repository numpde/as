
import math
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from progressbar import ProgressBar as Progress

IFILE = {
	'TCGA' : "OUTPUT/e_prepared/UV/tcga.pkl",
	'BCXX' : "OUTPUT/e_prepared/UV/bcxx.pkl",
}

OFILE = {
	'cos' : "./TCGA-v-BCXX.pkl",
}

# BCXX data
BCXX_DATA = pickle.load(open(IFILE['BCXX'], "rb"))
BCXX = BCXX_DATA['X']

# TCGA data
TCGA_DATA = pickle.load(open(IFILE['TCGA'], "rb"))
TCGA = TCGA_DATA['X']

# Keep only the symbols that are in both
(TCGA, BCXX) = TCGA.align(BCXX, join='inner', axis=0)

# Specialize to a subset of genes
# HALLMARK_APOPTOSIS
N = "CASP3	CASP9	DFFA	CASP7	CFLAR	BIRC3	PMAIP1	CASP8	JUN	BCL2L11	MCL1	IL1B	SPTAN1	DIABLO	BAX	BIK	IL1A	BID	CDKN1A	GADD45A	DDIT3	CDKN1B	TNF	GSN	TNFSF10	CASP6	SQSTM1	FASLG	EGR3	CD44	FAS	IL18	IGFBP6	PRF1	DAP	CCND1	BTG3	F2R	SATB1	BNIP3L	CASP4	TNFRSF12A	CREBBP	RHOB	GPX3	PDGFRB	TSPO	CCND2	XIAP	TIMP1	CTNNB1	IRF1	HSPB1	ADD1	TIMP2	BTG2	TIMP3	LEF1	CASP1	GPX1	BCL10	IGF2R	CDC25B	AIFM3	CD38	PPP3R1	HGF	CLU	ATF3	LGALS3	LUM	LMNA	GADD45B	CDK2	IFNB1	RETSAT	SMAD7	SOD1	PTK2	ENO2	HMOX1	IER3	BCL2L10	CD2	GCH1	MMP2	VDAC2	TAP1	PLAT	IFNGR1	APP	BRCA1	ROCK1	PSEN1	DCN	PSEN2	SOD2	BMF	EREG	KRT18	TGFB2	RELA	WEE1	RARA	CD14	CD69	PEA15	DNAJC3	CASP2	CTH	PLCB2	BMP2	HMGB2	LPPR4	H1F0	TGFBR3	EBP	TXNIP	ANKH	RHOT2	CYLD	GSTM1	GSR	BGN	BCL2L1	GNA15	MGMT	PPT1	F2	IL6	SC5DL	IFITM3	RNASEL	EMP1	CAV1	DNM1L	ANXA1	TOP2A	ISG20	SLC20A1	MADD	PPP2R5B	BCAP31	ERBB3	NEDD9	SAT1	PDCD4	BCL2L2	FEZ1	ERBB2	DNAJA1	DAP3	DPYD	NEFH	PAK1	FDXR	GPX4	ETF1	CCNA1	GUCY2D	AVPR1A ZZTOP"
#
N = N.split()
#
TCGA = TCGA.ix[N]
BCXX = BCXX.ix[N]

#print(TCGA)
#print(BCXX)
#exit()

C = TCGA_DATA['C']
C = C[ C['Type'] == "tumor" ]
C = C[ C['PAM50'] == "LumA" ]['Patient']
#print(C)
#exit()

M = np.zeros( (len(TCGA.columns), len(BCXX.columns)) )

for (m, a) in enumerate(TCGA) :
	
	if not a.startswith(tuple(C)) : continue
	#if not a.startswith(("TCGA-AN-A0FL", "TCGA-B6-A0RS", "TCGA-E2-A15A")) : continue
	
	for (n, b) in enumerate(BCXX) : 
		dfa = TCGA[a]
		dfb = BCXX[b]
		#df = pd.merge(dfa, dfb, left_index=True, right_index=True, how='inner')
		
		def dot(a, b) : return (a * b).sum()
		def cos(a, b) : return dot(a, b) / math.sqrt(dot(a, a) * dot(b, b))
		
		M[m, n] = cos(dfa, dfb)
	
	print(a, "x", b, " -- ", m, "out of", len(TCGA.columns))
	
	#I = sorted(range(M.shape[1]), key=(lambda n : -np.sum(M[:, n])), reverse=True)
	#M = M[:, I]
	
	J = sorted(range(M.shape[0]), key=(lambda m : +np.sum(M[m, :])), reverse=True)
	M = M[J, :]
	
	plt.ion()
	plt.clf()
	plt.imshow(M, aspect='auto')
	plt.show()
	plt.pause(0.1)

input()

#pickle.dump(
	#{
		#'M' : M,
		#'I' : TCGA.columns,
		#'J' : BCXX.columns,
	#},
	#open(OFILE['cos'], 'wb')
#)
