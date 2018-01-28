
# RA, 2018-01-28


import os
import math
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from multiprocessing import cpu_count
from joblib import Parallel, delayed
from progressbar import ProgressBar as Progress

IFILE = {
	'TCGA' : "OUTPUT/e_prepared/UV/tcga.pkl",
	'BCXX' : "OUTPUT/e_prepared/UV/bcxx.pkl",
}

OFILE = {
	'cos-img' : "OUTPUT/g_reclass_pam50/TCGA-v-BCXX.pdf",
	
	'cos-data' : "OUTPUT/g_reclass_pam50/TCGA-v-BCXX.pkl",
}

# Create output directories
for f in OFILE.values() :
	os.makedirs(os.path.dirname(f), exist_ok=True)

PARAM = {
	# Number of parallel computing processes
	'#proc' : min(12, math.ceil(cpu_count() / 1.5)),
}



def axb(a, b) :
	def dot(a, b) : return (a * b).sum()
	return dot(a, b) / math.sqrt(dot(a, a) * dot(b, b))


# BCXX data
BCXX_DATA = pickle.load(open(IFILE['BCXX'], "rb"))
BCXX = BCXX_DATA['X']

# TCGA data
TCGA_DATA = pickle.load(open(IFILE['TCGA'], "rb"))
TCGA = TCGA_DATA['X']

# Keep only the symbols that are in both
(TCGA, BCXX) = TCGA.align(BCXX, join='inner', axis=0)

# Specialize to a subset of genes
#
N = None
#
# HALLMARK_APOPTOSIS
N = "CASP3	CASP9	DFFA	CASP7	CFLAR	BIRC3	PMAIP1	CASP8	JUN	BCL2L11	MCL1	IL1B	SPTAN1	DIABLO	BAX	BIK	IL1A	BID	CDKN1A	GADD45A	DDIT3	CDKN1B	TNF	GSN	TNFSF10	CASP6	SQSTM1	FASLG	EGR3	CD44	FAS	IL18	IGFBP6	PRF1	DAP	CCND1	BTG3	F2R	SATB1	BNIP3L	CASP4	TNFRSF12A	CREBBP	RHOB	GPX3	PDGFRB	TSPO	CCND2	XIAP	TIMP1	CTNNB1	IRF1	HSPB1	ADD1	TIMP2	BTG2	TIMP3	LEF1	CASP1	GPX1	BCL10	IGF2R	CDC25B	AIFM3	CD38	PPP3R1	HGF	CLU	ATF3	LGALS3	LUM	LMNA	GADD45B	CDK2	IFNB1	RETSAT	SMAD7	SOD1	PTK2	ENO2	HMOX1	IER3	BCL2L10	CD2	GCH1	MMP2	VDAC2	TAP1	PLAT	IFNGR1	APP	BRCA1	ROCK1	PSEN1	DCN	PSEN2	SOD2	BMF	EREG	KRT18	TGFB2	RELA	WEE1	RARA	CD14	CD69	PEA15	DNAJC3	CASP2	CTH	PLCB2	BMP2	HMGB2	LPPR4	H1F0	TGFBR3	EBP	TXNIP	ANKH	RHOT2	CYLD	GSTM1	GSR	BGN	BCL2L1	GNA15	MGMT	PPT1	F2	IL6	SC5DL	IFITM3	RNASEL	EMP1	CAV1	DNM1L	ANXA1	TOP2A	ISG20	SLC20A1	MADD	PPP2R5B	BCAP31	ERBB3	NEDD9	SAT1	PDCD4	BCL2L2	FEZ1	ERBB2	DNAJA1	DAP3	DPYD	NEFH	PAK1	FDXR	GPX4	ETF1	CCNA1	GUCY2D	AVPR1A ZZTOP"
#
# HALLMARK_DNA_REPAIR
#N = "POLR2H	POLR2A	POLR2G	POLR2E	POLR2J	POLR2F	POLR2C	POLR2K	GTF2H3	POLR2D	ERCC3	DDB2	POLR1C	XPC	PCNA	POLR2I	SUPT4H1	POLD3	POLR3GL	POLR3C	GTF2B	POLR1D	NCBP2	RDBP	GTF2F1	ERCC5	LIG1	ERCC1	ERCC4	POLD4	COBRA1	RFC2	ELL	TAF10	RRM2B	SUPT5H	RPA3	SNAPC5	SSRP1	RFC3	RPA2	TCEB3	TAF12	TH1L	TAF13	TAF6	TAF9	GTF2A2	VPS37D	NME1	RNMT	ERCC2	POLE4	VPS37B	NT5C3	SNAPC4	AAAS	ZNRD1	RFC4	ITPA	POM121	BRF2	RFC5	SAC3D1	CLP1	NME4	PRIM1	VPS28	TSG101	USP11	TAF1C	TARBP2	POLH	CETN2	POLD1	CANT1	PDE4B	DGCR8	RAD51	SURF1	PNP	ADA	NME3	GTF3C5	NT5C	AK1	GTF2H1	HCLS1	APRT	ERCC8	IMPDH2	POLB	SDCBP	SF3A3	DAD1	UPF3B	GUK1	TP53	ADRM1	SEC61A1	POLA2	FEN1	ZNF707	NUDT9	PDE6G	TYMS	BCAP31	DDB1	NFX1	RAD52	ADCY6	ARL6IP1	DGUOK	POLL	SMAD5	MPG	DUT	POLA1	EIF2C4	RALA	ZWINT	BCAM	TK2	CSTF3	GTF2H5	HPRT1	BOLA2	GPX4	BRP44	CDA	THOC4	MRPL40	NPR2	REV3L	EDF1	DFNA5	TMED2	STX3	RAE1	UMPS	EIF1B	AK3	NUDT21	RBX1	SRSF6	GMPR2	DCTN4	COX17	CMPK2	CCNO"
#
# HALLMARK_ESTROGEN_RESPONSE_EARLY
#N = "GREB1	CA12	SLC9A3R1	MYB	ANXA9	IGFBP4	SYBU	NPY1R	PDZK1	NRIP1	MLPH	HSPB8	EGR3	KRT19	LRIG1	KDM4B	PGR	RHOBTB3	TPD52L1	ELOVL2	RET	TPBG	TFF1	MAPT	SCNN1A	ABAT	FLNB	XBP1	CELSR2	RAB31	MYBL1	MREG	FAM102A	MSMB	STC2	FAM134B	SIAH2	ZNF185	SLC19A2	SLC1A4	FHL2	BCL2	PMAIP1	AREG	OVOL2	TSKU	ADCY9	RASGRP1	MUC1	KAZN	SLC27A2	FKBP4	CXCL12	TMPRSS3	RARA	IL17RB	CBFA2T3	TFF3	UGCG	CCND1	SLC22A5	WFS1	PTGES	WWC1	WISP2	MYC	ITPK1	TMEM164	ARL3	MED13L	SEMA3B	KRT18	SLC16A1	TJP3	SLC26A2	FAIM3	SULT2B1	SNX24	TFAP2C	TTC39A	GJA1	PRSS23	OLFM1	RAPGEFL1	ASB13	TIPARP	ABCA3	FRK	DHRS2	AQP3	KCNK15	TGIF2	FOXC1	ELF3	REEP1	PEX11A	PODXL	KLF4	BAG1	CELSR1	PLA2G16	SLC7A5	MPPED2	TIAM1	CLDN7	MYOF	RBBP8	OLFML3	GFRA1	FARP1	SVIL	TGM2	DEPTOR	CYP26B1	PAPSS2	SLC1A1	DLC1	JAK2	AFF1	KLK10	P2RY2	BLVRB	CISH	GLA	ADD3	PDLIM3	FAM63A	FOS	KRT8	SLC37A1	B4GALT1	CALCR	ESRP2	IGF1R	NBL1	SFN	OPN3	ABHD2	AR	SLC39A6	SYT12	CD44	MED24	BCL11B	CANT1	KRT13	KRT15	TOB1	SLC7A2	LAD1	TUBB2B	TBC1D30	SEC14L2	ENDOD1	HR	SCARB1	NCOR2	RHOD	INPP5F	PPIF	DHRS3	FDFT1	GAB2	UNC119	KLF10	HES1	FKBP5	SLC2A1	AMFR	NADSYN1	INHBB	BHLHE40	CALB2	FASN	CHPT1	MYBBP1A	ELOVL5	DYNLT3	ABLIM1	SOX3	SLC24A3	RAB17	MAST4	KCNK5	ELF1	RPS6KA2	ISG20L2	IL6ST	SYNGR1	SH3BP5	ALDH3B1	THSD4	CLIC3	NXT1	NAV2	RRP12	ADCY1	DHCR7	MICB	AKAP1"
#
# HALLMARK_ESTROGEN_RESPONSE_LATE
#N = "TFF1	SLC9A3R1	TPD52L1	PRSS23	CA12	PDZK1	ANXA9	CELSR2	TJP3	PGR	RET	MYB	TPBG	EGR3	ARL3	OLFM1	NPY1R	SCNN1A	XBP1	AREG	IL17RB	NRIP1	ASS1	TFF3	FKBP4	SLC27A2	SEMA3B	GPER	LLGL2	AGR2	KRT19	WISP2	BLVRB	FLNB	PDCD4	CALCR	IGFBP4	DNAJC12	TIAM1	TSPAN13	CXCL12	RAB31	PKP3	CYP26B1	FKBP5	SIAH2	ISG20	TMPRSS3	SERPINA3	WFS1	MAPT	PDLIM3	RBBP8	GJB3	PRLR	SLC1A4	FOS	PLA2G16	SLC7A5	SERPINA5	IMPA2	DHCR7	MYOF	CDH1	EMP2	OVOL2	DLG5	SOX3	CHPT1	KLK10	ELOVL5	RAPGEFL1	JAK2	SLC26A2	SLC22A5	ITPK1	PCP4	PAPSS2	NAB2	FAM102A	BCL2	LSR	CACNA2D2	CA2	ASCL1	ACOX2	CISH	GLA	PTGES	PERP	OPN3	KRT13	HSPB8	UGDH	CLIC3	KLK11	PLAC1	ABHD2	SCARB1	DCXR	CCND1	SFN	ABCA3	SULT2B1	CCNA1	STIL	MICB	ZFP36	CAV1	NBL1	CD44	HR	HOMER2	BTG3	GAL	ETFB	BAG1	FRK	SLC16A1	AFF1	TFAP2C	IGSF1	HPRT1	CDC6	FARP1	AMFR	DHRS2	NXT1	S100A9	SLC29A1	SLC24A3	FOXC1	KIF20A	TOB1	FDFT1	DNAJC1	TPSAB1	TSTA3	FGFR3	SGK1	ID2	GALE	BATF	MAPK13	FABP5	MEST	JAK1	CYP4F11	KCNK5	CPE	XRCC3	CXCL14	SCUBE2	CDC20	IL6ST	GINS2	TRIM29	UNC13B	LAMC2	LARGE	SLC2A8	PLXNB1	PRKAR2B	RPS6KA2	HSPA4L	TFPI2	SERPINA1	TNNC1	HMGCS2	ALDH3A2	CD9	IDH2	SORD	MDK	ALDH3B1	PTGER3	RABEP1	KLF4	PPIF	SNX10	METTL3	PLK4	COX6C	ST14	NCOR2	MOCS2	NMU	TH	RNASEH2A	CHST8	TST	TOP2A	CKB	LTF	DUSP2	PTPN6	ATP2B4	ST6GALNAC2	ADD3	DYNLT3"
#
# HALLMARK_ANDROGEN_RESPONSE
#N = "KLK3	KLK2	ACSL3	PIAS1	CAMKK2	NKX3-1	TMPRSS2	APPBP2	CENPN	BMPR1B	MAF	FADS1	ZBTB10	HMGCR	SPCS3	INSIG1	NGLY1	UBE2J1	ELK4	ABCC4	ELOVL5	ALDH1A3	TARP	AZGP1	ABHD2	SAT1	DBI	DHCR24	SORD	STK39	TPD52	IDI1	B2M	MAP7	DNAJB9	FKBP5	HERC3	PGM3	HOMER2	ELL2	UAP1	SEC24D	LMAN1	PMEPA1	INPP4B	RRP12	SEPP1	NDRG1	KRT19	ITGAV	B4GALT1	MAK	SPDEF	GSR	KRT8	LIFR	PPAP2A	STEAP4	GUCY1A3	SRP19	CCND3	IQGAP2	RAB4A	SLC26A2	SMS	ZMIZ1	SGK1	MERTK	HPGD	ANKH	GNAI3	SRF	SLC38A2	HMGCS1	ARID5B	CDC14B	TNFAIP8	XRCC5	GPD1L	ACTN1	RPS6KA3	TSC22D1	PDLIM5	TMEM50A	ADAMTS1	NCOA4	AKAP12	SCD	PTK2B	CDK6	ADRM1	H1F0	HSD17B14	VAPA	CCND1	XRCC6	MYL12A	UBE2I	PTPN21	AKT1	PA2G4"
#
#
if N :
	N = N.split()
	TCGA = TCGA.ix[N]
	BCXX = BCXX.ix[N]

del N

#print(TCGA)
#print(BCXX)
#exit()


C = TCGA_DATA['subtype']
C = C[ C['tissue_type'] == "tumor" ]

# Check ID uniqueness
assert(len(set(C['aliquot_barcode'])) == len(C['aliquot_barcode']))

# Check that the PAM50 types are represented
PAM50 = ["Normal", "LumA", "LumB", "Her2", "Basal"]
assert(set(PAM50) == set(C['PAM50']))
# Make a new column with PAM50 type ID (0 to 4)
C = C.assign(**{ 'PAM50_ID' : C['PAM50'].apply(lambda x : PAM50.index(x)) })
# Sort by PAM50 type ID
C = C.sort_values(by='PAM50_ID')

# Remove those samples that are not in the merged TCGA table
C = C[ C['aliquot_barcode'].isin(TCGA.columns) ]
assert(set(C['aliquot_barcode']) <= set(TCGA.columns))

## Get the list of barcode
#C = list(C['aliquot_barcode'])

# Allocate matrix
#M = np.zeros( (len(C.index), len(BCXX.columns)) )
M = pd.DataFrame(index=list(C['aliquot_barcode']), columns=list(BCXX.columns))

# Aliquot barcode to PAM50 type
def aliq2pam(a) : return C[ C['aliquot_barcode'] == a ]['PAM50'].item()

# Sample to cluster
def bcxx_s2c(c) : return c[:-3]

# (TCGA PAM50 Type) x (BCXX tumor) matrix
#tt = pd.DataFrame(index=PAM50, columns=set(map(bcxx_s2c, BCXX.columns)))
#print(list(BCXX.columns))
#print(list(tt.columns))
#print(tt); exit()

plt.close('all')
plt.figure(figsize=(10, 5))

for (m, a) in np.random.permutation(list(enumerate(C['aliquot_barcode']))) :
	
	#if not a.startswith(tuple(C['aliquot_barcode'])) : continue
	#if not a.startswith(("TCGA-AN-A0FL", "TCGA-B6-A0RS", "TCGA-E2-A15A")) : continue
	
	M.loc[a, :] = Parallel(n_jobs=PARAM['#proc'])(
		delayed(axb)(TCGA[a], BCXX[b])
		for (n, b) in enumerate(BCXX.columns)
	)
	
	#print(M); exit()
	
	#print(a, t, "x", b)
	
	#I = sorted(range(M.shape[1]), key=(lambda n : -np.sum(M[:, n])), reverse=True)
	#M = M[:, I]
	
	#J = sorted(range(M.shape[0]), key=(lambda m : +np.sum(M[m, :])), reverse=True)
	#M = M[J, :]
	
	#P = np.asarray(M, dtype='f')
	#P = P[ np.isfinite(P.sum(axis=1)), : ]
	
	N = M.dropna().astype(float)
	if ('Normal' in N.index) : N = N.drop(index='Normal')
	N = N.groupby(aliq2pam, axis=0).mean()
	N = N.groupby(bcxx_s2c, axis=1).mean()
	N = N.reindex(pd.Categorical(N.index, categories=PAM50))
	N = N.sort_index()
	
	for b in N.columns :
		N[b] = N[b] / N[b].sum()

	P = np.asarray(N, dtype='f')

	#plt.ion(); plt.show()
	plt.clf()
	plt.imshow(P, aspect='auto', cmap=plt.cm.Blues)
	plt.xticks(range(len(N.columns)), list(N.columns), rotation=60)
	plt.yticks(range(len(N.index)), list(N.index))
	plt.tick_params(axis='x', which='major', labelsize=6)
	plt.tick_params(axis='y', which='major', labelsize=8)
	plt.ylabel("Classified as")
	#plt.colorbar()
	#plt.pause(0.1)
	
	plt.savefig(OFILE['cos-img'])


pickle.dump(
	{
		'M' : M,
		'N' : N,
	},
	open(OFILE['cos-data'], 'wb')
)
