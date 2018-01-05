
# RA, 2017-12-25

import networkx
import pickle

G = pickle.load(
	open(
		"../OUTPUT/0_go-graph/UV/go-graph.pkl",
		'rb'
	)
)

GO2T = {
	go : data['name']
	for (go, data) in G.nodes(data=True)
}

print(list(GO2T.items())[0:20])


top50 = [
	"spindle pole",
	"sister chromatid cohesion",
	"DNA replication",
	"endomembrane system",
	"centriole",
	"guanyl-nucleotide exchange factor activity",
	"peptidyl-serine phosphorylation",
	"microtubule cytoskeleton",
	"Rab GTPase binding",
	"microtubule organizing center",
]

last50 = [
	"multicellular organism development",
	"protein ubiquitination",
	"protein kinase binding",
	"viral process",
	"cadherin binding",
	"regulation of transcription from RNA polymerase II promoter",
	"intracellular",
	"cell differentiation",
	"ubiquitin protein ligase binding",
	"protein heterodimerization activity",
	"protein complex",
	"cell adhesion",
	"positive regulation of cell proliferation",
	"negative regulation of transcription, DNA-templated",
	"negative regulation of apoptotic process",
	"enzyme binding",
	"cell proliferation",
	"mitochondrial inner membrane",
	"neutrophil degranulation",
	"response to drug",
]

print("# Top 50 from TXP (2017.01.04)")
for name in top50 :
	go = [go for (go, t) in GO2T.items() if (t == name)][0]
	print('"{}", # {}'.format(go, name))

print("# Last 50 from TXP (2017.01.04)")
for name in last50 :
	go = [go for (go, t) in GO2T.items() if (t == name)][0]
	print('"{}", # {}'.format(go, name))

