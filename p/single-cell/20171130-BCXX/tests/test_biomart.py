#https://pypi.python.org/pypi/biomart/0.8.0
from biomart import BiomartServer

server = BiomartServer( "http://grch37.ensembl.org/biomart" )

# set verbose to True to get some messages
#server.verbose = True

## show server databases
#server.show_databases() # uses pprint behind the scenes

## show server datasets
#server.show_datasets() # uses pprint behind the scenes

# use the 'hsapiens_gene_ensembl' dataset
human = server.datasets['hsapiens_gene_ensembl']

#human.show_filters()
#input("Press enter...")
#human.show_attributes()
#input("Press enter...")

# 
response = human.search({
  'filters': {
      #'hgnc_symbol': 'TAOK3'
      "ensembl_gene_id" : "ENSG00000135090"
  },
  'attributes': [
      'go_id'
  ]
}, header = 0 )

print(response.content.decode("utf-8").split())
