The files

	BRCA.datafreeze.20120227.txt
	BRCA.547.PAM50.SigClust.Subtypes.txt

were downloaded from

	https://tcga-data.nci.nih.gov/docs/publications/brca_2012/



The meta files 

	gdc_manifest.2018-01-24T03_39_35.725692.txt
	files.2018-01-24T04_54_05.741627.json

were downloaded from

	https://portal.gdc.cancer.gov/repository?facetTab=files&filters={"op":"and","content":[{"op":"in","content":{"field":"cases.project.project_id","value":["TCGA-BRCA"]}},{"op":"in","content":{"field":"files.access","value":["open"]}},{"op":"in","content":{"field":"files.data_category","value":["Biospecimen","Clinical","Copy Number Variation","Transcriptome Profiling"]}}]}

Data selection is as follows

	Project Id IS TCGA-BRCA
	AND 
	Access IS open
	AND
	Data Category IN ( Biospecimen, Clinical, Copy Number Variation, Transcriptome Profiling )


Data files can be downloaded with

	cd UV; ./gdc-client download -m ../gdc_manifest.2018-01-24T03_39_35.725692.txt

(unless an error occurrs) or the a_download.py script.


For info about gdc-client see

	https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/


More about TCGA-BRCA:

	https://tcga-data.nci.nih.gov/docs/publications/brca_2012/
	https://wiki.cancerimagingarchive.net/display/Public/TCGA-BRCA
	https://wiki.nci.nih.gov/display/TCGA/File+Format+Specifications
	https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/


The file

	transcriptome_list.tsv

was obtained on 2018-01-28 using the command

	curl --request POST --header "Content-Type: application/json" --data @transcriptome_list_payload.json 'https://api.gdc.cancer.gov/files' > transcriptome_list.tsv




RA, 2017-01-24
