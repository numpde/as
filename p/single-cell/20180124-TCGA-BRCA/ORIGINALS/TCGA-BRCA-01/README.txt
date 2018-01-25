The meta files 

	gdc_manifest.2018-01-24T03_39_35.725692.txt
	files.2018-01-24T04_54_05.741627.json

were downloaded from

	https://portal.gdc.cancer.gov/repository?facetTab=files&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-BRCA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_category%22%2C%22value%22%3A%5B%22Biospecimen%22%2C%22Clinical%22%2C%22Copy%20Number%20Variation%22%5D%7D%7D%5D%7D&searchTableTab=files


Data selection is as follows

	Project Id IS TCGA-BRCA
	AND 
	Access IS open
	AND
	Data Category IN ( Biospecimen, Clinical, Copy Number Variation, Transcriptome Profiling )


Data files downloaded with

	cd UV; ./gdc-client download -m ../gdc_manifest.2018-01-24T03_39_35.725692.txt


For info about gdc-client see

	https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/


More about TCGA-BRCA:

	https://tcga-data.nci.nih.gov/docs/publications/brca_2012/
	https://wiki.cancerimagingarchive.net/display/Public/TCGA-BRCA
	https://wiki.nci.nih.gov/display/TCGA/File+Format+Specifications



RA, 2017-01-24
