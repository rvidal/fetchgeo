fetchgeo
========

Downloads gene expression raw data from NCBI's Gene Expression Omnibus (GEO)
Script runs in command line and takes two arguments:
- search query
- number of results to download

	$ python fetchgeo.py 'ovarian cancer' 20

Each result entries are downloaded & stored in a separate folder containing:
- Raw Data
- SOFT file
- Text file containing detailed metadata for each study