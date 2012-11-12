# fetchgeo #
========

## Requirements ##
* Python 2.7+
* Biopython
* urllib2

## Description ##
Downloads gene expression raw data from NCBI's Gene Expression Omnibus (GEO)
Script runs in command line and takes two arguments:
- search query
- number of results to download

## Usage ##
	$ python fetchgeo.py 'ovarian cancer' 20

## Output ##
Each result entries are downloaded & stored in a separately labelled folder containing:
- Raw Data
- SOFT file
- Text file containing detailed metadata for each study