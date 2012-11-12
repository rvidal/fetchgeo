#!/usr/bin/env python
# encoding: utf-8
"""
Downloads gene expression raw data from NCBI's Gene Expression Omnibus (GEO)
Script runs in command line and takes two arguments:
	- search query
	- number of results to download

Ex: #> python fetchgeo.py 'ovarian cancer' 20

Each result entries are downloaded & stored in a separate folder containing:
	- Raw Data
	- SOFT file
	- Text file containing detailed metadata for each study
"""

import sys
import os
from Bio import Entrez, Medline
from urllib2 import Request, urlopen, URLError
import time
import gzip


SEARCH_QUERY = sys.argv[1]
MAX_RESULTS = sys.argv[2]

EMAIL_ADDR = "EMAIL@HERE.COM"
HOMEDIR = os.getcwd()

def get_geo_data(s_terms, email, maxresults=1000):
	"""Get geo data for query term provided"""
	Entrez.email = email
	Entrez.tool = 'GetGeoMetaDataPythonScript'
	hits = []
	s_terms = s_terms + ' NOT GPL'
	handle = Entrez.esearch(db='gds', term=s_terms, retmax=maxresults, usehistory="y")
	results = Entrez.read(handle)
	newhandle = Entrez.esummary(db='gds', retmax = maxresults, webenv=results['WebEnv'], query_key=results['QueryKey'])
	summary = Entrez.read(newhandle)
	for i in range(len(summary)):
		samples = []
		geo_data = {'id':summary[i]['Id'],
					'n_samples':summary[i]['n_samples'],
					'pubdate':summary[i]['PDAT'],
					'platform':summary[i]['PlatformTitle'],
					'suppfile':summary[i]['suppFile'],
					'taxon':summary[i]['taxon'],
					'entry_type':summary[i]['entryType'],
					'gpl':summary[i]['GPL'],
					'gse':summary[i]['GSE'],
					'pubmed_ids':summary[i]['PubMedIds'],
					'title':summary[i]['title'],
					'gds_type':summary[i]['gdsType'],
					'summary':summary[i]['summary'],
					'soft_file':get_soft_url(summary[i]['GSE'], summary[i]['entryType'])
					}
		hits.append(geo_data)
	return hits


def get_pubmed_data(idlist):
	"""Takes a list of pubmed ids and returns title, auth, yr"""
	handle = Entrez.efetch(db='pubmed', id=idlist, rettype='medline', retmode='text')
	records = Medline.parse(handle)
	mypms = []
	for record in records:
		mypms.append((record["TI"], record["AU"], record["PMID"]))
	return mypms

def get_soft_url(file_id, datatype):
	"""Get soft file depending on file type and ID"""
	if datatype == 'GSE':
		soft_file = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/" + datatype + file_id + "/" + datatype + file_id + "_family.soft.gz"
	else:
		soft_file = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/" + "GSE" + file_id + "/" + "GSE" + file_id + "_family.soft.gz"
		#soft_file = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/" + datatype + "/" + datatype + file_id + ".soft.gz"
	return soft_file


def parse_soft(softzip):
	""" Take a SOFT gz file and parse it for specific data parse_soft
		Dictionary with platform and available suppfile(s) when available
		Type: results['platform'] (string)
		Type: results['suppfiles'] (list)
	"""
	z = gzip.open(softzip, 'rb')
	results = {}
	suppfiles = []
	platform = ''
	for line in z:
		# Get platform title
		if line.startswith('!Platform_title'):
			index = line.find('= ')
			platform = line[index+2:]
			continue
		else:
			pass
		# Get link to supplementary data file(s)
		if line.startswith('!Series_supplementary_file'):
			index = line.find('= ')
			suppfiles.append(line[index+2:].strip())
			continue
		else:
			pass
		# If platform and suppfiles found, stop reading (long!) file.
		if platform != '' and suppfiles != []:
			break
		else:
			continue
	results['platform'] = platform.strip()
	results['suppfiles'] = suppfiles	
	return results


def urlretrieve(urlfile, fpath):
	req = Request(urlfile)
	chunk = 4096
	try:
		rfile = urlopen(urlfile)
	except IOError, e:
		if hasattr(e, 'reason'):
			print 'Can\'t seem to reach the server.'
			print 'Reason: ', e.reason
		elif hasattr(e, 'code'):
			print 'The server couldn\'t fulfill the request.'
			print 'Error code: ', e.code
		return False
	else:
		print "Began downloading", fpath
		f = open(fpath, 'w')
		while 1:
			data = rfile.read(chunk)
			if not data:
				logfilepath = HOMEDIR + "/log.txt"
				logfile = open(logfilepath, 'a')
				donemsg = "Downloaded " + fpath + "\n"
				logfile.write(donemsg)
				logfile.close()
				print donemsg
				break
			f.write(data)
			#print "Read %s bytes"%len(data)


def dl_geo(data):
	"""Download geo data and log info in the following format
	"""
	for item in data:
		if item['suppfile'] != "":
			raw_data = 'ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/series/GSE' + str(item['gse']) + '/GSE' + str(item['gse']) + '_RAW.tar'
		else:
			raw_data = ''
		
		uniqueid =  "GSE" + str(item['gse'])
		newdir = HOMEDIR + "/" + uniqueid
		
		""" Check to see if the directory already exists.
			If so, skip to next item.
		"""
		if os.path.isdir(newdir) == True:
			print "Skipped " + uniqueid
			continue
		else:
			os.makedirs(newdir) # make new directory for sample
			os.chdir(newdir) # change to the created directory
			""" Download the SOFT file for series
			"""
			save_soft = uniqueid + "_family.soft.gz"
			if urlretrieve(item['soft_file'], save_soft) == False:
				print "No " + str(save_soft) + " file available."
				break
			else:
		
				""" Parse the downloaded SOFT file for platform & possible suppfile(s)
					Dictionary returned with platform and suppfiles (when available.)
				"""
				soft_res = parse_soft(save_soft)
		
				""" Download the raw file if available
					First attempts with constructed url, then falls back on SOFT file.
					If then not available, nothing is stored.
				"""
				if raw_data != '': # if generated _RAW.tar.gz generated (CEL, TXT, avail)
					save_raw = uniqueid + "_RAW.tar" # create local filename
					if urlretrieve(raw_data, save_raw) == False: # if dl attempt fails try soft for urls
						for supf in soft_res['suppfiles']: # check for multiple suppfile urls
							if supf != '': # verify if supplementary files entry is not empty or not available
								if urlretrieve(supf, os.path.basename(supf)) == False: # if dl attempt fails give up
									pass
								else:
									urlretrieve(supf, os.path.basename(supf)) # download the file, finally!
							else:
								pass
					else:
						pass
				else:
					pass
			
			""" Write text file with related information about the current
				downloaded series. Save to filename with GSEID (ex. GSE0000.txt)
			"""
			if item['platform'] == '':
				platf_name = soft_res['platform']
			else:
				platf_name = item['platform']
			geo_results = "GSE ID: " + str(item['gse']) + '\n' + \
						"PubDate: " + str(item['pubdate']) + '\n' + \
						"Number of Samples: " + str(item['n_samples']) + '\n' + \
						"Title: " + item['title'] + '\n' + \
						"Taxonomy: " + item['taxon'] + '\n' + \
						"PubMed: " + str(item['pubmed_ids']) + '\n' + \
						"SOFT downloaded from: " + item['soft_file'] + '\n' + \
						"Platform: " + platf_name + '\n' + \
						"Raw Data downloaded from: " + raw_data + "\n"
			filename = uniqueid + ".txt"
			f = open(filename, 'w')
			f.write(geo_results)
			f.close()
			os.chdir(HOMEDIR)

#run the query and pulling the data
dl_geo(get_geo_data(SEARCH_QUERY, EMAIL_ADDR, MAX_RESULTS))