"""
demolist2.txt is a demo of input files.
Running this scripts in the commmand line as: python DAVIDWebservice_Client.py demonlist2.txt 
"""

import pandas as pd
import sys
from suds.client import Client

url = 'https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl'
    
print ('url=%s' % url)

#
# create a service client using the wsdl.
#
client = Client(url)
client.wsdl.services[0].setlocation('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')
#
# print the service (introspection)
#
#print client

#authenticate user email 
# client.service.authenticate('yourRegisteredEmail@your.org')
client.service.autheticate("mcho8@bwh.harvard.edu")

# Define the input gene list file, identifier type, list type and a list name
inputFile = sys.argv[1]
idType = 'AFFYMETRIX_3PRIME_IVT_ID'
listType = 0
listName = "demolist2"

# Read input gene list file, convert ids to a comma-delimited string and upload the list to DAVID
df=pd.read_csv(inputFile, usecols=[0], delimiter='\t', index_col=False, names=['AFFYMETRIX_3PRIME_IVT_ID'])
inputIds = ",".join(df['AFFYMETRIX_3PRIME_IVT_ID'].astype(str).unique().tolist())
client.service.addList(inputIds, idType, listName, listType)
print (client.service.addList(inputIds, idType, listName, listType))

print (client.service.getDefaultCategoryNames())
print ('successfully finished!')
# setCategories
#categorySting = str(client.service.setCategories('BBID,BIOCARTA,COG_ONTOLOGY,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE'))
