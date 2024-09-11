"""
demolist2.txt is a demo of input files.
Running this scripts in the commmand line as: python tableReport.py demonlist2.txt
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
client.service.authenticate('yourRegisteredEmail@your.org')

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

#print client.service.getDefaultCategoryNames()
categorySting = str(client.service.setCategories('BBID,BIOCARTA,COG_ONTOLOGY,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE'))
categories=categorySting.split(',')

tableReport = client.service.getTableReport()
tableRow = len(tableReport)
print ('Total table records:',tableRow)

resF = 'list1.tableReport.txt'
with open(resF, 'w') as fOut:
    categoryConcat = '\t'.join(categories);
    fOut.write('ID\tGene Name\tSpecies\t'+categoryConcat)
    for tableRecord in tableReport:
        name = tableRecord.name
        species = tableRecord.species
        for arrayString in tableRecord.values:
            gene_id = ','.join(x for x in arrayString.array)
            rowList = [gene_id,name,species]
            fOut.write('\n'+'\t'.join(rowList))
        for annotationRecord in tableRecord.annotationRecords:
            default_value = ''
            category_dict = dict.fromkeys(categories,default_value)
            termsConcat = '';
            for term in annotationRecord.terms:
                termString = term.split("$")[1]
                termList = [termString,termsConcat]
                termsConcat = ','.join(termList)
                category_dict[str(annotationRecord.category)] = termsConcat;
            for key in category_dict:
                fOut.write('\t'+category_dict[key])	

print ('write file:', resF, 'finished!')