"""
demolist2.txt is a demo of input files.
Running this scripts in the commmand line as: python chartReport.py demonlist2.txt
"""

import pandas as pd
import sys
from suds.client import Client

# url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
url = 'https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl'
   
print ('url=%s' % url)

#
# create a service client using the wsdl.
#
client = Client(url)
client.wsdl.services[0].setlocation('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')
# client.wsdl.services[0].setlocation('https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')

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

# setCategories
#categorySting = str(client.service.setCategories('BBID,BIOCARTA,COG_ONTOLOGY,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE'))
#getChartReport
thd = 0.1
ct = 2
chartReport = client.service.getChartReport(thd,ct)
chartRow = len(chartReport)
print ('Total chart records:',chartRow)

#parse and print chartReport
resF = 'list1.chartReport.txt'
with open(resF, 'w') as fOut:
	fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n')
	for simpleChartRecord in chartReport:
            categoryName = simpleChartRecord.categoryName
            termName = simpleChartRecord.termName
            listHits = simpleChartRecord.listHits
            percent = simpleChartRecord.percent
            ease = simpleChartRecord.ease
            Genes = simpleChartRecord.geneIds
            listTotals = simpleChartRecord.listTotals
            popHits = simpleChartRecord.popHits
            popTotals = simpleChartRecord.popTotals
            foldEnrichment = simpleChartRecord.foldEnrichment
            bonferroni = simpleChartRecord.bonferroni
            benjamini = simpleChartRecord.benjamini
            FDR = simpleChartRecord.afdr
            rowList = [categoryName,termName,str(listHits),str(percent),str(ease),Genes,str(listTotals),str(popHits),str(popTotals),str(foldEnrichment),str(bonferroni),str(benjamini),str(FDR)]
            fOut.write('\t'.join(rowList)+'\n')

print ('write file:', resF, 'finished!')