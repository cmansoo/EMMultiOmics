"""
demolist2.txt is a demo of input files.
Running this scripts in the commmand line as: python termClusteringReport.py demonlist2.txt
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
# setCategories
#categorySting = str(client.service.setCategories('BBID,BIOCARTA,COG_ONTOLOGY,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE'))

#getTermClusteringReport
overlap=3
initialSeed = 3
finalSeed = 3
linkage = 0.5
kappa = 50
termClusteringReport = client.service.getTermClusterReport(overlap, initialSeed, finalSeed, linkage, kappa)

#parse and print report
totalRows = len(termClusteringReport)
print ('Total clusters:',totalRows)

resF = 'list1.termClusteringReport.txt'
with open(resF, 'w') as fOut:
    i = 0
    for simpleTermClusterRecord in termClusteringReport:
        i = i+1
        EnrichmentScore = simpleTermClusterRecord.score
        fOut.write('Annotation Cluster '+str(i) + '\tEnrichmentScore:'+str(EnrichmentScore)+'\n')
        fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n')
        for simpleChartRecord in  simpleTermClusterRecord.simpleChartRecords:
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

print ('write file:', resF, 'successfully finished!')