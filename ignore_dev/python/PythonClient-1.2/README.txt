DAVID Web Service client using Python
=============================================
The sampye client files contain Python source code showing you how to connect DAVID Web service 
and generate your reports. 

Prerequisites
==============
Register your email address in https://david.ncifcrf.gov/webservice/register.htm
Download and install suds-0.4 (orlater version) from https://fedorahosted.org/suds/

(https://fedorahosted.org/suds/ no longer available, so download pip install suds)

('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')
try changing this to SOAP12 endpoint?


Following the steps:

1. Download  PythonClient-1.1.zip from DAVID WebService site:
   https://davidbioinformatics.nih.gov/content.jsp?file=WS.html

2. Extract PythonClient-1.1.zip to your working directory: $PythonClient

3. go to your directory $PythonClient
   
   Edit DAVIDWebService_Client.py and run it to get you connected to the service 
  
4. check if the following web link returns a 'true', which means the server is ready for conection.
   https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService/authenticate?args0=YourRegisteredEmail@your.org

5. Run your client file and generate your own report text file 
