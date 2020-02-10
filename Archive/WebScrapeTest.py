#Importing Data from Matlab
import scipy.io as spio
mat = spio.loadmat('MetParsed.mat', squeeze_me=True)
mets = mat['mets']
print(mets[0])

#Test to understand how to scrape data from BiGG
import requests
from bs4 import BeautifulSoup as bs
URLbase="http://bigg.ucsd.edu/universal/metabolites/"
URLopt=mets[0]
URL=URLbase+URLopt
webpage=requests.get(URL)
webcontent = webpage.content
htmlcontent = bs(webcontent, "html.parser")
charge = htmlcontent.find("h4",text="Charges in BiGG models: ").find_next_sibling("p")
print(charge)
