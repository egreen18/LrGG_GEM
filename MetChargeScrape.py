#Importing Data from Matlab
import scipy.io as spio
mat = spio.loadmat('MetParsed.mat', squeeze_me=True)
mets = mat['mets']

#Test to understand how to scrape data from BiGG
import requests
import numpy as np
from bs4 import BeautifulSoup as bs
URLbase="http://bigg.ucsd.edu/universal/metabolites/"
charge=np.array([])
i = 0
while i<len(mets):
	URLopt=mets[i]
	URL=URLbase+URLopt
	webpage=requests.get(URL)
	webcontent = webpage.content
	htmlcontent = bs(webcontent, "html.parser")
	charge=np.append(charge,str(htmlcontent.find("h4",text="Charges in BiGG models: ").find_next_sibling("p")))
	i+=1
charge = charge.tolist()
spio.savemat('ChargeRaw.mat',mdict={'charge': charge})