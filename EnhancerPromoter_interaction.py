from common_functions import read_file,write_file,overlapped
from itertools import islice
import pdb



#--------------------------------------------------------------------------------------
def contact(x, y):

	if int(x[0]) > int(y[0]):
		return 1
	return -1


#---------------------------------------------------------------------------------------
def enhancer_promoter_interaction(EnhPromData, FIRST_CHR, SECOND_CHR):

	interactionsDetails = []
	
	for eachrow in EnhPromData:
		
		if(eachrow[3]!=eachrow[6]):
			#interactionsFound=interactionsFound+1
			interactionsDetails.append(eachrow)
	
	interactionsDetails = sorted(interactionsDetails, cmp=contact)

	#print "Total interactions found: " + str(interactionsFound)

	return interactionsDetails

#-----End of enhancer_promoter_interaction()-----------------------------------------------


