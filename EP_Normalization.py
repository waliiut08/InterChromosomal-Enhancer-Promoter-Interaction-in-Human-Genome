
from common_functions import read_file, write_file, overlapped, read_data
from itertools import islice
import pdb

#---------------------------------------------------------------------------



def normalizing(EP_interactionData, FIRST_CHR_NORMpath, SECOND_CHR_NORMpath):
	
	first_chr_NormData = read_file(FIRST_CHR_NORMpath)
	first_chr_NormData = [item.replace("\n","") for item in first_chr_NormData]

	second_chr_NormData = read_file(SECOND_CHR_NORMpath)
	second_chr_NormData = [item.replace("\n","") for item in second_chr_NormData]
	
	
	for item in EP_interactionData:

		chr1_NormLine = int(item[0])/5000 +1
		chr2_NormLine =	int(item[1])/5000 +1

		try:
			first_chr_NormData[chr1_NormLine-1] and second_chr_NormData[chr2_NormLine-1]
		except:
			pdb.set_trace()




		if first_chr_NormData[chr1_NormLine-1] and second_chr_NormData[chr2_NormLine-1]:
			try:
				norm = float(item[2])/(float(first_chr_NormData[chr1_NormLine-1])*float(second_chr_NormData[chr2_NormLine-1]))
				
			except:
				print str(chr1_NormLine) + " " + str(chr2_NormLine) + ": " + str(first_chr_NormData[chr1_NormLine-1]) + " " + str(second_chr_NormData[chr2_NormLine-1])
				norm = "None"
		else:
			norm = "None"

		item.append(str(norm))

	
	EP_InteractionsDetails = "\n".join(["\t".join(item) for item in EP_interactionData])
	return EP_InteractionsDetails

