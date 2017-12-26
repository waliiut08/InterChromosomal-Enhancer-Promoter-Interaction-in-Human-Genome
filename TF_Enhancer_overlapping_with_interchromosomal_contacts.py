from itertools import islice
from common_functions import read_file, write_file, overlapped, read_data, make_directory
import os, os.path
import pdb


#--------------------------------------------------------------------------------------------------------------------------

def getData_FromFile(filePath, FIRST_CHR, SECOND_CHR): 
	fileData = read_file(filePath)
	fileData = [item.replace("\n","").split("\t") for item in fileData]
	return [item for item in fileData if item[0] == FIRST_CHR or item[0] == SECOND_CHR]

#--------------------------------------------------------------------------------------------------------------------------

def contact(x, y):

	if int(x[0]) > int(y[0]):
		return 1
	return -1
#------------------------------------
def contactSecond(x, y):

	if int(x[1]) > int(y[1]):
		return 1
	return -1


#--------------------------------------------------------------------------------------------------------------------------




def binary_search_Overlapping(EnhancerData, contactData, CHROMOSOME, EOPstart_Indext, EOPend_Index,contactIndex):
	overlapsDetails = []
	contactSize = len(contactData)

	for enhancer in EnhancerData:
		enhancer_start = int(enhancer[EOPstart_Indext])
		enhancer_end   = int(enhancer[EOPend_Index])

		### Binary Search------------------------------------------------------------------------------------------------------------------------------------------------------------
		first = 0
		last = (contactSize-1)
		found = False
		midPoint = None
		contact_chrom_start = None

		while not found and first<= last:
			midPoint = (first + last)//2
			
			contact_chrom_start = int(contactData[midPoint][contactIndex])
			contact_chrom_end   = int(contactData[midPoint][contactIndex]) + 5000
			

			if (CHROMOSOME == enhancer[0] and overlapped(contact_chrom_start, contact_chrom_end, enhancer_start, enhancer_end)):
				overlapsDetails.append(contactData[midPoint] + enhancer)
				found = True

				### computing for different regions of FIRST CHROMOSOME:
				upStream = midPoint - 1
				downStream = midPoint + 1
				
				while upStream>=0 and int(contactData[upStream][contactIndex]) == contact_chrom_start:
					overlapsDetails.append(contactData[upStream] + enhancer)
					upStream -=1

				#################################################################################################################################################################
				### if an enhancer region is common in different regions. ex: enhancer region: 19000-21000 and two sequential contact regions: 15000-20000, 20001-25000. So,  ###
				### this enhancer will be picked as overlapping for both contact regionsself.                                                                                 ###
				#################################################################################################################################################################
				### upstream common region. ex: enhancer: 19000-21000. contact regions: 15000-20000

				if upStream>=0:
					second_chrom_start_prevRegion = int(contactData[upStream][contactIndex])
					second_chrom_end_prevRegion	  = int(contactData[upStream][contactIndex]) + 5000

					if (second_chrom_start_prevRegion != contact_chrom_start) and overlapped(second_chrom_start_prevRegion, second_chrom_end_prevRegion, enhancer_start, enhancer_end):
						
						overlapsDetails.append(contactData[upStream] + enhancer)

						upStream -=1
						while upStream>=0 and int(contactData[upStream][contactIndex]) == second_chrom_start_prevRegion:
							overlapsDetails.append(contactData[upStream] + enhancer)
							upStream -=1


				###-----------------

				while downStream<contactSize and int(contactData[downStream][contactIndex]) == contact_chrom_start :
					overlapsDetails.append(contactData[downStream] + enhancer)
					downStream +=1

				###downstream common region. ex : enhancer: 19000-21000. contact regions: 20001-25000
				if downStream < contactSize:
					second_chrom_start_nextRegion = int(contactData[downStream][contactIndex])
					second_chrom_end_nextRegion	  = int(contactData[downStream][contactIndex]) + 5000

					if (second_chrom_start_nextRegion != contact_chrom_start) and overlapped(second_chrom_start_nextRegion, second_chrom_end_nextRegion, enhancer_start, enhancer_end):
						overlapsDetails.append(contactData[downStream] + enhancer)

						downStream +=1
						while downStream<contactSize and int(contactData[downStream][contactIndex]) == second_chrom_start_nextRegion:
							overlapsDetails.append(contactData[downStream] + enhancer)
							downStream +=1

				###---------


			else:
				if enhancer_start < contact_chrom_start:
					last = midPoint - 1
				else:
					first = midPoint + 1

	return overlapsDetails






#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def get_OverlappingData(EnhancerData, contactFilePath, FIRST_CHR, SECOND_CHR):
	overlapDetails=[]
	first_chrom_overlapDetails = []
	second_chrom_overlapDetails = []
	enhancerStart_Index = 1
	enhancerEnd_Index = 2


	with open(contactFilePath) as contactFL:
	    while True:

	    	### File will read 100000 lines per processing--
	        contactData = list(islice(contactFL, 100000))
	        if not contactData:
	            break

	        contactData = [item.replace("\n","").split("\t") for item in contactData]

	        ################################################ SECOND_CHR ############################################################
	        secondContactIndex = 1
	        contactData_sorted_secondChrom = sorted(contactData, cmp=contactSecond)
	        second_chrom_overlapDetails += binary_search_Overlapping(EnhancerData, contactData_sorted_secondChrom, SECOND_CHR, enhancerStart_Index, enhancerEnd_Index, secondContactIndex)
	        
	        ################################################ FIRST_CHR##############################################################
	        firstContactIndex = 0
	        contactData_sorted_firstChrom = sorted(contactData, cmp=contact)
	        first_chrom_overlapDetails += binary_search_Overlapping(EnhancerData, contactData_sorted_firstChrom, FIRST_CHR,enhancerStart_Index, enhancerEnd_Index, firstContactIndex)

	        ########################################################################################################################


	overlapDetails = second_chrom_overlapDetails + first_chrom_overlapDetails
	return overlapDetails
	
				

#----End of get_OverlappingData()-----------------------------------------------------------------------

