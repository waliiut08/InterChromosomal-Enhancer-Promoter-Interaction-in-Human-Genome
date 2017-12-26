from itertools import islice
from common_functions import read_file, write_file, overlapped, read_data, make_directory



#----------------------------------------------------------------------------------------------
def binary_search_Overlapping(promoterData, contactData, CHROMOSOME, EOPstart_Indext, EOPend_Index,contactIndex):
	overlapsDetails = []
	contactSize = len(contactData)

	for promoter in promoterData:
		promoter_start = int(promoter[EOPstart_Indext])
		promoter_end   = int(promoter[EOPend_Index])

		### Binary Search--------------------------
		first = 0
		last = (contactSize-1)
		found = False
		midPoint = None
		contact_chrom_start = None

		while not found and first<= last:
			midPoint = (first + last)//2
			
			contact_chrom_start = int(contactData[midPoint][contactIndex])
			contact_chrom_end   = int(contactData[midPoint][contactIndex]) + 5000
			

			if (CHROMOSOME == promoter[0] and overlapped(contact_chrom_start, contact_chrom_end, promoter_start, promoter_end)):
				overlapsDetails.append(contactData[midPoint] + promoter)
				found = True

				### computing for different regions of FIRST CHROMOSOME:
				upStream = midPoint - 1
				downStream = midPoint + 1
				
				while upStream>=0 and int(contactData[upStream][contactIndex]) == contact_chrom_start:
					overlapsDetails.append(contactData[upStream] + promoter)
					upStream -=1

				###################
				### if a promoter region is common in different regions. ex: promoter region: 19000-21000 and two sequential contact regions: 15000-20000, 20001-25000. So, 
				### this promoter will be picked as overlapping for both contact regions.
				##################
				### upstream common region. ex: promoter: 19000-21000. contact regions: 15000-20000

				if upStream>=0:
					second_chrom_start_prevRegion = int(contactData[upStream][contactIndex])
					second_chrom_end_prevRegion	  = int(contactData[upStream][contactIndex]) + 5000

					if (second_chrom_start_prevRegion != contact_chrom_start) and overlapped(second_chrom_start_prevRegion, second_chrom_end_prevRegion, promoter_start, promoter_end):
						
						overlapsDetails.append(contactData[upStream] + promoter)

						upStream -=1
						while upStream>=0 and int(contactData[upStream][contactIndex]) == second_chrom_start_prevRegion:
							overlapsDetails.append(contactData[upStream] + promoter)
							upStream -=1


				###-----------------

				while downStream<contactSize and int(contactData[downStream][contactIndex]) == contact_chrom_start :
					overlapsDetails.append(contactData[downStream] + promoter)
					downStream +=1

				###downstream common region. ex : promoter: 19000-21000. contact regions: 20001-25000
				if downStream < contactSize:
					second_chrom_start_nextRegion = int(contactData[downStream][contactIndex])
					second_chrom_end_nextRegion	  = int(contactData[downStream][contactIndex]) + 5000

					if (second_chrom_start_nextRegion != contact_chrom_start) and overlapped(second_chrom_start_nextRegion, second_chrom_end_nextRegion, promoter_start, promoter_end):
						overlapsDetails.append(contactData[downStream] + promoter)

						downStream +=1
						while downStream<contactSize and int(contactData[downStream][contactIndex]) == second_chrom_start_nextRegion:
							overlapsDetails.append(contactData[downStream] + promoter)
							downStream +=1

				###---------


			else:
				if promoter_start < contact_chrom_start:
					last = midPoint - 1
				else:
					first = midPoint + 1

	return overlapsDetails

#--------------------------------------------------------------------------------------------------------------------------

def contact_second(x, y):

	if int(x[1]) > int(y[1]):
		return 1
	return -1

#--------------------------------------------------------------------------------------------------------------------------


def contact_first(x, y):

	if int(x[0]) > int(y[0]):
		return 1
	return -1

#--------------------------------------------------------------------------------------------------------------------------





def get_Promoter_OverlappingData(promoterData, enhancer_contactData, FIRST_CHR, SECOND_CHR):

	overlapDetails=[]
	first_chrom_overlapDetails = []
	second_chrom_overlapDetails = []
	promoterStart_index = 5
	promoterEnd_index = 6



	################################################ SECOND_CHR ############################################################################################
	second_chromIndex = 1
	contactData_sorted_secondChrom = sorted(enhancer_contactData, cmp = contact_second)
	second_chrom_overlapDetails += binary_search_Overlapping(promoterData, contactData_sorted_secondChrom, SECOND_CHR, promoterStart_index, promoterEnd_index, second_chromIndex)

	################################################ FIRST_CHR################################################################################################
	first_chromIndex = 0
	contactData_sorted_firstChrom = sorted(enhancer_contactData, cmp = contact_first)
	first_chrom_overlapDetails += binary_search_Overlapping(promoterData, contactData_sorted_firstChrom, FIRST_CHR, promoterStart_index, promoterEnd_index, first_chromIndex)

	##########################################################################################################################################################
	
	overlapDetails = second_chrom_overlapDetails + first_chrom_overlapDetails
	return overlapDetails

#--------------------------------------------------------------------------------------------------------------------------------------
	

