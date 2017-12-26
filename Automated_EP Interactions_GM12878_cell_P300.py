from common_functions import write_file, find_between, make_directory, elapsed_time, expected_finishing_time
from TF_Enhancer_overlapping_with_interchromosomal_contacts import getData_FromFile, get_OverlappingData
from EnhancerPromoter_interaction import enhancer_promoter_interaction
from Promoter_overlapping_with_interchromosomal_contacts import get_Promoter_OverlappingData
from EP_Normalization import normalizing
import time
import sys
import time
import os
import pdb
#---------------------------------------------------------------------------------------------------------------------
StartTime = time.time()

 
interchromosomalContactLoc = '/home/muhammad/PhD/Papers/Li(Nov16-2015)/Project(Enhancers)/Scripts/Efficient Search/GM 12878 Cell/Interchromosomal Contact data/GM12878_primary_interchromosomal/5kb_resolution_interchromosomal/'
contactPaths = []
contactFiles = []

## getting the contact file location and their corresponding contact and normalization files--------
i =0
j =0
for dirpath,dirname,filename in os.walk(interchromosomalContactLoc):	
	if (i+1)%3==0:
		contactPaths.append(dirpath)
		contactFiles.append(filename)
	i+=1
	j+=1

##-----------------------------------


total_contactDirectory = len(contactPaths)
print "'Cell: GM-12878 & TF: P300'. Number of interchromosomal contact direcotries are: " + str(total_contactDirectory) + "\n"

totalRequiredTime = total_contactDirectory ### as each directory requires approximately 1 minute to complete the whole computation result. 
print "Program started at:- " + str((time.strftime("%H:%M:%S")))
expected_finishing_time(totalRequiredTime)



for k in range(total_contactDirectory):

	iterative_StartTime = time.time()
	interactingChromosomes = find_between(contactPaths[k], "5kb_resolution_interchromosomal/" , "/MAPQG0")

	### Getting the interacting chromosome names. ex: FIRST_CHR = chr1, SECOND_CHR = chr2. ----------------------

	FIRST_CHR = interactingChromosomes[0:interactingChromosomes.index("_")]
	SECOND_CHR = interactingChromosomes[interactingChromosomes.index("_")+1:len(interactingChromosomes)]
	
	outputFile = "EP_Interaction_for_all_chromosome_GM-12878_P300/"+ FIRST_CHR + "_"+ SECOND_CHR +"/Enhancer-Promoter Interactions and normalized vector between " + FIRST_CHR + "_" + SECOND_CHR + ".txt"
	make_directory(outputFile)
	
	
	### Concatenating the contact directory path and Rawobserved file ----------------------------------------------

	rawContactFile = filter(lambda x: 'RAWobserved' in x, contactFiles[k])
	if rawContactFile==[]:
		print "No RAWobserved file for" + FIRST_CHR + "_" + SECOND_CHR
		continue

	rawContactFilePath = contactPaths[k] +"/" + rawContactFile[0]


	### Getting .KRnorm files path ---------------------------------------------------------------------------------
	normFiles = filter(lambda x: ".KRnorm" in x, contactFiles[k])
	first_chr_normFile = filter(lambda x: FIRST_CHR+"_" in x, normFiles)
	second_chr_normFile = filter(lambda x: SECOND_CHR+"_" in x, normFiles)
	
	FIRST_CHR_NORMpath = contactPaths[k] + "/"+ first_chr_normFile[0]
	SECOND_CHR_NORMpath = contactPaths[k] + "/"+ second_chr_normFile[0]
	
	print str(k+1) +". Calculating the Enhancer-Promoter interactions between: " +  FIRST_CHR + "_" + SECOND_CHR + "........... "
	#print "\n"+FIRST_CHR + ": norm file name  and path - " + FIRST_CHR_NORMpath
	#print SECOND_CHR + ": norm file name  and path - " + SECOND_CHR_NORMpath + "\n"

	


	################################ Main Computation Starts Here #################################################################################################################

	### TF(Enhancer) overlapping with Interchromosomal Contacts ---------------------------------------------------------------------------------------------------------------------------------

	print "step-1: Computing Enhancer_Contact Overlapping.."
	
	
	TF_EnhancerPath   = '/home/muhammad/PhD/Papers/Li(Nov16-2015)/Project(Enhancers)/Scripts/Efficient Search/GM 12878 Cell/GM 12878 Cell & P300 TF/Data/GSM935294_hg19_wgEncodeSydhTfbsGm12878P300IggmusPk.narrowPeak'
	TF_EnhancerData   = getData_FromFile(TF_EnhancerPath, FIRST_CHR, SECOND_CHR)
	
	for TF in TF_EnhancerData:
		del TF[3:10]
	
	Enhancer_Contact_OverlappedData = get_OverlappingData(TF_EnhancerData, rawContactFilePath, FIRST_CHR, SECOND_CHR)
	
	


	### Promoter Overlapping with Interchromosomal Contacts with Enhancer overlapping----------------------------------------------------------------------------------------------------
	print "step-2: Computing Promoter_Contact Overlapping.."
	
	
	promoterPath = '/home/muhammad/PhD/Papers/Li(Nov16-2015)/Project(Enhancers)/Scripts/Efficient Search/GM 12878 Cell/GM 12878 Cell & P300 TF/Data/promoter_regions(1k upstream of the Gene)_hg19_UCSC.txt' 
	promoterData = getData_FromFile(promoterPath, FIRST_CHR, SECOND_CHR)

	Promoter_Enhancer_ContactOverlappedData = get_Promoter_OverlappingData(promoterData, Enhancer_Contact_OverlappedData, FIRST_CHR, SECOND_CHR)

	
	### Enhancer-Promoter Interactions -------------------------------------------------------------------------------------------------------------------------------------------------
	print "step-3: Computing EP overlapping.."
	
	EP_interactionData = enhancer_promoter_interaction(Promoter_Enhancer_ContactOverlappedData , FIRST_CHR, SECOND_CHR)
	

	### Normalization --------------------------------------------------------------------------------------------------------------------------------------------------------------------
	print "step-4: Normalizing...."
	
	EP_normalizedData = normalizing(EP_interactionData, FIRST_CHR_NORMpath, SECOND_CHR_NORMpath)
	write_file(EP_normalizedData, outputFile, "w")
        #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
	iterative_endTime = time.time()
	elapsed_time(iterative_StartTime, iterative_endTime)

	print "Computation of EP interactions between " + FIRST_CHR + "_"+ SECOND_CHR + " is completed. Thanks! \n"
	
	
	
print "\nComputation of 'GM-12878 cell & P300 TF' has been completed."
EndTime = time.time()
elapsed_time(StartTime, EndTime)
