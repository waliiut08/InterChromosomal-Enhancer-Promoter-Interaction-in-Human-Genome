
from common_functions import write_file, find_between, make_directory, elapsed_time, expected_finishing_time, read_file
import pdb
import operator
import os



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def getData_FromFile(filePath, enhancerChrom): 
	fileData = read_file(filePath)
	fileData = [item.replace("\n","").split("\t") for item in fileData]
	return [item for item in fileData if item[3] == enhancerChrom]

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def contact(x, y):

	if int(x[4]) > int(y[4]):
		return 1
	return -1

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


def getting_EP_InteractionFiles(EP_interactionLoc, EnhancerChromosome):
	i =0
	matchedFiles=[]
	for dirpath,dirname,filename in os.walk(EP_interactionLoc):	
		if i>0:
			interactingChromosomes = find_between(filename[0], "Enhancer-Promoter Interactions and normalized vector between " , ".txt")
			interactingChromosomes = interactingChromosomes.split("_")
			
			if interactingChromosomes[0]==EnhancerChromosome or interactingChromosomes[1]==EnhancerChromosome:
				matchedFiles.append(dirpath +"/" + filename[0])

		i+=1
	print  EnhancerChromosome + ": Matched files - " + str(len(matchedFiles))
	return matchedFiles
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



def binary_search(enhancer, EP_interactedData_sorted):


	EP_size = len(EP_interactedData_sorted)
	first = 0
	last = (EP_size-1)
	found = False
	midPoint = None
	interactedChrom = []
	interactedDetails = []

	while not found and first<= last:
		midPoint = (first + last)//2
		if int(enhancer[1]) == int(EP_interactedData_sorted[midPoint][4]):
			found = True
			interactedChrom.append(EP_interactedData_sorted[midPoint][6] + "_" + EP_interactedData_sorted[midPoint][11] + "_" + EP_interactedData_sorted[midPoint][12])
			interactedDetails.append(enhancer + EP_interactedData_sorted[midPoint])


			upStream = midPoint - 1
			downStream = midPoint + 1

			while upStream>=0 and int(enhancer[1]) == int(EP_interactedData_sorted[upStream][4]):
				interactedChrom.append(EP_interactedData_sorted[upStream][6] + "_" + EP_interactedData_sorted[upStream][11] + "_" + EP_interactedData_sorted[upStream][12])
				interactedDetails.append(enhancer + EP_interactedData_sorted[upStream])
				upStream-=1

			while downStream<=EP_size and int(enhancer[1]) == int(EP_interactedData_sorted[downStream][4]):
				interactedChrom.append(EP_interactedData_sorted[downStream][6] + "_" + EP_interactedData_sorted[downStream][11] + "_" + EP_interactedData_sorted[downStream][12])
				interactedDetails.append(enhancer + EP_interactedData_sorted[downStream])
				downStream+=1


		else:
			if int(enhancer[1]) < int(EP_interactedData_sorted[midPoint][4]):
				last = midPoint - 1


			else:
				first = midPoint + 1

	interactedChrom = ",".join(interactedChrom)
	enhancerDetails = enhancer[0] +"\t" + enhancer[1]+"\t" + enhancer[2]+ "\t"
	interactedchromDetails = enhancerDetails + interactedChrom
	
	for item in interactedDetails:
		del item[10:14]
	#interactedDetails = "\n".join(["\t".join(item) for item in interactedDetails])
	#write_file(interactedDetails, "interacted.txt", "w")
	return interactedChrom

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





TF_EnhancerLoc = '/home/muhammad/PhD/Papers/Li(Nov16-2015)/Project(Enhancers)/Scripts/Efficient Search/GM 12878 Cell/GM 12878 Cell & P300 TF/Data/GSM935294_hg19_wgEncodeSydhTfbsGm12878P300IggmusPk.narrowPeak'
EP_interactionLoc = '/home/muhammad/PhD/Papers/Li(Nov16-2015)/Project(Enhancers)/Scripts/Efficient Search/GM 12878 Cell/GM 12878 Cell & P300 TF/EP_Interaction_for_all_chromosome_GM-12878_P300'


outputFile = 'GM12878_p300_EP statistics/GM12878_p300_EP statistics.txt'
make_directory(outputFile)

p300_enhancerData = read_file(TF_EnhancerLoc)
p300_enhancerData = [item.replace("\n","").split("\t") for item in p300_enhancerData]
for item in p300_enhancerData:
	del item[3:10]

p300_enhancerData.sort(key=operator.itemgetter(0))

for enhancer in p300_enhancerData[:2]:
	interactingPromoters = []	
	enhancerChromosome = enhancer[0]
	print "Enhancer details: "+ enhancerChromosome + " "+ enhancer[1]+"-" + enhancer[2]

	EP_MatchedFileLoc = getting_EP_InteractionFiles(EP_interactionLoc, enhancerChromosome)

	for ep_interactedLoc in EP_MatchedFileLoc:
		EP_interactedData = getData_FromFile(ep_interactedLoc, enhancerChromosome)
		EP_interactedData_sorted = sorted(EP_interactedData, cmp=contact)

		interactingPromoters.append(binary_search(enhancer, EP_interactedData_sorted))

	interactingPromoters = "----".join(interactingPromoters)
	enhancerDetails = enhancer[0] +"\t" + enhancer[1]+"\t" + enhancer[2]+ "\t"
	interactingPromoterDetails = enhancerDetails + interactingPromoters + "\n"

	write_file(interactingPromoterDetails, outputFile, "a")
	
print "--End--"