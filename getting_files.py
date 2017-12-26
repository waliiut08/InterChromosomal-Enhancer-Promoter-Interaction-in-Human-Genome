import os
from common_functions import write_file, find_between

EP_interactionLoc = '/home/muhammad/PhD/Papers/Li(Nov16-2015)/Project(Enhancers)/Scripts/Efficient Search/GM 12878 Cell/GM 12878 Cell & P300 TF/EP_Interaction_for_all_chromosome_GM-12878_P300'


matchedFiles = []
EnhancerChromosome = 'chr1'

## getting the EP file location that matches the Enhancer chromosome------------------------------------------------------------------------------------------------

def getting_EP_files(EP_interactionLoc, EnhancerChromosome):
	i =0

	for dirpath,dirname,filename in os.walk(EP_interactionLoc):	
		if i>0:
			interactingChromosomes = find_between(filename[0], "Enhancer-Promoter Interactions and normalized vector between " , ".txt")
			interactingChromosomes = interactingChromosomes.split("_")
			
			if interactingChromosomes[0]==EnhancerChromosome or interactingChromosomes[1]==EnhancerChromosome:
				matchedFiles.append(dirpath +"/" + filename[0])

		i+=1



