import tqdm
from tqdm import tqdm
import pandas as pd
from pandas import DataFrame
cols = ['Number', 'Name' , 'Residue' ,'ResidueID', 'X', 'Y', 'Z', 'Time' ,'Area' ,' Thick']

data = open("DPPCDUPCCHOLaplvoro.txt", 'r')
lines = data.readlines()
N = 1124
outputData = []
indices = [] 
for index, line in tqdm(enumerate(lines)):
	data = []
	if line.split(' ')[0] == '#Framenumber':
		indices.append(index)


for index in range(1,len(indices)):
	data = []
	for ind in range(indices[index-1]+1, indices[index]):
		print(indices[index-1], indices[index])
		data.append(lines[ind])
	outputData.append(data)

	
for ind1, data in enumerate(outputData):
	replacementdata = []
	for line in data:
		line = [i for i in line.split(' ') if i is not '']
		replacementdata.append(line)
	outputData[ind1] = pd.DataFrame(replacementdata, columns = cols)
	
		
