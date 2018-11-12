import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

mess = []

with open('./netBenefitsMESS.txt') as f:
	lines = f.readlines()
	for i in range(len(lines)):
		mess.append([int(line) for line in lines[i].split(" ")[:-1]])

greed = []

with open('./netBenefitsGREEDY.txt') as f:
	lines = f.readlines()
	for i in range(len(lines)):
		greed.append([int(line) for line in lines[i].split(" ")[:-1]])

x1 = []
x2 = []

for i in range(2,52,2):
	x1.append(i) #DP
	x2.append(i) #greedy


avgM = []
currAM = 0
avgG = []
currAG = 0

if len(mess) == 1:
	plt.plot(x1, mess[0], color = 'xkcd:raspberry')
	plt.plot(x2, greed[0], color = 'xkcd:cobalt blue')

elif len(mess) > 1:
	for i in range(0,len(mess)):
		plt.plot(x1, mess[i], color = 'xkcd:raspberry', alpha = 0.3)
		plt.plot(x2, greed[i], color = 'xkcd:cobalt blue', alpha = 0.3)

	for j in range(0,len(mess[0])):
		for i in range(0,len(mess)):
			currAM += mess[i][j]
			currAG += greed[i][j]
		currAM /= len(mess)
		currAG /= len(mess)
		avgM.append(currAM)
		avgG.append(currAG)
	
	plt.plot(x1, avgM, color = 'xkcd:wine red', label = 'MESS algorithm')
	plt.plot(x2, avgG, color = 'xkcd:royal blue', label = 'simple greedy algorithm')

plt.xlabel('number of batteries')
plt.ylabel('total net cost (cost + benefit)')
plt.legend()
plt.title('Benefit of 2-50 batteries for 10 microgrids over 10 days')
plt.show()
