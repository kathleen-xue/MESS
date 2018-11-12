import matplotlib.pyplot as plt
with open('./netBenefits.txt') as f:
	lines = f.readlines()
	y1 = [int(line.split()[0]) for line in lines]
	y2 = [int(line.split()[1]) for line in lines]


x1 = [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50] #DP

x2 = [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50] #greedy
#greedy

plt.plot(x1,y1, label = "MESS algorithm")
plt.plot(x2,y2, label = "simple greedy algorithm")

plt.xlabel('number of batteries')
plt.ylabel('total benefit')
plt.legend()
plt.title('Benefit of 2-50 batteries for 10 microgrids over 10 days')
plt.show()