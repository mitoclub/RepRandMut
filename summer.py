import time
import timeit
import os




file = open("subrep.csv",'r')



reps3 = open("subrep_sum.csv", "w")


name = file.readline().strip()



reps3.write("pos;sub;-;=;+\n")



mc = []
for i in range(16600):
	mc.append([[],[],[],[]])

#mc = [[[],[],[],[]]]*16569

while True:


	#print(mc[:10])

	# 6;6;A>C;-

	read = file.readline().strip()

	if read == "":
		break

	row = read.split(";")


	print(row)


	if row[2][2] == "A":
		if mc[int(row[0])][0] == []:
			mc[int(row[0])][0] = [row[2],0,0,0]
		if row[3] == "-":
			mc[int(row[0])][0][1]+=1
		elif row[3] == "=":
			mc[int(row[0])][0][2]+=1		
		elif row[3] == "+":
			mc[int(row[0])][0][3]+=1	
		
	elif row[2][2] == "T":
		if mc[int(row[0])][1] == []:
			mc[int(row[0])][1] = [row[2],0,0,0]
		if row[3] == "-":
			mc[int(row[0])][1][1]+=1
		elif row[3] == "=":
			mc[int(row[0])][1][2]+=1		
		elif row[3] == "+":
			mc[int(row[0])][1][3]+=1

	elif row[2][2] == "G":
		if mc[int(row[0])][2] == []:
			mc[int(row[0])][2] = [row[2],0,0,0]
		if row[3] == "-":
			mc[int(row[0])][2][1]+=1
		elif row[3] == "=":
			mc[int(row[0])][2][2]+=1		
		elif row[3] == "+":
			mc[int(row[0])][2][3]+=1

	elif row[2][2] == "C":
		if mc[int(row[0])][3] == []:
			mc[int(row[0])][3] = [row[2],0,0,0]
		if row[3] == "-":
			mc[int(row[0])][3][1]+=1
		elif row[3] == "=":
			mc[int(row[0])][3][2]+=1		
		elif row[3] == "+":
			mc[int(row[0])][3][3]+=1




for i in range(len(mc)): 
	if mc[i] != [[],[],[],[]]:
		for q in mc[i]:	
			if q != []:
				reps3.write("%s;%s;%s;%s;%s\n" % (i+1,q[0],q[1],q[2],q[3]))