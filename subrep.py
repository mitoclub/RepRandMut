#coding=utf8


# 20230207: Fix it for work




# Importing of libraries
import time
import timeit
import os
import sqlite3
import sys

try:
	from StringIO import StringIO ## for Python 2
except ImportError:
	from io import StringIO ## for Python 3


"""
Structure of working tables:
repeats (
id INTEGER PRIMARY KEY,  
spece_id INTEGER,			  - id of species in table species
first_start INTEGER,		   - position of first sequence in genome
second_start INTEGER,		  - position of second sequence in genome
length INTEGER,				- length of repeat
id_type INTEGER,			   - type of repeat
first_seq VARCHAR(1000),	   - first sequence
second_seq VARCHAR(1000),	  - second sequence
alt_second_seq VARCHAR(1000),  - second sequence transformed in accordance to the type of repeat
errs (INTEGER)				 - number of errors
)
species (
id INTEGER PRIMARY KEY,		- id of species
spece TEXT					 - species name
)		  
types (
id INTEGER PRIMARY KEY,		- id of type of repeat
type TEXT					  - type of repeat
)
"""

front=[0,99999] # !!!!!!!!!!!!!!!!!!!!!!!!!

#AncestorsOfHominida.muscle.sequester.fasta
#Mitochondria_GreatApes.muscle.sequester.fasta
#SupraHominida.muscle.sequester.fasta
pgen_path = 'Homo_sapiens.fa'

folder = "REPS"


dirpath = os.getcwd() + "/" + folder
# if __name__ == "__main__":

#	 if len (sys.argv) > 2:
#		 pgen_path = sys.argv[1]
#		 dirpath = os.getcwd() + "/" + sys.argv[2]
#	 elif len (sys.argv) > 1:
#		 pgen_path = sys.argv[1]
pgen = open(pgen_path,'r')
nucAlphabet = {'a', 't', 'g', 'c', 'x'}

gcount=sum(1 for line in pgen)/2
print("We have %s genomes for analysis" % gcount)

pgen = open(pgen_path,'r')


try:
	os.mkdir(dirpath)
except OSError:
	print ("FOUND directory is exist")
else:
	print ("FOUND directory created")


# creating of empty dot-plot (matrix)
def zeros(shape):
	retval = []
	for x in range(shape[0]):
		retval.append([])
		for y in range(shape[1]):
			retval[-1].append(0)
	return retval

# Slicing sequenses from genome of length 10 with a shift of 1
def extractWordsTen(s, alphabet):
	arr = []
	n=10
	for i in range(len(s)-n+1):
		arr.append(''.join(s[i:i+n]))
	return arr

# Counter of searching time
class Profiler(object):
	def __enter__(self):
		self._startTime = time.time()	  
	def __exit__(self, type, value, traceback):
		print("Elapsed time: {:.3f} sec".format(time.time() - self._startTime))






# Computation of Levenshtein distance between two sequence
def lev(s, t):
    if s == t: return 0
    elif len(s) == 0: return len(t)
    elif len(t) == 0: return len(s)
    v0 = [None] * (len(t) + 1)
    v1 = [None] * (len(t) + 1)
    for i in range(len(v0)):
        v0[i] = i
    for i in range(len(s)):
        v1[0] = i + 1
        for j in range(len(t)):
            cost = 0 if s[i] == t[j] else 1
            v1[j + 1] = min(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost)
        for j in range(len(v0)):
            v0[j] = v1[j]
    return v1[len(t)]










# Matches checking
# def  match_check (words,invert,score,oi,oj,i,j,e,reptype,origlen):
# 	lengthj=len(score[oi])
# 	lengthi=len(score)
# 	match=1 # matches
# 	dmatch=e # mismatches
# 	# ten cells passage
# 	for k in range(1, 10-dmatch): 
# 		ci=i+k
# 		cj=j+k
# 		if ((ci>(lengthi-1)) or (cj>(lengthj-1))): # Output to the end of the genome	
# 			return 0
# 		if score[ci][cj]==10: # Match of the nucleotides marked in the cell
# 			match=match+1  
# 			if (match>8): # If there is 9-10 matches record a repeat
# 				if invert:
# 					point=lengthj-oj-11 # Position for inverted sequence
# 				else:
# 					point=oj
# 				if point>origlen-1:
# 					point=point-origlen   
# 				if oi>origlen-1:
# 					oi=oi-origlen
# 				q = "INSERT INTO repeats (id,first_start, second_start, length, id_type, first_seq, second_seq) VALUES(NULL,?,?,?,?,?,?)"
# 				cur.execute(q,(oi,point,10,reptype,words[oi],words[point]))
# 				return 1
# 		else:				 
# 			dmatch=dmatch+1
# 			if (dmatch>1) :		  
# 				return 0
# 	if (match>8):
# 		if invert:
# 			point=lengthj-oj-11
# 		else:
# 			point=oj
# 		if point>origlen-1:
# 			point=point-origlen   
# 		if oi>origlen-1:
# 			oi=oi-origlen
# 		q = "INSERT INTO repeats (first_start, second_start, length, id_type, first_seq, second_seq) VALUES(?,?,?,?,?,?)"
# 		cur.execute(q,(oi,point,10,reptype,words[oi],words[point]))
# 		return 1
# 	else: 	
# 		return 0




path = 'subrep.csv'
reps = open(path, "w")



def substits(nuc):
	arr = []
	pot = ['A','T','G','C']
	for p in pot:
		if p is not nuc:
			arr.append(p)
	return arr




# The procedure for working with the matrix
def master (seq1, seq2, invert,reptype,origlen):
	seq2=seq2.replace("X","N")
	m, n = len(seq1), len(seq2)
	words=extractWordsTen(seq1, nucAlphabet)
	

	for i in range(len(words)-1):
		print(i)
		for j in range(len(words)-10):

			if j<i-10 or j>i+9:

				errs=lev(words[i],words[j])

				if errs < 3: #####

					orierrs = errs

					if errs < 2: #####
						q = "INSERT INTO repeats (first_start, second_start, length, id_type, first_seq, second_seq, errs, subt) VALUES(?,?,?,?,?,?,?,?)"
						cur.execute(q,(i,j,10,1,words[i],words[j],errs,"orig"))


			
					for x in range(10):
						#print(words[i])
						nucs = substits(words[i][x])
						for n in nucs:	
							neword = words[i][:x] + n + words[i][x+1:]
							errs = lev(neword,words[j])

							if errs < 2: #####
								q = "INSERT INTO repeats (first_start, second_start, length, id_type, first_seq, second_seq, errs, subt) VALUES(?,?,?,?,?,?,?,?)"
								cur.execute(q,(i,j,10,1,words[i],words[j],errs,"alt:%s%s>%s" % (x,words[i][x],n) ))

							gloch = i + x
						
							if errs < orierrs:
								reps.write("%s;%s;%s>%s;+\n" % (gloch, x, words[i][x], n))
							elif errs == orierrs:
								reps.write("%s;%s;%s>%s;=\n" % (gloch, x, words[i][x], n))
							elif errs > orierrs:
								reps.write("%s;%s;%s>%s;-\n" % (gloch, x, words[i][x], n))


		con.commit()



	# score = zeros((m+1, n+1)) 


	 
	# print("matching")
	# for i in range(0, m):	   
	# 	print("%s\r" % i),
	# 	for j in range(i+10, n):
	# 		if seq1[i] == seq2[j]:
	# 			score[i][j]=10
	# print("checking")
	# lengthj=lengthi=len(score)
	# for i in range(0, m-10  +1):
	# 	print("%s\r" % i),
	# 	for j in range(i+10, n-10): 
	# 		match_check(words,invert,score,i,j,i,j,0,reptype,origlen)




# Preprocedure of genome preparation for finding complimentary repeats
# def mastercomp (seq1, seq2,origlen):
#	 reptype=2	  
#	 print("\nCompliment:")
#	 arr1=[]
#	 s2=""
#	 for char in seq1:
#		 if char=='A':
#			 arr1.append('T') 
#		 elif char=='T':
#			 arr1.append('A')
#		 elif char=='G':
#			 arr1.append('C')
#		 elif char=='C':
#			 arr1.append('G')
#		 else:
#			 arr1.append('X')
#	 s2=(''.join(arr1))
#	 master (seq1, s2,False,reptype,origlen)

# Preprocedure of genome preparation for finding mirror repeats
# def masterpoli (seq1, seq2,origlen):
#	 reptype=3	  
#	 print("\nMirror:")
#	 arr2=[]
#	 s3=""
#	 for char in seq1:
#		 arr2.append(char)
#	 s3=(''.join(reversed(arr2)))
#	 master (seq1, s3,True,reptype,origlen)

# Preprocedure of genome preparation for finding inverted repeats
# def masterinvert (seq1, seq2,origlen):
#	 reptype=4
#	 print("\nInvert:")	
#	 arr1=[]
#	 s2=""
#	 for char in seq1:  
#		 if char=='A':
#			 arr1.append('T')
#		 elif char=='T':
#			 arr1.append('A')
#		 elif char=='G':
#			 arr1.append('C')
#		 elif char=='C':
#			 arr1.append('G')
#		 else:
#			 arr1.append('X')
#	 arr3=reversed(arr1)
#	 s4=""
#	 s4=(''.join(arr3))
#	 master (seq1, s4,True,reptype,origlen)



# Alternate start searching for different types of repeats in each genome

pr=1
while True:
	if pr>front[0]:
		break
	name = pgen.readline().strip()[1:]
	genome = pgen.readline().strip()

	pr+=1


while True:
	# if pr>front[1]:
	#	 break



	name = pgen.readline().strip()[1:]

	if name == "":
		break


	genome = pgen.readline().strip()

	print("GENOME: " + name)
	length=len(genome)
	genome=genome[:-1]+genome[:10] # genome expansion 
	print('%s/%s_%s.sqlite' % (dirpath,pr,name))
	con = sqlite3.connect('%s/%s_%s.sqlite' % (dirpath,pr,name))
	cur = con.cursor()
	cur.execute('CREATE TABLE species (id INTEGER PRIMARY KEY, name VARCHAR(100))')
	con.commit()
	cur.execute('CREATE TABLE repeats (id INTEGER PRIMARY KEY, spece_id INTEGER, first_start INTEGER, second_start INTEGER, length INTEGER, id_type INTEGER, first_seq VARCHAR(50), second_seq VARCHAR(50), errs INTEGER, subt VARCHAR(50))')
	con.commit()
	reptype=1
	q = "INSERT INTO species (id, name) VALUES(NULL, \""+name+"\")"
	cur.execute(q)
	
	with Profiler() as p:   
		with Profiler() as a:
			print("\nSuper direct  for "+name+":")
			master(genome,genome,False,1,length)
			con.commit()
		# with Profiler() as b:
		#	 print("\nSuper compliment for "+name+":")
		#	 mastercomp(genome,genome,length)
		#	 con.commit()		  
		# with Profiler() as c:
		#	 print("\nSuper mirror for "+name+":")
		#	 masterpoli(genome,genome,length)
		#	 con.commit()
		# with Profiler() as d:
		#	 print("\nSuper invert for "+name+":")
		#	 masterinvert(genome,genome,length)
		#	 con.commit()
		# print("\nSuper all for "+name+":")
	con.commit()
	
	pr=pr+1
	if pr>gcount:
		break

