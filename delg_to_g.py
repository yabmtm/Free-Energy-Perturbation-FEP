#!/usr/bin/python

i = 0
mylist=[]

file="del_g.txt"

with open(file,'r') as f:
    while True:
	content = f.readline()
	if not content:
		break
	mylist.append(content.replace('\n',''))

for j in mylist:
    i+=float(j)
    print i
