#/usr/bin/python

#     Message Passing Algoritm Implementation for Max-Sum:
#     Written by: Ramezan Paravi

import numpy as np 
import sys
import itertools
import randompermutation
from numpy import prod
from collections import deque
import math
from math import log
from operator import sub
from operator import add

import sage
from sage.all import *
# class of factor graph
class FactorGraph:
    def __init__(self):
        self.FacNodeList=[] #factor nodes list
        self.VarNodeList=[] #variable node list

    def AddVarNode(self, label):
    	newNode=VarNode(label)
        self.VarNodeList.append(newNode)
        return newNode

    def AddFacNode(self, label):
    	newNode=FacNode(label)
        self.FacNodeList.append(newNode)
        return newNode
    def FindVarNode(self,id):
        for w in self.VarNodeList:
            if w.NodeLabel==id:
                return w
    def FindFacNode(self,id):
        for w in self.FacNodeList:
            if w.NodeLabel==id:
                return w
    def FindFacSoketID(self,soket):
        for w in self.FacNodeList:
            if soket in w.socket_id:
                return w
                
# class of variable node
class VarNode:
    def __init__(self, label ):

        self.NodeLabel=label
        self.NeighList=[]
        self.socket_id=[]
        self.message_in=[]
        self.message_out=[]

#***********************************************************************************************************************************************
# class of factor node
class FacNode:
    def __init__(self, label):
        self.NodeLabel=label
        self.NeighList=[]
        self.socket_id=[]
        self.message_in=[]
        self.message_out=[]



#For Definitions of polynomials used in simulation such as Lambda and P, see T. Richardson and R. Urbanke, Modern Coding Theory. Cambridge University Press, 2008.
# Length of codeword
N=106; 

#Variable Node Polynomial
Lambda=[(64,2),(24,3),(5,4),(13,8)] # Lambda(x)=64x^2 + 24*x^3+5*x^4+13*x^8
P=[(54,6)] # P(x)=54*x^6

num_of_edges=0
num_of_checks=54
for i in range(0,len(Lambda)):
	num_of_edges=num_of_edges+Lambda[i][0]*Lambda[i][1]

x=randompermutation.randomperm(num_of_edges)

#Form Bipartite Factor Graph
T=FactorGraph()
for i in range(1,N+1):
	T.AddVarNode(i)
for i in range(1,num_of_checks+1):
	T.AddFacNode(i)

counter1=0
counter2=0
for i in range(0,len(Lambda)):
	s=Lambda[i][0]
	for j in range(1,s+1):
		counter1=counter1+1
		w=T.FindVarNode(counter1)
		for k in range(1,Lambda[i][1]+1):
			counter2=counter2+1
			w.socket_id.append(counter2)

counter1=0
counter2=0
for i in range(0,len(P)):
	s=P[i][0]
	for j in range(1,s+1):
		counter1=counter1+1
		w=T.FindFacNode(counter1)
		for k in range(1,P[i][1]+1):
			counter2=counter2+1
			w.socket_id.append(counter2)

for w in T.VarNodeList:
	for z in w.socket_id:
		q=T.FindFacSoketID(x[z-1])
		w.NeighList.append(q.NodeLabel)
		q.NeighList.append(w.NodeLabel)

for w in T.VarNodeList:
	w.message_in=[0]*len(w.NeighList)
	w.message_out=[0]*len(w.NeighList)
for w in T.FacNodeList:
	w.message_in=[0]*len(w.NeighList)
	w.message_out=[0]*len(w.NeighList)

#===========================================================================================================	
# Channel Simulation: Using Sage to create valid codeword

M=MatrixSpace(GF(2),num_of_checks,N,sparse=True)
H=M.matrix()
for w in T.FacNodeList:
	i=T.FindFacNode(w.NodeLabel)
	i=i.NodeLabel-1
	for q in w.NeighList:
		j=T.FindVarNode(q)
		j=j.NodeLabel-1
		H[i,j]=1

C=codes.LinearCodeFromCheckMatrix(H)
G=C.gen_mat_systematic()
I=MatrixSpace(GF(2),1,N-num_of_checks)
m=I.random_element()
#print m[0]
codeword=m*G
test=H*codeword.transpose()
print "Test result if the sequence is actually a codeword:"
print test.is_zero()

X=[]
for w in codeword[0]:
	if w==1:
		X.append(1)
	else:
		X.append(0)
print "Codeword is:"
print X

#Erasure Probability
epsilon=0.05
Y=[0]*len(X)
for i in range(0,len(X)):
	 u=np.random.uniform(0,1,1)
	 if u<epsilon:
	 	Y[i]='?'
	 else:
	 	Y[i]=X[i]
print "Recieved sequence is:"
print Y
#===========================================================================================================	

#Message Passing Algorithm:
iteration=3000


#Initialization:
temp=0
for w in T.VarNodeList:
	w.message_out=[Y[temp]]*len(w.NeighList)
	for v in w.NeighList:
		v=T.FindFacNode(v)
		index=v.NeighList.index(w.NodeLabel)
		v.message_in[index]=Y[temp]
	temp=temp+1


for iter in range(1,iteration+1):
	#print iter
	#Message sent out from check nodes:
	for w in T.FacNodeList:
		for s in w.NeighList:
			index=w.NeighList.index(s)
			r=T.FindVarNode(s)
			index3=r.NeighList.index(w.NodeLabel)
			my_indexes=range(0,len(w.NeighList))
			my_indexes.remove(index)

			indicator=0
			for l in my_indexes:
				q=T.FindVarNode(w.NeighList[l])
				index2=q.NeighList.index(w.NodeLabel)
				if q.message_out[index2]=='?':
					indicator=1
					w.message_out[index]='?'
					r.message_in[index3]='?'
			if indicator==0:
				temp3=0
				for l in my_indexes:
					q=T.FindVarNode(w.NeighList[l])
					index2=q.NeighList.index(w.NodeLabel)
					temp3=temp3+q.message_out[index2]
				if temp3%2==0:
					w.message_out[index]=0
				else:
					w.message_out[index]=1

	#Message sent out from variable nodes:
	for w in T.VarNodeList:
		if Y[w.NodeLabel-1]=='?':
			for s in w.NeighList:
				index=w.NeighList.index(s)
				r=T.FindFacNode(s)
				index3=r.NeighList.index(w.NodeLabel)
				my_indexes=range(0,len(w.NeighList))
				my_indexes.remove(index)

				indicator=0
				for l in my_indexes:
					q=T.FindFacNode(w.NeighList[l])
					index2=q.NeighList.index(w.NodeLabel)
					if q.message_out[index2]!='?':
						indicator=1
						zz=q.message_out[index2]
						#w.message_out[index]='?'
						#r.message_in[index3]='?'
				if indicator==0:
					w.message_out[index]='?'
					r.message_in[index3]='?'
				else:
					w.message_out[index]=zz
					r.message_in[index3]=zz
					Y[w.NodeLabel-1]=zz



print Y
if Y==X :
	print "Decode Correctly."
else:
	print "Decode Incorrectly. Increase iterations numbers."