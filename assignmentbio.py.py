
# coding: utf-8

# In[118]:


from __future__ import division
import random
import re
import csv
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from array import *
import decimal
from scipy import stats



filepath = 'ATPbindingProtien.txt'

sequence_names = []
sequences = []

def Read_datset(filepath):
    print("Reading sequences ")
    
    f = open(filepath, 'r')
    counter = 0
    for i in f:
        i = i.strip()
        #print(i)
        if (i[0] == '>'):
            counter += 1
            sequence_names.append(i[1:].replace('\n', ''))
            sequences.append(str())
        else:
            sequences[counter - 1] = sequences[counter - 1] + i.replace('\n', '')
    f.close()
    #print(sequence_names,sequences)
    #return (sequence_names, sequences)
    
Read_datset(filepath)   



#print (len(sequences))
AA="GPAVILMCFYHKRQNEDST"
AAtotal_dict={'G':0,
                 'P':0,
                 'A':0,
                 'V':0,
                 'L':0,
                 'I':0,
                 'M':0,
                 'C':0,
                 'F':0,
                 'Y':0,
                 'W':0,
                 'H':0,
                 'K':0,
                 'R':0,
                 'Q':0,
                 'N':0,
                'E':0,
                 'D':0,
                 'S':0,
                'T':0,
                 }

AANItotal_dict={'G':0,
                 'P':0,
                 'A':0,
                 'V':0,
                 'L':0,
                 'I':0,
                 'M':0,
                 'C':0,
                 'F':0,
                 'Y':0,
                 'W':0,
                 'H':0,
                 'K':0,
                 'R':0,
                 'Q':0,
                 'N':0,
                'E':0,
                 'D':0,
                 'S':0,
                'T':0,
                 }







def AA_frequency(sequence_names,sequences):
    totalnfresidue=0

    for i in range (int(len(sequences)/2)):
        AA_dict={'G':0,
                 'P':0,
                 'A':0,
                 'V':0,
                 'L':0,
                 'I':0,
                 'M':0,
                 'C':0,
                 'F':0,
                 'Y':0,
                 'W':0,
                 'H':0,
                 'K':0,
                 'R':0,
                 'Q':0,
                 'N':0,
                'E':0,
                 'D':0,
                 'S':0,
                'T':0,
                 }
        AANI_dict={'G':0,
                 'P':0,
                 'A':0,
                 'V':0,
                 'L':0,
                 'I':0,
                 'M':0,
                 'C':0,
                 'F':0,
                 'Y':0,
                 'W':0,
                 'H':0,
                 'K':0,
                 'R':0,
                 'Q':0,
                 'N':0,
                'E':0,
                 'D':0,
                 'S':0,
                'T':0,
                 }
        print("////Individual sequence analysis /////")
        print (sequence_names[2*i])
        print(sequences[2*i])
        totalnfresidue+=len(sequences[2*i])
        print (sequences[2*i+1])
        
        j=0
        for j in range (len(sequences[2*i])):
        
            if (sequences[2*i+1][j]=='1'):
                letter=sequences[2*i][j]
                #print (letter)
            
                if letter in AA:
                    AA_dict[letter]+=1
                    AAtotal_dict[letter]+=1
                #print(AA_dict[letter])
            else :
                
                letter=sequences[2*i][j]
                if letter in AA:
                    AANI_dict[letter]+=1
                    AANItotal_dict[letter]+=1
                
                #AA_dict(sequences[i][j])+=1
                #print(AA_dict(sequences[i][j]))
        #j+=1
        
        
                #print(AA_dict[letter])
        #print (AA_dict)			
        with open('test.csv', 'w') as f:
            fields=['AA','Interacting frequency']
            writer=csv.DictWriter(f,fieldnames=fields)
            writer.writeheader()
            for key in AA_dict.keys():
                f.write("%s,%s\n"%(key,int(AA_dict[key])))  		 
        
       # print (AANI_dict)			
        
        
        with open('test1.csv', 'w') as g:
            fields=['AA','Non interacting frequency']
            writer=csv.DictWriter(g,fieldnames=fields)
            writer.writeheader()
            for key1 in AANI_dict.keys():
                #print(AANI_dict[key1])
                g.write("%s,%s\n"%(key1,int(AANI_dict[key1])))  		 

        
        
        i+=1
        dataset=pd.read_csv("test.csv", index_col=None)
        dataset1=pd.read_csv("test1.csv",index_col=None)
        
        #col1=dataset.iloc[:,[0]].values
        #col2=dataset.iloc[:,[1]].values
        #plt.plot(col1,col2)
        #print(col2)
        
        
        col1=dataset.iloc[:,[0]].values
        col2=dataset.iloc[:,[1]].values
        col3=dataset1.iloc[:,[1]].values
        #print(col3)
        #print(col2)
        
        
        print (dataset)
        
        for k in range(20):
             
            plt.bar(col1[k],int(col2[k]), 
               width = 0.7, color = 'yellow')
            #print(col1[k],int(col2[k]))
            plt.xlabel('AA') 
            # frequency label 
            plt.ylabel('Frequency') 
            # plot title 
            plt.title('Amino acid frequency plot') 
          
        if (i<=int(len(sequences)/2)):    
            plt.show()
            # x-axis label 
            plt.xlabel('AA') 
            # frequency label 
            plt.ylabel('Frequency') 
            # plot title 
            plt.title('Amino acid frequency plot') 
        
        
        print (dataset1)
        
        for k in range(20):
          
            plt.bar(col1[k],int(col3[k]), 
               width = 0.7, color = 'magenta')
            #print(col1[k],int(col3[k]))
           
        if (i<=int(len(sequences)/2)):    
            plt.show()
            # x-axis label 
            plt.xlabel('AA') 
            # frequency label 
            plt.ylabel('Frequency') 
            # plot title 
            plt.title('Amino acid frequency plot') 
        

    return totalnfresidue

        
        
totalnfresidue= AA_frequency(sequence_names,sequences)







with open('test2.csv', 'w') as g:
            fields=['AA','total frequency interacting']
            writer=csv.DictWriter(g,fieldnames=fields)
            writer.writeheader()
            for key1 in AAtotal_dict.keys():
                #print(AANI_dict[key1])
                g.write("%s,%s\n"%(key1,int(AAtotal_dict[key1])))  		 

with open('test3.csv', 'w') as g:
            fields=['AA','total frequency non interacting']
            writer=csv.DictWriter(g,fieldnames=fields)
            writer.writeheader()
            for key1 in AANItotal_dict.keys():
                #print(AANI_dict[key1])
                g.write("%s,%s\n"%(key1,int(AANItotal_dict[key1])))  		 

                
                
dataset2=pd.read_csv("test2.csv",index_col=None)
dataset3=pd.read_csv("test3.csv",index_col=None)
col1=dataset2.iloc[:,[0]].values
col2=dataset2.iloc[:,[1]].values
col3=dataset3.iloc[:,[1]].values
#col4=dataset3.iloc[:,[1]].values
p=[]
AAcomposition=[]
for k in range(20):
    if (((int(col2[k]))+(int(col3[k])))!=0):
        #col4[k]=k
        num=int(col2[k])*100
        denom=int(col2[k])+int(col3[k])
        p.append(round(num/denom,3))
    else :
         p.append(0.0)
            
            
total=sum(col2)+sum(col3)
            
#print (total)  
totalint=sum(col2)
#print (totalint)
#print(totalnfresidue)

for k in range(20):
        
        num=int(col2[k])*100
        AAcomposition .append(round(num/totalint[0],3))


 



        
        

        
        

print("\\\\\Total Frequency Plot\\\\")


print (dataset2)
        
for k in range(20):
           
    plt.bar(col1[k],int(col2[k]), 
        width = 0.7, color = 'blue')
            #print(col1[k],int(col2[k]))
plt.show()        
plt.xlabel('AA') 
            # frequency label 
plt.ylabel('Frequency') 
            # plot title 
plt.title('Amino acid frequency plot') 
          

            
print (dataset3)
        
for k in range(20):
          
    plt.bar(col1[k],int(col3[k]), 
        width = 0.7, color = 'red')
            #print(col1[k],int(col3[k]))
           
plt.show()


print ("\\\\Analysis Plot\\\\")

print("Percentage of amino acid composition for the all interactive residues")
for k in range(20):
    print (col1[k],"=",AAcomposition[k])      
    plt.bar(col1[k],AAcomposition[k], 
        width = 0.7, color = 'orange')
            #print(col1[k],int(col3[k]))
    plt.xlabel('AA') 
            # frequency label 
    plt.ylabel('AA composition percentage') 
            # plot title 
    plt.title('Amino acid vs composition plot')  



plt.show()
print("Propensity  score for the interacting residues")
for k in range(20):
    print (col1[k],"=",p[k])      
    plt.bar(col1[k],p[k], 
        width = 0.7, color = 'Green')
            #print(col1[k],int(col3[k]))
    plt.xlabel('AA') 
            # frequency label 
    plt.ylabel('Propensity score') 
            # plot title 
    plt.title('Amino acid vs  Propensity score plot')        

          









