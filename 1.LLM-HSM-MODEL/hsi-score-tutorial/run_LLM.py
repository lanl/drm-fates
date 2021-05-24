#!/usr/bin/env python
import numpy as np
import LLM_model as llm
import hviscore as hvi


p = llm.LLM()     # assign p to the llm class
# Maturity age of the HW tree, for high fire prob it gets killed to often so start with small age

p.dim=15
p.verbose=1       # 0: do not print out scores, 1: print scores on the screen
p.mature_hw_age= 30  
p.fire_prob = 0.10 # yearly fire probability, Ich: 0.35, Ord: 0.25
p.mast_prob = 0.15 # seed masting prob. for LP, based on Boyer's long-term data
    
# To generate fire and mast randomly 'readfireprobfromfile' and 'readmastprobfromfile' have to be turned off 
p.readfireprobfromfile=0
p.readmastprobfromfile=0

sc_rcw=[]
sc_gt=[]
sc_sq=[]

nyears=2
   
p.instantiate(1) # 1: reads input data from file, 0: generate inputs internally
p.run(nyears) # runs only 2 year
sc_rcw.append((np.asarray(p.age_sc)+np.asarray(p.hw_sc)+np.asarray(p.ageHW_sc)+np.asarray(p.hwHW_sc))*0.25)
sc_gt.append(p.gt_sc)
sc_sq.append(p.sq_sc)
