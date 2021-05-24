import numpy as np
import hviscore as hvi
import LLM_model as llm


def multiple_runs(hw_mature,fprob,mprob):
    p = llm.LLM()     # assign p to the llm class
    p.verbose=0       # 0: do not print out scores, 1: print scores on the screen
    p.dim=15
    # Maturity age of the HW tree, for high fire prob it gets killed to often so start with small age
    p.mature_hw_age=hw_mature
    
    p.fire_prob = fprob # yearly fire probability, Ich: 0.35, Ord: 0.25
    p.mast_prob = mprob # seed masting prob. for LP, based on Boyer's long-term data
    
    # To generate fire and mast randomly 'readfireprobfromfile' and 'readmastprobfromfile' have to be turned off 
    p.readfireprobfromfile=0
    p.readmastprobfromfile=0
    sc_rcw=[]
    sc_gt=[]
    sc_sq=[]
    n=10
    nyears=200
    for i in range(n):
        p.instantiate(0) # 1: reads input data from file, 0: generate inputs internally
        p.run(nyears) # runs only 2 year
        sc_rcw.append((np.asarray(p.age_sc)+np.asarray(p.hw_sc)+
                       np.asarray(p.dbh_sc)+np.asarray(p.ageHW_sc)+np.asarray(p.hwHW_sc))*0.2)
        sc_gt.append(p.gt_sc)
        sc_sq.append(p.sq_sc)
        
    #sc_gt=np.reshape(sc_gt, (nyears, n))
    #sc_sq=np.reshape(sc_sq, (nyears, n))
    #sc_rcw=np.reshape(sc_rcw, (nyears, n))
    return sc_rcw,sc_gt,sc_sq
    

rcw=[]
gt=[]
sq=[]
f_prob=[0.0,0.1,0.2,0.3]
m_prob=[0.0,0.05,0.1,0.15]

for f in f_prob:
    for m in m_prob:
        print f,m
        [sc_rcw,sc_gt,sc_sq]=multiple_runs(30,f,m)
        rcw.append(np.mean(sc_rcw, axis=0))
        gt.append(np.mean(sc_gt, axis=0))
        sq.append(np.mean(sc_sq, axis=0))

#[sc_rcw,sc_gt,sc_sq]=multiple_runs(30,f_prob[0],m_prob[0])
#rcw.append(np.mean(sc_rcw, axis=0))
#gt.append(np.mean(sc_gt, axis=0))
#sq.append(np.mean(sc_sq, axis=0))

#rcw=np.transpose(rcw)
#sq=np.transpose(sq)
#gt=np.transpose(gt)
    
np.savetxt('results_rcw.txt', rcw, fmt='%1.3e')
np.savetxt('results_gt.txt', gt, fmt='%1.3e')
np.savetxt('results_sq.txt', sq, fmt='%1.3e')
