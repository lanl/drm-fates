import numpy as np
#import LLM_model as llm
import LLM_model_cpy as llm
#import LLM_model_class as llm
import time
import matk

def multiple_runs(p):
    # First we run the first run and save the results 
    nyears=200
    m = llm.LLM()     # assign p to the llm class
    m.dim=25
    m.instantiate(0)  # 1: reads input data from file, 0: generate inputs internally
    m.verbose=0       # 0: do not print out scores, 1: print scores on the screen
    m.fire_prob = p['par1'] # yearly fire probability, Ich: 0.35, Ord: 0.25
    m.mast_prob = p['par2'] # seed masting prob. for LP, based on Boyer's long-term data
    m.tree_mature_age=p['par3'] # Maturity age of the HW tree, for high fire prob it gets killed to often so start with small age
    # To generate fire and mast randomly 'readfireprobfromfile' and 'readmastprobfromfile' have to be turned off 
    m.readfireprobfromfile=0
    m.readmastprobfromfile=0
    m.run(nyears)
    m.save_pickle()
    
    # Then we use them as an input for sensitivity runs
    sc_rcw=[]
    sc_gt=[]
    sc_sq=[]
    n=5
    for i in range(n):
        m.instantiate(1) # 1: reads input data from file, 0: generate inputs internally
        m.run(nyears) # runs only 2 year
        sc_rcw.append((np.asarray(m.age_sc)+np.asarray(m.hw_sc)+
                       np.asarray(m.ageHW_sc)+np.asarray(m.hwHW_sc))*0.2)
        sc_gt.append(m.gt_sc)
        sc_sq.append(m.sq_sc)
        
    rcw=np.asarray(np.mean(sc_rcw,0))
    gt=np.asarray(np.mean(sc_gt,0))
    sq=np.asarray(np.mean(sc_sq,0))
    
    y=np.zeros(3)
    z=np.zeros(3)
    
    zz = np.polyfit(range(nyears), rcw, 1)
    z[0] = zz[0]
    zz = np.polyfit(range(nyears), gt, 1)
    z[1] = zz[0]
    zz = np.polyfit(range(nyears), sq, 1)
    z[2] = zz[0]
    
    y[0]=np.mean(rcw[-30:-1])
    y[1]=np.mean(gt[-30:-1])
    y[2]=np.mean(sq[-30:-1])
    print y
    
    ydict = dict([('obs'+str(i+1), v)  for i,v in enumerate(y)])
    ydict1 = dict([('slp'+str(i+1), v)  for i,v in enumerate(z)])
    print ydict

    with open("slopes.txt", "a") as myfile:
       myfile.write(str(ydict1))

    return ydict
 
with open('slopes.txt','wb') as f:
   f.write('')    

start_time = time.time()

p = matk.matk(model=multiple_runs)

p.add_par('par1',min=0.0,max=0.3) # fite prob
p.add_par('par2',min=0.0,max=0.15) #mast prob
p.add_par('par3',min=10,max=10,value=10) # hw=mature

s = p.parstudy(nvals=[4,4,1])
print s.samples.values


s.run( cpus=4, outfile='results.dat', logfile='log.dat',verbose=False)
print("--- %s seconds ---" % (time.time() - start_time))
