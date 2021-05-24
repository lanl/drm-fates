#!/usr/bin/env python
import run_LLM as frun
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

xx,yy,zz,xx1,yy1= frun.call_LLM()
plt.plot(xx,'r')
plt.plot(yy,'b')
plt.plot(zz,'g--')
plt.plot(xx1,'k')
plt.plot(yy1,'k--')

'''
fig, ax1 = plt.subplots()

ax1.plot(xx, 'b-',linewidth=2)
ax1.set_xlabel('time (yr)')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('score age', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(yy, 'r*')
ax2.set_ylabel('num', color='r')
ax2.tick_params('y', colors='r')
'''
#fig.tight_layout()


#xx,x,y,z=frun.call_LLM()
#np.savetxt('test_LLM.out', [x,y,z], fmt='%1.4e') 

#x=frun.call_LLM()

#fig, ax = plt.subplots(figsize=(15, 12))
#im = ax.imshow(xx, interpolation = 'nearest')#, vmin=0, vmax=200) 
#fig.colorbar(im, ax=ax)
#plt.title('LLP tree age')
## Loop over data dimensions and create text annotations.
#for i in range(len(xx[:,0])):
#    for j in range(len(xx[0,:])):
#        text = ax.text(j, i, xx[i, j],
#                       ha="center", va="center", color="w")
                       
#Y = sp.spatial.distance.pdist(x, 'euclidean')
#print Y

plt.show()

