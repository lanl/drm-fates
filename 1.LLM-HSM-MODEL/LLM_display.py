
# Display Class for LLM_final_2010.py
# Longleaf pine and hardwood model, the "LLM", Loudermilk et al. 2011, Ecological Modeling paper.

from pylab import *

class display:

    def __init__(self):
        
        self.old_ht = None # not necessary to run with LLM 
        self.old_age = None
        self.litter = None
        self.litterHW = None
        self.litterWG = None
        self.tot_biomass = None
        self.seeds = None
        self.germ_seeds = None
        self.new_ht = None
        self.sum_ht = None
        self.old_LPcount = None
        self.old_HWcount = None
        self.fireflag = None
        self.mortflag = None
        self.old_htHW = None
        self.new_LPcount = None
        self.new_HWcount = None
        self.dead_htLP = None
        self.dead_htHW = None
        self.dead_ctHW = None
        self.dead_ctLP = None
        self.new_htHW = None
    

    def display_stand(self, yr = 5):

        '''graphing spatial variables
        self.displayed is the chosen variable in the init() method'''

        '''        self.g_types = {'ht':self.old_ht,
                        'age':self.old_age,
                        'seeds':self.seeds,
                        'sum_ht':self.sum_ht,
                        'litter':self.litter,
                        'f_fire':self.fireflag}'''
       
        try:
            d = self.displayed
        except:
            pass

        if d == 'ht':
            
            vm = 0.0
            vx = 40.0
            t = 'Longleaf Height (m)'
            imshow(self.old_ht, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)
            
        elif d == 'htHW':
            
            vm = 0.0
            vx = 20.0
            t = 'Hardwood Height (m)'
            imshow(self.old_htHW, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)            

        elif d == 'age':
            
            vm = 0.0
            vx = self.fintim + 75.0
            t = 'Longleaf Age (yr)'
            imshow(self.old_age, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)

        elif d == 'ageHW':
            
            vm = 0.0
            vx = self.fintim + 75.0
            t = 'Hardwood Age (yr)'
            imshow(self.old_ageHW, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)
            
        elif d == 'litter':
            
            vm = 0.0
            vx = 25.0 # kg
            t = 'Pine Needle Litter - kg/cell'
            imshow(self.litter, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)

        elif d == 'litterHW':
            
            vm = 0.0
            vx = 5.0 # kg
            t = 'HW Leaf Litter - kg/cell'
            imshow(self.litterHW, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)

        elif d == 'WG':
            
            vm = 0.0
            vx = self.max_wtWG # kg
            t = 'Wiregrass - kg/cell'
            imshow(self.litterWG, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)

        elif d == 'biomass':
            
            vm = 0.0
            vx = 15.0 # kg
            t = 'Total Biomass - kg/cell'
            imshow(self.tot_biomass, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)
            
        elif d == 'seeds':
            
            vm = 0.0
            vx = 200.0
            t = 'Seeds'
            imshow(self.seeds, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)

        elif d == 'germ_seeds':
            
            vm = 0.0
            vx = 10.0
            t = 'Germinated Seedlings'
            imshow(self.seeds, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)

        elif d == 'germ': # delete?
            
            vm = 0.0
            vx = 15.0
            t = 'Germ Seedlings and Existing Trees'
            imshow(self.old_ht, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)

        elif d == 'sum_ht':
            
            vm = 0.0
            vx = 2500.0
            t = 'Height Sum (m)'
            imshow(self.sum_htLP1, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)


        elif d == 'LPcount':
            
            vm = 0.0
            vx = 10.0
            t = 'LP Count'
            imshow(self.old_LPcount, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)

        elif d == 'HWcount':
            
            vm = 0.0
            vx = 10.0
            t = 'HW Count'
            imshow(self.old_HWcount, vmin = vm, vmax = vx,
               interpolation = 'nearest')
            title(t)
            
        elif d == 'f_fire':
            
            gray()
            t = 'Fire Flag'
            imshow(self.fireflag, interpolation = 'nearest')
            title(t)

        elif d == 'f_mort':
            
            gray()
            t = 'Mortality Flag LP'
            imshow(self.mortflag, interpolation = 'nearest')
            title(t)



        elif d == 'tree_histLP': # LP tree distribution
            
            self.ycount = self.ycount + 1
            
            if self.ycount < 5:
                return
            
            th = []
            for row in range(len(self.new_ht)):
                for col in range(len(self.new_ht)):
                    if self.new_ht[row][col] > 2.5:
                        for i in range(self.new_LPcount[row][col]): # incorp. no. trees w/in cell
                            th.append(self.new_ht[row][col])                     
            clf()
            try:
                h, bins, patches = hist(th, 25)
                setp (patches, 'facecolor', 'g', 'alpha', 0.75)
                axis([0.0,30.0,0.0,150.0])
                title('LP Tree Height Distribution (> 2.5 m)    Year = ' + str(yr + 1))
                xlabel('m')
                ylabel('no. LLP trees')
                
            except:
                pass
            self.ycount = 0

        elif d == 'dtree_histLP': # dead LP tree distribution
            
            self.ycount = self.ycount + 1
            
            if self.ycount < 5:
                return
            
            th = []
            for row in range(len(self.new_ht)):
                for col in range(len(self.new_ht)):
                    if self.dead_htLP[row][col] > 2.5:
                        for i in range(self.dead_ctLP[row][col]): # incorp. no. trees w/in cell
                            th.append(self.dead_htLP[row][col])                     
            clf()
            try:
                h, bins, patches = hist(th, 25)
                setp (patches, 'facecolor', 'g', 'alpha', 0.75)
                axis([0.0,30.0,0.0,50.0])
                title('LP Dead Tree Height Distribution (> 2.5 m)    Year = ' + str(yr + 1))
                xlabel('m')
                ylabel('no. LLP trees')
                
            except:
                pass
            self.ycount = 0
            
        elif d == 'tree_histHW': # HW tree distribution
            
            self.ycount = self.ycount + 1
            
            if self.ycount < 5:
                return
            
            th = []
            for row in range(len(self.new_htHW)):
                for col in range(len(self.new_htHW)):
                    if self.new_htHW[row][col] > 2.0:
                        for i in range(self.new_HWcount[row][col]): # incorp. no. trees w/in cell
                            th.append(self.new_htHW[row][col])                     
            clf()
            try:
                h, bins, patches = hist(th, 25)
                setp (patches, 'facecolor', 'g', 'alpha', 0.75)
                axis([0.0,30.0,0.0,150.0])
                title('Tree Height Distribution (> 2.0 m)    Year = ' + str(yr + 1))
                xlabel('m')
                ylabel('no. HW trees')
                
            except:
                pass
            self.ycount = 0

        elif d == 'dtree_histHW': # dead LP tree distribution
            
            self.ycount = self.ycount + 1
            
            if self.ycount < 5:
                return
            
            th = []
            for row in range(len(self.new_htHW)):
                for col in range(len(self.new_htHW)):
                    if self.dead_htHW[row][col] > 2.0:
                        for i in range(self.dead_ctHW[row][col]): # incorp. no. trees w/in cell
                            th.append(self.dead_htHW[row][col])                     
            clf()
            try:
                h, bins, patches = hist(th, 25)
                setp (patches, 'facecolor', 'g', 'alpha', 0.75)
                axis([0.0,30.0,0.0,50.0])
                title('HW Dead Tree Height Distribution (> 2.5 m)    Year = ' + str(yr + 1))
                xlabel('m')
                ylabel('no. HW trees')
                
            except:
                pass
            self.ycount = 0
            
        else: # nothing graphs
            return

    def multi_display(self, n=2, yr = 5):
        

        ''' displaying multiple graphs '''

        m = self.displayed_multi
        
        subplot(221)
        self.displayed = m[0]
        self.display_stand()
        if self.c == 0:
            colorbar()

        subplot(222)
        self.displayed = m[1]
        self.display_stand()
        if self.c == 0:
            colorbar()       

        subplot(223)
        self.displayed = m[2]
        self.display_stand()
        xlabel('Year = ' + str(yr))
        if self.c == 0:
            colorbar()

        subplot(224)
        self.displayed = m[3]
        self.display_stand()
        if self.c == 0:
            colorbar()  
