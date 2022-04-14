'''    LLM: longleaf pine - hardwood simulation model      '''

#  developed by Louise Loudermilk and Wendell P. Cropper, Jr. c. 2008-2010
#  see Loudermilk Dissertation (2010) and Loudermilk et al. (2011) for details
#  citation: Loudermilk, E. L., W. P. Cropper Jr, R. J. Mitchell, and H. Lee. 2011. Longleaf pine (Pinus palustris) and hardwood dynamics in a fire-maintained ecosystem: A simulation approach. Ecological Modelling 222:2733-2750.

''' GENERAL COMPONENTS: 1 yr. time step, 5 m x 5 m spatial resolution, 1.56 ha extent,
longleaf pine and southern hardwoods modeled, up to 10 trees per cell,
seed masting (pine) and clonal spread (HW) functions, mortality: fire, 'natural causes', and competition (for pine only),
site or stand index drive asymptotic growth function, fire frequency is a probability,
each fire burns the entire landscape, fuel accumulation (pine litter, HW litter, wiregrass) - a function of tree density,
fire intensity is a function of fuel accumulation,
calibration parameters for: competition index (impacts growth and mortality),
tree density index (litter accumulation), fire mortality '''

'''model initialization: HW: random (10%) landscape of 0.5 m HWs,
LP: imports trees (pickle) from a previous LLM run, model stabilization: after ~ 100 years,
This model was originally calibrated for two longleaf pine research sites:
Ichauway (Newton, GA) & Ordway Swisher Biological Station (Melrose, FL) '''

import matplotlib
from matplotlib import pylab
from pylab import *
import LLM_display # script for visually graphing (displaying) LLM outputs
#import LLM_print_shell1 # script for printing out some model details
import random
import math
import copy
import pickle
import hsiscore_class as HSI
#import hsiscore as HSI #cython code
##import psyco
##psyco.full()

def count_trees(tree_age,tree_count,num):
    #count number of trees > num
    idh=np.nonzero(tree_age > num)
    num_of_trees=np.sum(tree_count[idh[0][:],idh[1][:]])
    return num_of_trees

class LLM:

    def __init__(self):
        # read HSI input file
        [self.lb, \
         self.ub, \
         self.bn, \
         self.xsize, \
         self.ysize, \
         self.sm] = np.loadtxt('HVI_inputs.txt',delimiter=',',usecols=(1,2,3,4,5,6), skiprows=7 ,unpack=True)
        self.noHWyear=-1 
        self.hw_gt_age=[]
        self.hw_tot_age=[]
        self.lp_gt_age=[]
        self.lp_tot_age=[]
        # user INPUT values
        self.HSI=HSI.hsi_score()
        self.verbose=1 #EJ modified
        self.version = 'LLM_final.py'
        self.dim = 25 # dimensions
        self.fintim = 3 # no. of time steps

        # asymptotic height, based on stand index (SI * 1.3)
        self.StandIndexLP = 45.5 # Ich: 45.5, Ord: 32.5; stand index for LP (m), Ich = 35, Ord = 25
        self.StandIndexHW =32.5 # Ich: 32.5, Ord: 26.0; stand index for HW (m), Ich = 25, Ord = 20
        
        ''' Calibration Parameters'''
        # parameter for 'competition index ' LP & HW, VERY sensitive parameter
        self.CI_constLP = 0.07 # BASE Ich:0.070 (w/less fire: 0.074), Ord: 0.075
         # used for calibration, lower value = higher LP growth, less comp mortality, less dispersion b/w age/ht         
        self.CI_constHW = 0.095 # BASE w/ fire = Ich:0.095, w/less fire = 0.08, Ord: 0.1, 0.07 (w/less fire) # constant for 'competition index HW', 
         # used for calibration, lower value = higher HW growth, less dispersion b/w age/ht, no mortality from comp on HW

        # parameter for 'tree density index LP', used for calibration, higher = lower litter accumulation
        self.TDI_constLP = 0.26 # Ich: 0.26, Ord: 0.26
        self.TDI_constHW = 0.26  # Ich: 0.26, Ord: 0.26
        # TDI_const - not a very sensitive constant, slightly influencing fire mortality rates

        # Fire Mortality parameter, used for calibration (0-1), lower value = lower mortality from fire
        self.fire_constHW = 0.66 # Ich: 0.66, Ord: 0.7
        self.fire_constLP = 0.99 # Ich: 0.99, Ord: 0.95

        ''' other constants '''
        self.fire_prob = 0.25 # yearly fire probability, Ich: 0.35, Ord: 0.25
        self.mast_prob = 0.15 # seed masting prob. for LP, based on Boyer's long-term data
        self.max_LPcount = 10.0 # max no. of LLP seedlings that can germinate in a cell (hence max no. of LP trees in cell)
        self.max_HWcount = 10.0 # max no. of HW sprouts that can establish in a cell (hence max no. of HW trees in cell)
        self.adult = self.StandIndexLP * 0.51 # adult size for LP
        self.subadult = self.StandIndexLP * 0.25 # subadult size for LP
        self.subadultHW = self.StandIndexHW * 0.25 # subadult size for HW (no adult size for HW)
        self.seeds_cone = 32.0 # no. of seeds/cone
        self.germ_prob = 0.03 # max germination prob. for LP seeds
        self.germ_size = 0.1
        self.clonal_rate = 0.05 # spread rate of HW sprouts
        self.size_LP = 10 # no. of size classes for LP
        self.size_HW = 7 # no. of size classes for HW
        self.critical_TSF = 20 # critical time since fire - impacts mort_fire() and germLP() 
        self.high_mort = 0.05 # max 'natural' mortality rate for HWs
        self.tsf_mort = 0.0 # replaced with years time since fire when fires occur
        self.tree_mature_age = 7
        self.readfireprobfromfile=1
        self.readmastprobfromfile=1
        self.area_of_HW_resprouting=.1
        random.seed()
        self.rand_value=[]
        self.randint_value=[]
        self.randpos_value=[]
        self.rcount=0
        self.rcount_int=0
        self.rcount_pos=0
        self.randfromfile=0

    def instantiate(self,readf):
        # LLP & HW inputs
        # readf: 0 do not read from file

        self.old_age = zeros((self.dim,self.dim), 'float')+1
        self.old_ht = zeros((self.dim,self.dim), 'float')+2 # FOR NOW, assume all trees in cell are same ht and age!! 
        self.old_LPcount = zeros((self.dim,self.dim), 'float')+1 # tracking no. of all LP trees in cell
        ''' LLP initial values'''
        self.hts = []
        self.ages = [] # delete, to record age of 1.5 m LP
        self.new_age = zeros((self.dim,self.dim), 'float')
        self.new_ht = zeros((self.dim,self.dim), 'float')
        self.dead_htLP = zeros((self.dim,self.dim), 'float') # tracking dead trees
        self.dead_htLPc = zeros((self.dim,self.dim), 'float')
        self.dead_htLPn = zeros((self.dim,self.dim), 'float')
        self.dead_htLPf = zeros((self.dim,self.dim), 'float')
        self.seeds = zeros((self.dim,self.dim), 'float')
        self.germ_seeds = zeros((self.dim,self.dim), 'float')
        self.dorm_tag = zeros((self.dim,self.dim), 'float')
        self.new_LPcount = zeros((self.dim,self.dim), 'float')
        self.dead_ctLP = zeros((self.dim,self.dim), 'float') # dead, total
        self.dead_ctLPc = zeros((self.dim,self.dim), 'float') # dead from comp
        self.dead_ctLPn = zeros((self.dim,self.dim), 'float') # dead from nat
        self.dead_ctLPf = zeros((self.dim,self.dim), 'float') # dead from fire
        self.sum_htLP1 = zeros((self.dim,self.dim), 'float')
        self.sum_htLP2 = zeros((self.dim,self.dim), 'float')
        self.litter = zeros((self.dim,self.dim), 'float')
        self.burntlit = zeros((self.dim,self.dim), 'float')
        self.tot_biomass = zeros((self.dim,self.dim), 'float') # tracking total biomass (HW, LP, WG) after burn
        self.old_litter = zeros((self.dim,self.dim), 'float') # added for final map     
        self.old_seeds = zeros((self.dim,self.dim), 'float')

        self.Pg = zeros((self.dim,self.dim), 'float') # a cell's germination prob.
        self.fireflag = zeros((self.dim,self.dim), 'float') # 0 = no fire, 1 = fire
        self.time_since_fire = zeros((self.dim,self.dim), 'float') # tracking time since fire by cell: for germ prob., litter
        self.mortflag = zeros((self.dim,self.dim)) # '0' = has not been checked for mortality; 1 = checked for mortality
        self.lp_dbh = zeros((self.dim,self.dim), 'float') #Elchin added dbh
        self.lp_CR = zeros((self.dim,self.dim), 'float') #Elchin added crown radius

        ''' Wiregrass initial values '''
        self.litterWG = zeros((self.dim,self.dim), 'float') 
        self.old_litterWG = zeros((self.dim,self.dim), 'float') # added for final map
        self.burntlitWG = zeros((self.dim,self.dim), 'float')
        self.max_wtWG = 6.25 # max weight of wiregrass within a cell - kg/cell, asymptote
        
        ''' HW initial values '''
        self.old_HWcount = zeros((self.dim,self.dim), 'float')
        ##        self.old_HWcount = pickle.load(f4) # tracking no. of all HW trees in cell
        self.new_HWcount = zeros((self.dim,self.dim), 'float')
        self.dead_ctHW = zeros((self.dim,self.dim), 'float')
        self.dead_ctHWn = zeros((self.dim,self.dim), 'float')
        self.dead_ctHWf = zeros((self.dim,self.dim), 'float')
        self.surv_smHW_fire = [] # survive fire < 2m, HW
        self.mort_smHW_fire = [] # survive fire < 2m, HW       
        ##        self.old_ageHW = pickle.load(f5)
        self.old_ageHW = zeros((self.dim,self.dim), 'float')
        self.new_ageHW = zeros((self.dim,self.dim), 'float')
        ##        self.old_htHW = pickle.load(f3)
        self.old_htHW = zeros((self.dim,self.dim), 'float')
        self.new_htHW = zeros((self.dim,self.dim), 'float')
        self.dead_htHW = zeros((self.dim,self.dim), 'float')
        self.dead_htHWn = zeros((self.dim,self.dim), 'float')
        self.dead_htHWf = zeros((self.dim,self.dim), 'float')
        self.sum_htHW1 = zeros((self.dim,self.dim), 'float') # HW
        self.sum_htHW2 = zeros((self.dim,self.dim), 'float') # HW
        self.litterHW = zeros((self.dim,self.dim), 'float')
        self.old_litterHW = zeros((self.dim,self.dim), 'float') # for final map
        self.burntlitHW = zeros((self.dim,self.dim), 'float') # burnt litter in cell
        self.mortflagHW = zeros((self.dim,self.dim)) # '0' = has not been checked for mortality
        self.total_burntlit = zeros((self.dim,self.dim), 'float') # burnt litter in cell from both HW & LLP & WG, used in mort_fire
        self.hw_dbh = zeros((self.dim,self.dim), 'float') # Elchin added dbh
        self.hw_CR = zeros((self.dim,self.dim), 'float') # Elchin added crown radius

        ''' time-elapsed graphs '''

        self.D = LLM_display.display() # file.class name
        self.display_stand = self.D.display_stand # function within file.class
        self.multi_display = self.D.multi_display
        self.D.fintim = self.fintim

        # can only display single or multiple (time-elapsed) graphs during a run, not both
        # slows down processing time signficantly

        # single graphs, also see LLM_display.py
        self.D.displayed = 'LPcount' # choose which spatial variable to plot
        self.graph_yn = 'off'  ##### set from on to off to stop real time graphics #####
        if self.graph_yn == 'off': self.D.displayed = 'none'
        self.te = 0
        self.D.ycount = 5 # for histogram graphing

        # multiple graphs, this runs a bit slower than the single graphs, but good for visualization through time
        self.age_graph = ['age', 'ageHW', 'ht', 'htHW']
        self.cnt_graph = ['LPcount', 'HWcount', 'age', 'ageHW']
        #self.hist_graph = ['tree_histLP', 'dtree_histLP', 'tree_histHW', 'dtree_histHW']# not working yet
        self.lit_graph = ['litter', 'litterHW', 'WG','biomass']
        self.comp_graph = ['ht','htHW', 'WG', 'LPcount']

        self.D.displayed_multi = self.comp_graph ##### chose 4 graph combo above
        self.graph_yn_m = 'off' # on or off
        if self.graph_yn_m == 'off': self.D.displayed_multi = 'none'
        self.D.c = 0 # for colorbar
        self.D.max_wtWG = self.max_wtWG

        ''' other starting values '''
        self.mast_yr = 'no' # mast year
        self.copius_yr = 'no' # copius litter drop year
        self.recruit = 'yes'
        self.dead_comp = 0 
        self.dead_nat = 0
        self.dead_natHW = 0
        self.dead_fire = 0
        self.dead_fireHW = 0 # HW
        self.TK_HW = 0 # HW, top-killed count
        self.surv1 = 0 # for survival < 2m, HW FIRE
        self.dead1 = 0 # for mort < 2m, HW FIRE
        self.alive1 = 0

        ''' ending graphs, etc. '''

        self.graph_lg = self.subadult # tree size restrictions for outputs graphs
        self.graph_lgHW = self.subadultHW
        self.graph_sm = 5
        self.tracking_CI = []
        self.tracking_CI2 = []

        self.disperse = [] 
        self.tree_number = [] # adults
        self.HWtree_number = []
        self.tree_density = []
        self.HWtree_density = []
        self.relA_LP_all = []
        self.relA_HW_all = []
        self.HWdead = [] # keeping track of total HWs dead
        self.LPdead = []
        self.tracking_HWsprouts = []
        self.total_HWsprouts = []
        self.tracking_deadnatLP = [] # tracking dead by nat, fire, comp (in that order)
        self.tracking_deadfireLP = []
        self.tracking_deadcompLP = [] 
        self.tracking_deadnatHW = []
        self.tracking_deadfireHW = []
        self.tracking_seedlingsLP = []
        self.tracking_masting = []
        self.tracking_burnyrs = []
        self.tracking_mort_rates = []
        self.tracking_mortc_rates = []
        self.tracking_mortn_rates = []
        self.tracking_mortf_rates = []
        self.tracking_morth_rates = [] # HW
        self.tracking_mortnh_rates = []
        self.tracking_mortfh_rates = []
        self.max_yrht = []
        self.max_yrhtHW = []
        self.tlt8 = []
        self.tlt8HW = []   #number less than 8 m ht. (or whichever cut-off height is chosen)
        self.burn_area =[]  # number of cells burned in a time step (currently, all cells)
        self.all_cones_sub = [] # cones per tree through time
        self.all_cones_ad = []
        self.class_LPcountsA = [0] * self.size_LP # size-class distribution for alive LP
        self.class_LPcountsD = [0] * self.size_LP # size-class distribution for dead LP
        self.class_LPcountsDc = [0] * self.size_LP
        self.class_LPcountsDn = [0] * self.size_LP
        self.class_LPcountsDf = [0] * self.size_LP

        self.class_HWcountsA = [0] * self.size_HW # size-class distribution for alive HW
        self.class_HWcountsD = [0] * self.size_HW # size-class distribution for dead HW
        self.class_HWcountsDn = [0] * self.size_HW
        self.class_HWcountsDf = [0] * self.size_HW
        #score variables
        self.hw_sc = []
        self.age_sc = []
        self.hwHW_sc = []
        self.ageHW_sc = []
        self.dbh_sc = []
        self.ba_sc = []
        self.wire_sc = []
        self.gt_sc = []
        self.sq_sc = []

        self.graph_set()
        # EJ added beg_stand here which initialize LP and HW trees when we start over
        # previoulsy it was at the beginning of the run function 

        
        if self.randfromfile:
           #self.randnumberfromfile=np.loadtxt('LLM_random_number.txt')
           #self.randintnumberfromfile=int(np.loadtxt('LLM_random_int_number.txt'))
           
           #intrandfile = open("LLM_random_int_number.txt", "r")
           #self.randintnumberfromfile = []
           #for val in intrandfile.read().split():
           #  self.randintnumberfromfile.append(int(val))
           #intrandfile.close()

           #intrandfile = open("LLM_random_pos.txt", "r")
           #self.randposnumberfromfile = []
           #for val in intrandfile.read().split():
           #  self.randposnumberfromfile.append(int(val))
           #intrandfile.close()

           f1 = open('LLM_random_number.txt')   # used random HW sprouts, see beg_stand()
           f2 = open('LLM_random_int_number.txt')
           f3 = open('LLM_random_pos.txt')
           self.randnumberfromfile = pickle.load(f1)
           self.randintnumberfromfile = pickle.load(f2)
           self.randposnumberfromfile = pickle.load(f3) # input age data (from 100 yr. run)
           f1.close()
           f2.close()
           f3.close()

           print ('RANDOMFILEREAD', np.shape(self.randnumberfromfile))
           print ('RANDOMINTFILEREAD', np.shape(self.randintnumberfromfile))
           print ('RANDOMINTFILEREAD', np.shape(self.randposnumberfromfile))
           print ('')

        #NOTE: order beg_stand has to go after readrandfile
        self.beg_stand()

        if readf:
           f = open('data/LP_ht.txt', 'rb') # from 100 year simulation (output from a previous version)
           f1 = open('data/LP_count.txt', 'rb')
           f2 = open('data/LP_age.txt', 'rb')
           f3 = open('data/lit_WG.txt', 'rb')
           f4 = open('data/HW_ht.txt', 'rb')   # used random HW sprouts, see beg_stand()
           f5 = open('data/HW_count.txt', 'rb')
           f6 = open('data/HW_age.txt', 'rb')
           f7 = open('data/lit_LLP.txt', 'rb')
           f8 = open('data/lit_HW.txt', 'rb')
        
           self.old_ht = pickle.load(f)
           self.old_LPcount = pickle.load(f1)
           self.old_age = pickle.load(f2) # input age data (from 100 yr. run)
           self.old_litterWG = pickle.load(f3)
           self.old_htHW = pickle.load(f4)
           self.old_HWcount = pickle.load(f5)
           self.old_ageHW = pickle.load(f6)
           self.litter = pickle.load(f7)
           self.litterHW = pickle.load(f8)
           f.close()
           f1.close()
           f2.close()
           f3.close()
           f4.close()
           f5.close()
           f6.close()
           f7.close()
           f8.close()
           print ('no. of starting LLP trees read from file:', np.sum(self.old_LPcount))
           print ('no. of starting HW trees read from file:', np.sum(self.old_HWcount))
           print ('Data read from files!')
   
    def save_randfromfile(self):
        if not self.randfromfile:
           #np.savetxt('LLM_random_number.txt', self.rand_value, fmt='%1.4f')
           #np.savetxt('LLM_random_int_number.txt', self.randint_value, fmt='%i')
           #np.savetxt('LLM_random_pos.txt', self.randpos_value, fmt='%i')
           
           pickle.dump( self.rand_value, open( "LLM_random_number.txt", "wb" ) )
           pickle.dump( self.randint_value, open( "LLM_random_int_number.txt", "wb" ) )
           pickle.dump( self.randpos_value, open( "LLM_random_pos.txt", "wb" ) )
           
           print ('randfromfile SAVED', np.shape(self.rand_value))
           print ('LLM_random_int_number SAVED', np.shape(self.randint_value))
           print ('LLM_random_pos SAVED', np.shape(self.randpos_value))

    def save_mast_fire_to_files(self):
        np.savetxt('LLM_fire_years.out', self.myears, fmt='%i') 
        np.savetxt('LLM_mast_years.out', self.fyears, fmt='%i') 

    def save_pickle(self):
        pickle.dump( self.old_ht, open( "data/LP_ht.txt", "wb" ) )
        pickle.dump( self.old_LPcount, open( "data/LP_count.txt", "wb" ) )
        pickle.dump( self.old_age, open( "data/LP_age.txt", "wb" ) )
        pickle.dump( self.old_htHW, open( "data/HW_ht.txt", "wb" ) )
        pickle.dump( self.old_HWcount, open( "data/HW_count.txt", "wb" ) )
        pickle.dump( self.old_ageHW, open( "data/HW_age.txt", "wb" ) )
        pickle.dump( self.litterHW, open( "data/lit_HW.txt", "wb" ) )
        pickle.dump( self.litter, open( "data/lit_LLP.txt", "wb" ) )
        pickle.dump( self.old_litterWG, open( "data/lit_WG.txt", "wb" ) )

   
    def graph_set(self):

        # see file 'LLM_display.py' for details

        self.D.old_ht = self.old_ht
        self.D.old_age = self.old_age
        self.D.litter = self.litter
        self.D.litterHW = self.litterHW
        self.D.litterWG = self.litterWG
        self.D.tot_biomass = self.tot_biomass
        self.D.seeds = self.seeds
        self.D.germ_seeds = self.germ_seeds
        self.D.new_ht = self.new_ht 
        self.D.sum_htLP1 = self.sum_htLP1
        self.D.old_LPcount = self.old_LPcount
        self.D.old_HWcount = self.old_HWcount
        self.D.fireflag = self.fireflag
        self.D.mortflag = self.mortflag
        self.D.old_htHW = self.old_htHW
        self.D.old_ageHW = self.old_ageHW
        self.D.new_LPcount = self.new_LPcount
        self.D.new_HWcount = self.new_HWcount
        self.D.dead_htLP = self.dead_htLP
        self.D.dead_htHW = self.dead_htHW
        self.D.dead_ctHW = self.dead_ctHW
        self.D.dead_ctLP = self.dead_ctLP
        self.D.new_htHW = self.new_htHW
        
    def rand_calc(self): # random no. 0-1

        '''calculate random values from 0-1'''
        
        if self.randfromfile:
           self.rnd_value = self.randnumberfromfile[self.rcount]
           self.rcount+=1
           #print self.rcount
        else:
           self.rnd_value = random.random()
           self.rand_value.append(self.rnd_value)


    def beg_stand(self): # to make tree landscape

        ''' to initialize stand '''
##
        no_old_ht = 0
        for row in range(self.dim):
            for col in range(self.dim):
                self.rand_calc()
                if self.rnd_value < 0.20:
                    no_old_ht = no_old_ht + 1
                    self.old_ht[row][col] = 20.0
                    self.old_age[row][col] = 100.0
                    self.old_LPcount[row][col] = 1.0
        
        print ('no. of starting LLP trees', no_old_ht)
        
        no_old_ht2 = 0
        for row in range(self.dim):
            for col in range(self.dim):
                self.rand_calc()
                if self.rnd_value < 0.10: # random landscape of hardwoods, within ~10% of the area (cells)
                    no_old_ht2 = no_old_ht2 + 1
                    self.old_htHW[row][col] = 0.5
                    self.old_ageHW[row][col] = 1.0 # difficult to determine age of young HW with frequent fires and resprouting
                    self.old_HWcount[row][col] = 1.0
        
        print ('no. of starting HW trees', no_old_ht2)


    ''' HARDWOOD section '''


    def reproHW(self, row,col): # done

        ''' Hardwood reproduction, mainly by clonal spread, only seed-out with adults '''

        '''INPUTS: old_ht, rnd_value, ring_id, old_HWcount;
        OUTPUTS: new_htHW, new_agHW, new_HWcount'''

        trackingHW = 0
        
        self.rand_calc()
        if self.old_htHW[row][col] > 0.0:  # if HW is present in cell
            self.ring_id(row,col,1)
            len_c = len(self.cells) # should be 8
            
            if self.randfromfile:
               pos = self.randposnumberfromfile[self.rcount_pos]
               self.rcount_pos+=1
            else:
               pos = random.randint(0,len_c-1)
               self.randpos_value.append(pos)
            #print pos 
            #pos = random.randint(0,len_c-1)

            # added by EJ
            cs = self.cells[pos] # randomly chosen cell for clonal spread
            if self.old_htHW[cs[0],cs[1]] == 0.0: # if no HW present
                if self.rnd_value < self.clonal_rate: #  chance of establishing

                    if self.randfromfile:
                       r = self.randintnumberfromfile[self.rcount_int]
                       self.rcount_int+=1
                    else:
                       r = random.randint(1,self.max_HWcount) # up to 'x' sprouts
                       self.randint_value.append(r) 

                    # added by EJ
                    self.new_htHW[cs[0],cs[1]] = 0.5 # new sprout establishes
                    self.new_ageHW[cs[0],cs[1]] = 1.0
                    self.new_HWcount[cs[0],cs[1]] = self.old_HWcount[cs[0],cs[1]] + float(r)
                    trackingHW = trackingHW + self.new_HWcount[cs[0],cs[1]]
                    
            else: pass

            self.tracking_HWsprouts.append(trackingHW)

        else: pass          
    

    def HtSumHW(self, row, col): # done

        ''' calculates sum of adult and subadult tree heights OR heights taller
        than the focal tree in the 2 cell neighborhood, adjusts for no. of trees within cell, both HW & LP
        impacts HW growth'''

        '''INPUTS: ring cells, old_ht, subadult; OUTPUTS: sum_htHW1 or sum_htHW2 (incorporating LPcounts) '''

        ''' Height Sum for Competition Index - IndexHW()'''

        self.ring_id(row, col, 1)
        self.cells_h1 = self.cells
        self.ring_id(row, col, 2)
        self.cells_h2 = self.cells + self.cells_h1
        
        if self.old_htHW[row,col] > 0.0: # if a tree is present
            for k in self.cells_h2:
                if self.old_htHW[row,col] <= self.old_ht[k[0],k[1]] or self.old_ht[k[0],k[1]] >= self.subadult: # effects from LP
                    self.sum_htHW1[row,col] = self.sum_htHW1[row,col] + self.old_ht[k[0],k[1]] * math.exp(self.old_LPcount[k[0],k[1]]/self.max_LPcount)# adjusts for no. LP trees present
                if self.old_htHW[row,col] <= self.old_htHW[k[0],k[1]] or self.old_htHW[k[0],k[1]] >= self.subadultHW:
                    self.sum_htHW1[row,col] = self.sum_htHW1[row,col] + self.old_htHW[k[0],k[1]] * math.exp(self.old_HWcount[k[0],k[1]]/self.max_HWcount)# adjusts for no. HW trees present
        if self.sum_htHW1[row,col] > self.old_htHW[row,col]: # restricts for neg. sum-height values
            self.sum_htHW1[row,col] = self.sum_htHW1[row,col] - self.old_htHW[row,col]# still allows for sum_ht values less than old_ht value

        ''' Height Sum for Tree Density Index - IndexHW(), used for biomass accumulation of HW litter '''
        
        # focal cell trees
        self.sum_htHW2[row][col] = (self.old_htHW[row][col] * math.exp(self.old_HWcount[row][col]/self.max_HWcount))*3
        # *3 = more influence for litter in focal cell than surrounding cells
 
        for k in self.cells_h2: # nearby trees, subadult and adult HW only
            if self.old_htHW[k[0],k[1]] >= self.subadultHW:
                self.sum_htHW2[row,col] = self.sum_htHW2[row,col] + self.old_htHW[k[0],k[1]] * math.exp(self.old_HWcount[k[0],k[1]]/self.max_HWcount)# adjusts for no. HW trees present
        
      
    def IndexHW(self, row, col): # done

        '''INPUTS: sum_htHW1 or sum_htHW2, maxht, CI,Nr; OUTPUTS: CI_HW or TDI_HW '''

        # Competition Index (LP & HW), impacts HW growth, LP growth and mortality from comp.
            
        self.CI_HW = self.CI_constHW * self.sum_htHW1[row,col]/self.StandIndexHW
        
        if self.CI_HW > 1.0: self.CI_HW = 1.0 # eliminates CI > 1 values
        
        if self.CI_HW > 0.6: # tracking CI values
            self.tracking_CI.append(self.old_htHW[row][col]) 
            self.tracking_CI2.append(self.CI_HW)
                
        # Tree Density Index (HW), impacts HW fuel distribution and accumulation
        
        self.TDI_HW = self.sum_htHW2[row][col]/self.StandIndexHW * self.TDI_constHW
        if self.TDI_HW > 1.0: self.TDI_HW = 1.0 # eliminates TDI > 1 values

        
    def HWgrowth(self, row, col):

        ''' HW growth.  Equation from Greenberg & Simons, age/ht relationship, Fig. 5 QULA'''

        if self.mortflagHW[row][col] == 1: return # if dead or top-killed, no (additional) growth

        if self.old_ageHW[row][col] > 0.0 and self.old_ageHW[row][col]<= 3.0:# younger HWs
                
            growthHW = 0.5 # linear growth for 1st 3 years, no effects from 'CI'
            
            self.new_htHW[row][col] = self.old_htHW[row][col] + growthHW
            self.new_ageHW[row][col] = self.old_ageHW[row][col] + 1.0
                
        elif self.old_ageHW[row][col] > 3.0 and self.old_htHW[row][col] < self.StandIndexHW: # older HWs
            self.IndexHW(row, col)

            ht = self.StandIndexHW * (1.0 - math.exp(-0.02 * self.old_ageHW[row][col]))
            ht2 = self.StandIndexHW  * (1.0 - math.exp(-0.02 * (self.old_ageHW[row][col] + 1.0)))  # same as LP
            growthHW = (ht2 - ht) * (1.0 - self.CI_HW)

            self.new_htHW[row][col] = self.old_htHW[row][col] + growthHW
            self.new_ageHW[row][col] = self.old_ageHW[row][col] + 1.0

        elif self.old_htHW[row][col] >= self.StandIndexHW:

            self.new_htHW[row][col] = self.StandIndexHW # no ht growth, but still ages 
            self.new_ageHW[row][col] = self.old_ageHW[row][col] + 1.0            

        else:
            return # no hardwood to grow
       

    def litHW(self,row,col): # working
        
        ''' Hardwood litter is a function of biomass (from literature and field data)
        measured as kg/cell (or 25 sq. m), done as biomass increment accumulation'''

        if self.time_since_fire[row][col] > 0.0: # Need to check/verify more!!
            self.IndexHW(row,col)

            bio = 0.5 * (self.time_since_fire[row][col] - 1) 
            bio2  = 0.5 * self.time_since_fire[row][col]
            bio_inc = (bio2 - bio) * self.TDI_HW
               
            self.litterHW[row][col] = self.litterHW[row][col] + bio_inc
            
        else:
            return
    
    def mort_compHW(self): # inactive function

        ''' we assume there is no direct mortality on HW from competition w/ LP, 
            only growth suppression from LP '''
        pass


    def mort_natHW(self,row, col):
        
        # linearly increasing, caps at 0.05 

        dead_cell = 0
        s = 0.002 
                
        if self.mortflagHW[row][col] == 0:
            if self.old_htHW[row][col] > 0.0: # checks all present HWs
                for i in range(int(self.old_HWcount[row][col])):
                    self.rand_calc()
                    prob_nat = s * self.old_htHW[row][col]
                    if prob_nat > 0.05: # mortality limited to 0.05
                        prob_nat = 0.05
                    if self.rnd_value < prob_nat:
                     
                        self.dead_natHW = self.dead_natHW + 1
                        dead_cell = dead_cell + 1
                    
                if dead_cell == self.old_HWcount[row][col]:
                    self.new_ageHW[row][col] = 0.0 # all trees die and age & ht is reset to 0
                    self.new_htHW[row][col] = 0.0
                    self.new_HWcount[row][col] = 0.0
                    self.mortflagHW[row][col] = 1.0
                    self.dead_htHW[row][col] = self.old_htHW[row][col] # tracking dead HW hts
                    self.dead_htHWn[row][col] = self.old_htHW[row][col]
                    
                else: # if not all trees in cell died
                    self.new_HWcount[row][col] = self.old_HWcount[row][col] - dead_cell # removes dead trees
                    self.old_HWcount[row][col] = self.new_HWcount[row][col] # these trees will be checked for other mortality                    
                    self.dead_htHW[row][col] = self.old_htHW[row][col]
                    self.dead_htHWn[row][col] = self.old_htHW[row][col]
                    
                self.dead_ctHW[row][col] = dead_cell # tracking no. dead HWs w/in cells
                self.dead_ctHWn[row][col] = dead_cell

    def mort_fireHW(self,row,col):

        # function of: time since fire (fuel accum.), surrounding fuels (both affecting fire intensity), size (ht) of HW
        # use data collected from Ichauway (for smaller HWs at least): threshold ht ~ 2 m for re-sprouting
        # could use crown survival/mortality from Rebertus EA 1989 for QULA ?
        
        dead_cell = 0

        tot_burntlit = self.total_burntlit # HW & LLP & WG burned litter in cell, see fire_spread()
        
        if self.fireflag[row][col] == 1: # if fires occurred
            if self.mortflagHW[row][col] == 0: # tree has not been killed yet
                if 0.0 < self.old_htHW[row][col] < 2.0: # HW present and < 2m ht, (also < 4yrs. old since establishment or last fire)

                    self.alive1 = self.alive1 + self.old_HWcount[row][col]
                    #print 'total alive', self.alive1
                    self.rand_calc()
                    r1 = self.rnd_value
                    self.rand_calc()
                    r2 = self.rnd_value
                    if r1 < 0.02: # 2% pure survival prob. ( 2% from Ichauway data)
                        self.surv1 = self.surv1 + self.old_HWcount[row][col]
                        return # nothing changes, all HWs in cell survive

                    elif r2 < 0.008: # 0.6% pure mortality prob. (0.6% from Ichauway data)
                        self.dead1 = self.dead1 + self.old_HWcount[row][col]
                        self.dead_htHW[row][col] = self.old_htHW[row][col]
                        self.dead_htHWf[row][col] = self.old_htHW[row][col]
                        
                        self.new_htHW[row][col] = 0.0
                        self.new_ageHW[row][col] = 0.0
                        self.mortflagHW[row][col] = 1
                        dead_cell = dead_cell + self.old_HWcount[row][col]
                        self.dead_fireHW = self.dead_fireHW + dead_cell
                        
                    else: # top-killed otherwise
                        self.TK_HW = self.TK_HW + self.old_HWcount[row][col]
                        self.new_htHW[row][col] = 0.5
                        self.new_ageHW[row][col] = 1.0 # resets age to 1.0
                        self.mortflagHW[row][col] = 1

                    ''' incorporate mortality prob? Rebertus 1989
                        incorporate survival prob farther from adult LPs? Williamson & Black '''
                        
                elif self.old_htHW[row][col] >= 2.0:     # HW present > 2m ht

                    fire_int = 1.0 - math.exp(-tot_burntlit[row][col])
                    size_prob = math.exp(0.2 * (-self.old_htHW[row][col] + 1.0))
                    mf = (fire_int * size_prob) * self.fire_constHW
                    
                    for i in range(int(self.old_HWcount[row][col])): # check all trees in cell
                        self.rand_calc()
                        if self.rnd_value < mf: # tree dies
                            dead_cell = dead_cell + 1

                    if dead_cell == self.old_HWcount[row][col]:                       
                        self.dead_fireHW = self.dead_fireHW + dead_cell
                        self.new_htHW[row][col] = 0.0  
                        self.new_ageHW[row][col] = 0.0  
                        self.mortflagHW[row][col] = 1
                        self.dead_htHW[row][col] = self.old_htHW[row][col]
                        self.dead_htHWf[row][col] = self.old_htHW[row][col]

                    else: # if not all trees in cell died
                        self.dead_fireHW = self.dead_fireHW + dead_cell 
                        self.new_HWcount[row][col] = self.old_HWcount[row][col] - dead_cell # removes dead trees, but hts and age remain the same for remaining trees
                        self.old_HWcount[row][col] = self.new_HWcount[row][col] # trees can be checked for other mortality                        
                        self.dead_htHW[row][col] = self.old_htHW[row][col]
                        self.dead_htHWf[row][col] = self.old_htHW[row][col]
                                               
                else: pass # no HWs in cell

                self.dead_ctHW[row][col] = dead_cell # tracking no. dead HWs w/in cells
                self.dead_ctHWf[row][col]  = dead_cell

                
    ''' WIREGRASS Section '''


    def litWG(self, row, col):

        ''' wiregrass is only modeled as a fuel component, not competition, etc. '''

        if self.time_since_fire[row][col] > 0.0: 
            
            self.IndexHW(row,col)
            self.IndexLP(row,col)
            self.TDI_WG = (self.TDI_HW + self.TDI_LP)/2 
            if self.TDI_WG > 1.0: self.TDI_WG = 1.0

            bio = self.max_wtWG * (1.0 - math.exp( -0.5 * (self.time_since_fire[row][col]-1)))
            bio2  = self.max_wtWG * (1.0 - math.exp( -0.5 * self.time_since_fire[row][col]))
            bio_inc = (bio2 - bio) * (1.0 - self.TDI_WG)
    
            self.litterWG[row][col] = self.litterWG[row][col] + bio_inc

            
        else:
            return
              

    ''' LONGLEAF PINE Section '''

    def seed_intLP(self):
        
        ''' LLP seeding intensity function '''

        if self.randfromfile:
           self.seed_int = self.randintnumberfromfile[self.rcount_int]
           self.rcount_int+=1
        else:
           self.seed_int = random.randint(1,4) # level of seeding/masting intensity (1-4)
           self.randint_value.append(self.seed_int)



        # added by EJ 

    def seed_prodLP(self):

        ''' LLP seed dispersal intensity function: for mast or normal year
        and 4-levels of varying cone-production intensity '''

        ''' INPUTS: disperse, seeds_int, seeds_cone, mast_yr OUTPUTS: self.disperse '''

        self.seed_intLP()  
        if self.mast_yr == 'yes':
            
            if self.seed_int == 1:
                if self.randfromfile:
                   cones_tree1 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree1 = random.randint(10,20) # subadult cones/tree
                   self.randint_value.append(cones_tree1)
                # added by EJ
                if self.randfromfile:
                   cones_tree2 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree2 = random.randint(40,50) # subadult cones/tree
                   self.randint_value.append(cones_tree2)
                # added by EJ
                
            elif self.seed_int == 2:
                if self.randfromfile:
                   cones_tree1 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree1 = random.randint(20,30) # subadult cones/tree
                   self.randint_value.append(cones_tree1)
                # added by EJ
                if self.randfromfile:
                   cones_tree2 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree2 = random.randint(50,70) # subadult cones/tree
                   self.randint_value.append(cones_tree2)
                # added by EJ
                
            elif self.seed_int == 3:
                if self.randfromfile:
                   cones_tree1 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree1 = random.randint(30,50) # subadult cones/tree
                   self.randint_value.append(cones_tree1)
                # added by EJ
                if self.randfromfile:
                   cones_tree2 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree2 = random.randint(70,90) # subadult cones/tree
                   self.randint_value.append(cones_tree2)
                # added by EJ
                
            else: # = 4
                if self.randfromfile:
                   cones_tree1 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree1 = random.randint(50,70) # subadult cones/tree
                   self.randint_value.append(cones_tree1)
                # added by EJ
                if self.randfromfile:
                   cones_tree2 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree2 = random.randint(90,125) # subadult cones/tree
                   self.randint_value.append(cones_tree2)
                # added by EJ
                
        else: # normal seed year
            
            if self.seed_int == 1:
                if self.randfromfile:
                   cones_tree1 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree1 = random.randint(0,5) # subadult cones/tree
                   self.randint_value.append(cones_tree1)
                # added by EJ
                if self.randfromfile:
                   cones_tree2 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree2 = random.randint(0,10) # subadult cones/tree
                   self.randint_value.append(cones_tree2)
                # added by EJ
                
            elif self.seed_int == 2:
                if self.randfromfile:
                   cones_tree1 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree1 = random.randint(5,10) # subadult cones/tree
                   self.randint_value.append(cones_tree1)
                # added by EJ
                if self.randfromfile:
                   cones_tree2 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree2 = random.randint(10,20) # subadult cones/tree
                   self.randint_value.append(cones_tree2)
                # added by EJ
                
            elif self.seed_int == 3:
                if self.randfromfile:
                   cones_tree1 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree1 = random.randint(10,15) # subadult cones/tree
                   self.randint_value.append(cones_tree1)
                # added by EJ
                if self.randfromfile:
                   cones_tree2 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree2 = random.randint(20,30) # subadult cones/tree
                   self.randint_value.append(cones_tree2)
                # added by EJ
                
            else: # = 4
                if self.randfromfile:
                   cones_tree1 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree1 = random.randint(15,20) # subadult cones/tree
                   self.randint_value.append(cones_tree1)
                # added by EJ
                if self.randfromfile:
                   cones_tree2 = self.randintnumberfromfile[self.rcount_int]
                   self.rcount_int+=1
                else:
                   cones_tree2 = random.randint(30,40) # subadult cones/tree
                   self.randint_value.append(cones_tree2)
                # added by EJ

        self.all_cones_sub.append(cones_tree1) # store for graph
        self.all_cones_ad.append(cones_tree2)
        d1 = cones_tree1 * self.seeds_cone
        d2 = cones_tree2 * self.seeds_cone
        self.disperse.append(d1)
        self.disperse.append(d2) # stores no. seeds per tree for both subadult & adult
    
    def seed_dispLP(self, row, col): # 
        
        ''' LLP seed dispersal during mast and normal years, assumed 
        only adults and subadults produce seeds, seeds only dispersed in cells
        with no tree present'''

        '''INPUTS: mast_yr, old_ht, disperse, ring cells, adult, subadult;
        OUTPUTS: seeds'''
              
        if self.mast_yr == 'yes': # if mast year
            
            if self.old_ht[row][col] >= self.adult: # adult present, mast year

                div = self.disperse[1] / 80.0               
                if div == 0.0: return # if no seeds, return
                elif div < 1: div = 1 # have at least 1 seed/cell if seeds are produced
                else: pass
                div = int(div)

                self.ring_id(row,col,3) # 3nd ring around focal cell
                for k in self.cells: 
                    if self.old_ht[k[0],k[1]] == 0.0:
                        self.seeds[k[0],k[1]] = div + self.seeds[k[0],k[1]]  

                    else:
                        self.seeds[k[0],k[1]] = 0 + self.seeds[k[0],k[1]]
                        
                self.ring_id(row,col,2) # 2nd ring
                for k in self.cells: 
                    if self.old_ht[k[0],k[1]] == 0.0:
                        self.seeds[k[0],k[1]] = div * 2 + self.seeds[k[0],k[1]]  

                    else:
                        self.seeds[k[0],k[1]] = 0 + self.seeds[k[0],k[1]]  

                self.ring_id(row,col,1) # 1st ring
                for k in self.cells: 
                    if self.old_ht[k[0],k[1]] == 0.0:
                        self.seeds[k[0],k[1]] = div * 3 + self.seeds[k[0],k[1]]   

                    else:
                        self.seeds[k[0],k[1]] = 0 + self.seeds[k[0],k[1]]  
                        

            elif self.old_ht[row][col] < self.adult and self.old_ht[row][col] > self.subadult: # subadult present
                
                div = self.disperse[0] / 32.0
                if div == 0.0: return # if no seeds, return
                elif div < 1: div = 1
                else: pass
                div = int(div)

                self.ring_id(row,col,2) # 2nd ring
                for k in self.cells: 
                    if self.old_ht[k[0],k[1]] == 0.0:
                        self.seeds[k[0],k[1]] = div + self.seeds[k[0],k[1]]  

                    else:
                        self.seeds[k[0],k[1]] = 0 + self.seeds[k[0],k[1]]  

                self.ring_id(row,col,1) # 1st ring
                for k in self.cells: 
                    if self.old_ht[k[0],k[1]] == 0.0:
                        self.seeds[k[0],k[1]] = div * 2 + self.seeds[k[0],k[1]]   

                    else:
                        self.seeds[k[0],k[1]] = 0 + self.seeds[k[0],k[1]]
            else:
                return

                        
        else: # normal year
            
            if self.old_ht[row][col] >= self.adult: # normal year, adult tree
                div = self.disperse[1] / 32.0
                if div == 0.0: return # if no seeds, return
                elif div < 1: div = 1
                else: pass
                div = int(div)

                self.ring_id(row,col,2) # 2nd ring
                for k in self.cells: 
                    if self.old_ht[k[0],k[1]] == 0.0:
                        self.seeds[k[0],k[1]] = div + self.seeds[k[0],k[1]]  

                    else:
                        self.seeds[k[0],k[1]] = 0 + self.seeds[k[0],k[1]]  

                self.ring_id(row,col,1) # 1st ring
                for k in self.cells: 
                    if self.old_ht[k[0],k[1]] == 0.0:
                        self.seeds[k[0],k[1]] = div * 2 + self.seeds[k[0],k[1]]   

                    else:
                        self.seeds[k[0],k[1]] = 0 + self.seeds[k[0],k[1]]  

            elif self.old_ht[row][col] < self.adult and self.old_ht[row][col] > self.subadult: # normal year, subadult tree

                div = self.disperse[0] / 8.0
                if div == 0.0: return # if no seeds, return
                if 0 < div < 1: div = 1
                div = int(div)
                
                self.ring_id(row,col,1) # 1st ring
                for k in self.cells: 
                    if self.old_ht[k[0],k[1]] == 0.0:
                        self.seeds[k[0],k[1]] = div + self.seeds[k[0],k[1]]                      
 
                    else:
                        self.seeds[k[0],k[1]] = 0 + self.seeds[k[0],k[1]]
            else:
                return
                
    def germLP(self): # check time since fire influence on germination
        
        ''' Germination of LLP seeds, prob. is a function of time since fire'''
        
        '''INPUTS: seeds, time_since_fire, rnd_value, fireflag; OUTPUTS: LPcount, new_ht, new_age'''

        ''' starting height = 0.1 m, can establish up to 10 seedlings in a cell '''

        if self.recruit == 'yes':
            maxlp = self.max_LPcount
            #print 'shape(old_ht):',np.shape(self.old_ht)    
            for i in range(len(self.seeds)):
                for j in range(len(self.seeds)):
   
                    if self.old_ht[i,j] == 0.0 and self.old_htHW[i,j]== 0.0 and self.time_since_fire[i,j] < 10.0:
                        # no trees present and time since fire < 10 yrs
                            
                        if self.seeds[i][j] >= 1.0: # if seeds are in cell, up to 10 can germinate            
                            self.Pg[i][j] = self.germ_prob * math.exp(-1.0 * self.time_since_fire[i][j])
                            # germ prob, based on time since fire
                            
                            for k in range(int(self.seeds[i][j])): # check all seeds in cell for germination
                                if self.new_LPcount[i][j] < maxlp: # only max allowed to germinate                        
                                    self.rand_calc()
                                    
                                    if self.rnd_value < self.Pg[i][j]: # if germinated
                                        self.new_LPcount[i][j] = self.new_LPcount[i][j] + 1
                                        self.germ_seeds[i][j] = self.germ_seeds[i][j] + 1

                                else:
                                    break # exit loop, max count has been reached

                            if self.new_LPcount[i][j] > 0.0:  # if 1+ seeds germinated
                                if self.randfromfile:
                                   self.dorm_tag[i][j] = self.randintnumberfromfile[self.rcount_int]
                                   self.rcount_int+=1
                                else:
                                   self.dorm_tag[i][j] = random.randint(1,10) # seedlings remain 'dormant' for 1-10 years
                                   self.randint_value.append(self.dorm_tag[i][j])
                                # added by EJ
                                # they stay at 0.1 m and age in LPgrowth()
                                
        else:
            pass # no recruitment past tsf_mort
    
    def HtSumLP(self, row, col): # see adjustment for counts

        ''' calculates sum of adult and subadult tree heights OR heights taller
        than the focal tree in the 2 cell neighborhood, used in competition equation: eq. 5
        RA() function, CHANGED to adjust for no. of trees within cell, both HW & LP '''

        '''INPUTS: ring cells, old_ht, subadult; OUTPUTS: sum_ht (incorporating LPcounts) '''
    
        self.ring_id(row, col, 1)
        self.cells_h1 = self.cells
        self.ring_id(row, col, 2)
        self.cells_h2 = self.cells + self.cells_h1

        ''' Height Sum for Competition Index '''
        
        if self.old_ht[row,col] > 0.0: # if a tree is present
            for k in self.cells_h2:
                if self.old_ht[row,col] <= self.old_ht[k[0],k[1]] or self.old_ht[k[0],k[1]] >= self.subadult:
                    self.sum_htLP1[row,col] = self.sum_htLP1[row,col] + self.old_ht[k[0],k[1]] * math.exp(self.old_LPcount[k[0],k[1]]/self.max_LPcount)# adjusts for no. LP trees present
                if self.old_ht[row,col] <= self.old_htHW[k[0],k[1]] or self.old_htHW[k[0],k[1]] >= self.subadultHW:
                    self.sum_htLP1[row,col] = self.sum_htLP1[row,col] + self.old_htHW[k[0],k[1]] * math.exp(self.old_HWcount[k[0],k[1]]/self.max_HWcount)# adjusts for no. HW trees present
        if self.sum_htLP1[row,col] > self.old_ht[row,col]: # NEW, restricts for neg. sum-height values
            self.sum_htLP1[row,col] = self.sum_htLP1[row,col] - self.old_ht[row,col]# still allows for sum_ht values less than old_ht value


        ''' Height Sum for Tree Density Index, used for biomass accumulation of LP litter '''

        self.sum_htLP2[row][col] = (self.old_ht[row][col] * math.exp(self.old_LPcount[row][col]/self.max_LPcount))*3 # focal trees
        # *3 = more influence for litter in focal cell than surrounding cells

        for k in self.cells_h2: # nearby trees
            if self.old_ht[k[0],k[1]] >= self.subadult:
                self.sum_htLP2[row,col] = self.sum_htLP2[row,col] + self.old_ht[k[0],k[1]] * math.exp(self.old_LPcount[k[0],k[1]]/self.max_LPcount)# adjusts for no. LP trees present

        
    def IndexLP(self, row, col): 
        
        '''Competition from LP&HW: used for competition mortality, growth suppression,
        Tree Density Index of LP: used for litter accumulation of '''

        '''INPUTS: sum_htLP1 or sum_htLP2, maxht, CI,Nr; OUTPUTS: CI_LP or TDI_HW '''

        # Competition Index (LP & HW)
        
        self.CI_LP = self.CI_constLP * self.sum_htLP1[row,col]/self.StandIndexLP        
        
        if self.CI_LP > 1.0: self.CI_LP = 1.0 # eliminates neg. CI values

        # Tree Density Index (LP)

        self.TDI_LP = self.sum_htLP2[row,col]/ self.StandIndexLP * self.TDI_constLP

        if self.TDI_LP > 1.0: self.TDI_LP = 1.0 # eliminates TDI > 1 values           
            
        
    def litLP(self,row,col): 

        ''' Pine litter is a function of biomass (from literature and field data)
        measured as kg/cell (or 25 sq. m), done as biomass increment accumulation'''

        if self.time_since_fire[row][col] > 0.0: 
            self.IndexLP(row,col)

            if self.time_since_fire[row][col] == 1.0: 
                bio = 0.0 
            else:
                bio = 2.5 + 2.5 * (self.time_since_fire[row][col] -1)
            bio2  = 2.5 + 2.5 * self.time_since_fire[row][col]
            bio_inc = (bio2 - bio) * self.TDI_LP
             
            self.litter[row][col] = self.litter[row][col] + bio_inc
            
        else:
            return
            
    def LPgrowth(self, row, col):

        '''  LP growth from Platt & curve expert; VERY similar to Greenberg and Simons '99 '''

        self.IndexLP(row, col)

        if self.mortflag[row,col] == 1: return # if dead, no growth

        if self.dorm_tag[row,col] == 0.0: # if seedlings are not in 'dormant grass stage'
            
            if self.old_ht[row,col] > 0.0 and self.old_ht[row,col] < self.StandIndexLP:  # check to see if tree is in a cell (and < max ht)

                h1 = self.StandIndexLP * (1.0 - math.exp(-0.02 * self.old_age[row][col])) # Platt data
                h2 = self.StandIndexLP * (1.0 - math.exp(-0.02 * (self.old_age[row][col] + 1.0)))

                growthLP = (h2 - h1) * (1.0 - self.CI_LP) # CI impacts growth increment

                self.new_ht[row,col] = self.old_ht[row,col] + growthLP
                self.new_age[row,col] = self.old_age[row,col] + 1 # one year older

            elif self.old_ht[row,col] >= self.StandIndexLP:  # max ht 

                self.new_ht[row,col] = self.StandIndexLP
                self.new_age[row,col] = self.old_age[row,col] + 1 # one year older
                 
            else: # if no tree or has germinated seedling < 1m
                pass
        else: # if dorm_tag > 0
            self.new_age[row,col] = self.old_age[row,col] + 1 # one year older even if dormant
            self.new_ht[row,col] = self.germ_size # stays at 0.1 - germination size


    def ring_id(self, row, col, depth=1): 

        ''' recording the neighborhood cells that influence the focal cell,
        assuming a TORUS landscape, a ring (depth) of 1,2,or 3 adjacent cells
        from the focal cell is output '''

        '''INPUTS: depth value (1-3), dim; OUTPUTS: ring cells (self.cells) '''

        self.cells = []
        self.offset = range(-depth, depth + 1)
        self.x = -depth
        for y in self.offset:

            self.xc = self.x + row
            if self.xc < 0: self.xc = self.xc + self.dim
            if self.xc >= self.dim: self.xc = self.xc - self.dim

            self.yc = y + col
            if self.yc < 0: self.yc = self.yc + self.dim
            if self.yc >= self.dim: self.yc = self.yc - self.dim

            self.cells.append([self.xc, self.yc])  # top row of cells

        self.x = depth
        for y in self.offset:

            self.xc = self.x + row
            if self.xc < 0: self.xc = self.xc + self.dim
            if self.xc >= self.dim: self.xc = self.xc - self.dim

            self.yc = y + col
            if self.yc < 0: self.yc = self.yc + self.dim
            if self.yc >= self.dim: self.yc = self.yc - self.dim
            
            self.cells.append([self.xc, self.yc])  # bottom row of cells

        self.y = -depth
        self.xo = self.offset[1:depth * 2]
        for x in self.xo:

            self.xc = x + row
            if self.xc < 0: self.xc = self.xc + self.dim
            if self.xc >= self.dim: self.xc = self.xc - self.dim

            self.yc = self.y + col
            if self.yc < 0: self.yc = self.yc + self.dim
            if self.yc >= self.dim: self.yc = self.yc - self.dim

            self.cells.append([self.xc, self.yc])   # left column, no corners

        self.y = depth
        for x in self.xo:

            self.xc = x + row
            if self.xc < 0: self.xc = self.xc + self.dim
            if self.xc >= self.dim: self.xc = self.xc - self.dim

            self.yc = self.y + col
            if self.yc < 0: self.yc = self.yc + self.dim
            if self.yc >= self.dim: self.yc = self.yc - self.dim

            self.cells.append([self.xc ,self.yc])   # right column, no corners   

        return self.cells
    
    def mort_comp(self, row,col): 
        
        ''' mortality in cell due to intra- and inter-specific competition (mc), a function of CI_LP
        it seems that this is a very sensitive function, adjusted function in curve expert'''

        '''INPUTS: CI_LP,  mortflag, old_ht, rnd_value; OUTPUTS: new_age, new_ht, mortflag '''

        self.IndexLP(row, col)       
        
        dead_cell = 0 # no. of dead trees in cell
        
        if self.mortflag[row][col] == 0:
            if self.old_ht[row][col] > 0.0: # if tree present
                
                if self.old_ht[row][col] < 10.0 and self.CI_LP > 0.5: self.CI_LP = 1.0
                
                mc = -0.1/(1.0 - 45.0 * math.exp(-3.7 * self.CI_LP)) # logistic function
                if self.recruit == 'no': mc = mc * 2 # if past the 'critical time since fire'
                #if self.fireflag[row][col] == 0.0:  mc = mc * self.tsf_mort # this does not work!
                
                for i in range(int(self.old_LPcount[row][col])):# check all trees in cell                 
                    self.rand_calc()
                    if self.rnd_value < mc:
                        self.dead_comp = self.dead_comp + 1
                        dead_cell = dead_cell + 1
                        
                if dead_cell == self.old_LPcount[row][col]:
                    self.new_age[row][col] = 0.0 # all trees die and age & ht is reset to 0
                    self.new_ht[row][col] = 0.0
                    self.new_LPcount[row][col] = 0.0
                    self.mortflag[row][col] = 1.0
                    self.dead_htLP[row][col] = self.old_ht[row][col] # tracking dead LP hts
                    self.dead_htLPc[row][col] = self.old_ht[row][col]
                    
                else: # if not all trees in cell died
                    self.new_LPcount[row][col] = self.old_LPcount[row][col] - dead_cell # removes dead trees
                    self.old_LPcount[row][col] = self.new_LPcount[row][col] # these trees will be checked for other mortality and can grow
                    self.dead_htLP[row][col] = self.old_ht[row][col]
                    self.dead_htLPc[row][col] = self.old_ht[row][col]
                    
                self.dead_ctLPc[row][col] = self.dead_ctLPc[row][col] + dead_cell # recording dead tree count
                self.dead_ctLP[row][col] = self.dead_ctLP[row][col] + dead_cell 
        
    def mort_nat(self, row, col):
        
        '''mortality due to natural causes, simply probability based'''

        ''' adult mortality: a logistic function, with threshold max of ~ 0.03 prob of dying.
        'younger than adult' mortality: a weibull function.  
        From MPM mortality values: Platt EA 1988 & Loudermilk and Cropper 2007.
        Still higher mortality when young, decreases quickly, then increases
        at very large sizes, also see Kaiser dissertation. 'Bathtub mortality curve'   '''

        '''INPUTS: mortflag, old_ht, adult, subadult, rnd_value; OUTPUTS: new_age, new_ht, mortflag'''
        
        dead_cell = 0
        
        if self.mortflag[row][col] == 0:
            if 0.0 < self.old_ht[row][col] <= self.adult:
                young_test = 0.177 - 0.193 * math.exp(-0.21 * self.old_ht[row][col] ** -0.224) # NEW; weibull function
                
                for i in range(int(self.old_LPcount[row][col])): # check all trees in cell
                    self.rand_calc()
                    if self.rnd_value < young_test: 
                        self.dead_nat = self.dead_nat + 1
                        dead_cell = dead_cell + 1
                        
                if dead_cell == self.old_LPcount[row][col]:
                    self.new_age[row][col] = 0.0 # all trees die and age & ht is reset to 0
                    self.new_ht[row][col] = 0.0
                    self.new_LPcount[row][col] = 0.0
                    self.mortflag[row][col] = 1.0
                    self.dead_htLP[row][col] = self.old_ht[row][col] # tracking dead LP hts
                    self.dead_htLPn[row][col] = self.old_ht[row][col]
                    
                else: # if not all trees in cell died
                    self.new_LPcount[row][col] = self.old_LPcount[row][col] - dead_cell # removes dead trees
                    self.old_LPcount[row][col] = self.new_LPcount[row][col] # these trees will be checked for other mortality
                    self.dead_htLP[row][col] = self.old_ht[row][col]
                    self.dead_htLPn[row][col] = self.old_ht[row][col]
                    
                self.dead_ctLPn[row][col] = self.dead_ctLPn[row][col] + dead_cell
                self.dead_ctLP[row][col] = self.dead_ctLP[row][col] + dead_cell

            if self.old_ht[row][col] > self.adult:
                adult_test = 0.03/(1.0 + 593063.0 * math.exp(-0.43 * self.old_ht[row][col])) # NEW; logistic function
                
                for i in range(int(self.old_LPcount[row][col])): # check all trees in cell
                    self.rand_calc()
                    if self.rnd_value < adult_test: 
                        self.dead_nat = self.dead_nat + 1
                        dead_cell = dead_cell + 1
                        
                if dead_cell == self.old_LPcount[row][col]:
                    self.new_age[row][col] = 0.0 # all trees die and age & ht is reset to 0
                    self.new_ht[row][col] = 0.0
                    self.new_LPcount[row][col] = 0.0
                    self.mortflag[row][col] = 1.0
                    self.dead_htLP[row][col] = self.old_ht[row][col]
                    self.dead_htLPn[row][col] = self.old_ht[row][col]
                    
                else: # if not all trees in cell died
                    self.new_LPcount[row][col] = self.old_LPcount[row][col] - dead_cell # removes dead trees, but hts and age remain the same for remaining trees
                    self.old_LPcount[row][col] = self.new_LPcount[row][col] # these trees will be checked for other mortality      
                    self.dead_htLP[row][col] = self.old_ht[row][col]
                    self.dead_htLPn[row][col] = self.old_ht[row][col] 
                    
                self.dead_ctLPn[row][col] = self.dead_ctLPn[row][col] + dead_cell
                self.dead_ctLP[row][col] = self.dead_ctLP[row][col] + dead_cell

    def mort_fire(self, row, col):
        
        ''' mortality due to fire, a function of
        litter accumulation (fire intensity) and tree height,
        use burnt litter (HW litter, LP litter, and wiregrass) for calc fire intensity '''

        '''INPUTS: fireflag, mortflag, old_ht, burntlit, rnd_value; OUTPUTS: new_age, new_ht, mortflag'''
        
        self.rand_calc()
        burntlit = self.total_burntlit # HW, LP, WG
        dead_cell = 0

        if self.fireflag[row][col] == 1: # if fires occurred
            
            if self.mortflag[row][col] == 0: # tree has not been killed yet
                if self.old_ht[row][col] > 0.0: # tree present
                    
                    fire_int = 1.0 - math.exp(-burntlit[row][col]) 
                    size_prob = math.exp(0.2 * (-self.old_ht[row][col] + 1.0))
                    #mf = (fire_int + size_prob)/2 # killed all LP???
                    mf = (fire_int * size_prob) * self.fire_constLP
                    
                    if self.tsf_mort >= self.critical_TSF and mf < self.high_mort: # 
                        
                        mf = self.high_mort # no mort rate < high_mort 
                         
                    else: # follow normal fire mortality probability
                        pass
                                             
                    for i in range(int(self.old_LPcount[row][col])): # check all trees in cell
                        self.rand_calc()
                        if self.rnd_value < mf: # tree dies
                            self.dead_fire = self.dead_fire + 1
                            dead_cell = dead_cell + 1
                            
                    if dead_cell == self.old_LPcount[row][col]:
                        self.new_age[row][col] = 0.0 # all trees die and age & ht is reset to 0
                        self.new_ht[row][col] = 0.0
                        self.new_LPcount[row][col] = 0.0
                        self.mortflag[row][col] = 1.0
                        self.dead_htLP[row][col] = self.old_ht[row][col]
                        self.dead_htLPf[row][col] = self.old_ht[row][col]
                        
                    else: # if not all trees in cell died
                        self.new_LPcount[row][col] = self.old_LPcount[row][col] - dead_cell # removes dead trees, but hts and age remain the same for remaining trees
                        self.old_LPcount[row][col] = self.new_LPcount[row][col] # these trees will be checked for other mortality and can grow
                        self.dead_htLP[row][col] = self.old_ht[row][col]
                        self.dead_htLPf[row][col] = self.old_ht[row][col]

                    self.dead_ctLPf[row][col] = self.dead_ctLPf[row][col] + dead_cell
                    self.dead_ctLP[row][col] = self.dead_ctLP[row][col] + dead_cell
                    

    ''' FIRE section '''

    def fire_ign(self):

        ''' fire ignition probability; assuming that fires start from one random cell  '''

        '''INPUTS: dim, rnd_value; OUTPUTS: ig_cell'''
        
        ignite = zeros((self.dim, self.dim))
        for row in range(self.dim):
            for col in range(self.dim):
                self.rand_calc()
                if self.rnd_value < 0.25: # 25% chance of being an ignition cell (arbitrary, not in fortran)
                    ignite[row][col] = 1.0 # fire COULD ignite in this cell
        
        ign_pos = []
        for i in range(len(ignite)):
            for j in range(len(ignite)):
                if ignite[i][j] == 1.0:
                    ign_pos.append([i,j])

        if self.randfromfile:
           a = int(self.randintnumberfromfile[self.rcount_int])
           self.rcount_int+=1
        else:
           a = random.randint(0,len(ignite))
           self.randint_value.append(a)

        # added by EJ      
        self.ig_cell = ign_pos[a] # randomly chosen ignition cell
       
    def fire_spread_flag(self):
        
        ''' fire spread module, based simply on area desired to burn '''

        '''INPUTS: fireflag, dim, ring cells; OUTPUTS: fireflag'''
        
        a = self.ig_cell
        self.fireflag[a[0],a[1]] = 1.0
        burndist = int(self.dim/2 + 1)  #random.randint(1,self.dim/2+1) 
        bd_range = range(1,burndist+1)
        self.check_cells = []
        num_burned_cells = 0

        if burndist > (self.dim/2): # all cells are checked for fire spread
            for row in range(self.dim):
                for col in range(self.dim):                    
                    self.fireflag[row,col] = 1.0 # cell burns
                    num_burned_cells = num_burned_cells + 1
                        
        else: # check for particular area (based on burn distance) to be burned
            for i in bd_range:
                b = self.ring_id(a[0],a[1],i)
                for j in range(len(b)):
                    self.check_cells.append(b[j]) # save all cells in each ring around focal cell

            for k in self.check_cells:
                self.fireflag[k[0],k[1]] = 1.0
                    
    def fire_spread_burn(self):

        ''' burns the litter where the fire has spread '''

        ''' INPUTS: fireflag, litter, litterHW; OUTPUTS: total_burntlit '''
        
        nburned = 0

        for row in range(self.dim):
            for col in range(self.dim):
                if self.fireflag[row,col] == 1.0:
                    nburned = nburned + 1
                    self.burntlit[row,col]= self.litter[row,col] # save LLP litter that was burnt
                    self.burntlitHW[row,col]= self.litterHW[row,col] # HW
                    self.burntlitWG[row,col]= self.litterWG[row,col] # wiregrass
                    self.total_burntlit[row,col]= self.burntlitHW[row,col] + self.burntlit[row,col] +  self.burntlitWG[row,col]
                    # total litter: HW, LLP, WG, first thing run when a fire occurs, used for mort_fire, both HW, LP

        self.burn_area.append(nburned)

    def size_class_A(self,row,col):
     
        ''' Split the alive and dead LLP into size classes, based on MPM size classes
            HW size classes are arbitrary '''

        ''' inputs: new_LPcount, new_ht, dead_htLP, dead_ctLP; outputs: class_LPcountA, class_LPcountD '''

        ''' ALIVE size classes'''
    
        if self.old_ht[row][col] > 0.0: # LP tree present

            if self.old_ht[row][col] >= 1.5 and self.old_ht[row][col] < 1.55: # not necessary, monitoring ages of 1.5m LP
                self.hts.append(self.old_ht[row][col])
                self.ages.append(self.old_age[row][col])
           
            if self.old_ht[row][col] < 1.0:
                self.class_LPcountsA[0] = self.class_LPcountsA[0] + self.old_LPcount[row][col]
            elif 1.0 <= self.old_ht[row][col] < 3.0:
                self.class_LPcountsA[1] = self.class_LPcountsA[1] + self.old_LPcount[row][col]
            elif 3.0 <= self.old_ht[row][col] < 11.0:
                self.class_LPcountsA[2] = self.class_LPcountsA[2] + self.old_LPcount[row][col]
            elif 11.0 <= self.old_ht[row][col] < 18.0:
                self.class_LPcountsA[3] = self.class_LPcountsA[3] + self.old_LPcount[row][col]
            elif 18.0 <= self.old_ht[row][col] < 23.0:
                self.class_LPcountsA[4] = self.class_LPcountsA[3] + self.old_LPcount[row][col]                       
            elif 23.0 <= self.old_ht[row][col] < 27.0:
                self.class_LPcountsA[5] = self.class_LPcountsA[4] + self.old_LPcount[row][col]                  
            elif 27.0 <= self.old_ht[row][col] < 30.0:
                self.class_LPcountsA[6] = self.class_LPcountsA[5] + self.old_LPcount[row][col]  
            elif 30.0 <= self.old_ht[row][col] < 32.0:
                self.class_LPcountsA[7] = self.class_LPcountsA[7] + self.old_LPcount[row][col]
            elif 32.0 <= self.old_ht[row][col] < 34.0:
                self.class_LPcountsA[8] = self.class_LPcountsA[8] + self.old_LPcount[row][col]  
            else: # >= 34 m ht
                self.class_LPcountsA[9] = self.class_LPcountsA[9] + self.old_LPcount[row][col]
                
        if self.old_htHW[row][col] > 0.0: # HW tree present
            
            if self.old_htHW[row][col] <= 1.5: # those from 'linear growth'
                self.class_HWcountsA[0] = self.class_HWcountsA[0] + self.old_HWcount[row][col]
            elif 1.5 < self.old_htHW[row][col] <= 3.0:
                self.class_HWcountsA[1] = self.class_HWcountsA[1] + self.old_HWcount[row][col]
            elif 3.0 < self.old_htHW[row][col] <= 5.0:
                self.class_HWcountsA[2] = self.class_HWcountsA[2] + self.old_HWcount[row][col]                       
            elif 5.0 < self.old_htHW[row][col] <= 10.0:
                self.class_HWcountsA[3] = self.class_HWcountsA[3] + self.old_HWcount[row][col]                  
            elif 10.0 < self.old_htHW[row][col] <= 15.0:
                self.class_HWcountsA[4] = self.class_HWcountsA[4] + self.old_HWcount[row][col]  
            elif 15.0 < self.old_htHW[row][col] <= 20.0:
                self.class_HWcountsA[5] = self.class_HWcountsA[5] + self.old_HWcount[row][col] 
            else: # >= 20 m ht
                self.class_HWcountsA[6] = self.class_HWcountsA[6] + self.old_HWcount[row][col]

    def size_class_D(self,row,col):
        
        ''' DEAD LP size classes'''
    
        if self.dead_htLP[row][col] > 0.0: # dead LP tree present, all 
            
            if self.dead_htLP[row][col] < 1.0:
                self.class_LPcountsD[0] = self.class_LPcountsD[0] + self.dead_ctLP[row][col]
            elif 1.0 <= self.dead_htLP[row][col] < 3.0:
                self.class_LPcountsD[1] = self.class_LPcountsD[1] + self.dead_ctLP[row][col]
            elif 3.0 <= self.dead_htLP[row][col] < 11.0:
                self.class_LPcountsD[2] = self.class_LPcountsD[2] + self.dead_ctLP[row][col]                
            elif 11.0 <= self.dead_htLP[row][col] < 18.0:
                self.class_LPcountsD[3] = self.class_LPcountsD[3] + self.dead_ctLP[row][col]
            elif 18.0 <= self.dead_htLP[row][col] < 23.0:
                self.class_LPcountsD[4] = self.class_LPcountsD[4] + self.dead_ctLP[row][col]                       
            elif 23.0 <= self.dead_htLP[row][col] < 27.0:
                self.class_LPcountsD[5] = self.class_LPcountsD[5] + self.dead_ctLP[row][col]                  
            elif 27.0 <= self.dead_htLP[row][col] < 30.0:
                self.class_LPcountsD[6] = self.class_LPcountsD[6] + self.dead_ctLP[row][col]  
            elif 30.0 <= self.dead_htLP[row][col] < 32.0:
                self.class_LPcountsD[7] = self.class_LPcountsD[7] + self.dead_ctLP[row][col]
            elif 32.0 <= self.dead_htLP[row][col] < 34.0:
                self.class_LPcountsD[8] = self.class_LPcountsD[8] + self.dead_ctLP[row][col]  
            else: # >= 34 m ht
                self.class_LPcountsD[9] = self.class_LPcountsD[9] + self.dead_ctLP[row][col]
                
        if self.dead_htLPc[row][col] > 0.0: # dead tree present, from competition
            
            if self.dead_htLPc[row][col] < 1.0:
                self.class_LPcountsDc[0] = self.class_LPcountsDc[0] + self.dead_ctLPc[row][col]
            elif 1.0 <= self.dead_htLPc[row][col] < 3.0:
                self.class_LPcountsDc[1] = self.class_LPcountsDc[1] + self.dead_ctLPc[row][col]
            elif 3.0 <= self.dead_htLPc[row][col] < 11.0:
                self.class_LPcountsDc[2] = self.class_LPcountsDc[2] + self.dead_ctLPc[row][col]                
            elif 11.0 <= self.dead_htLPc[row][col] < 18.0:
                self.class_LPcountsDc[3] = self.class_LPcountsDc[3] + self.dead_ctLPc[row][col]
            elif 18.0 <= self.dead_htLPc[row][col] < 23.0:
                self.class_LPcountsDc[4] = self.class_LPcountsDc[4] + self.dead_ctLPc[row][col]                       
            elif 23.0 <= self.dead_htLPc[row][col] < 27.0:
                self.class_LPcountsDc[5] = self.class_LPcountsDc[5] + self.dead_ctLPc[row][col]                  
            elif 27.0 <= self.dead_htLPc[row][col] < 30.0:
                self.class_LPcountsDc[6] = self.class_LPcountsDc[6] + self.dead_ctLPc[row][col]  
            elif 30.0 <= self.dead_htLPc[row][col] < 32.0:
                self.class_LPcountsDc[7] = self.class_LPcountsDc[7] + self.dead_ctLPc[row][col]
            elif 32.0 <= self.dead_htLPc[row][col] < 34.0:
                self.class_LPcountsDc[8] = self.class_LPcountsDc[8] + self.dead_ctLPc[row][col]  
            else: # >= 34 m ht
                self.class_LPcountsDc[9] = self.class_LPcountsDc[9] + self.dead_ctLPc[row][col]
                
        if self.dead_htLPn[row][col] > 0.0: # dead tree present, natural mortality
            
            if self.dead_htLPn[row][col] < 1.0:
                self.class_LPcountsDn[0] = self.class_LPcountsDn[0] + self.dead_ctLPn[row][col]
            elif 1.0 <= self.dead_htLPn[row][col] < 3.0:
                self.class_LPcountsDn[1] = self.class_LPcountsDn[1] + self.dead_ctLPn[row][col]
            elif 3.0 <= self.dead_htLPn[row][col] < 11.0:
                self.class_LPcountsDn[2] = self.class_LPcountsDn[2] + self.dead_ctLPn[row][col]                
            elif 11.0 <= self.dead_htLPn[row][col] < 18.0:
                self.class_LPcountsDn[3] = self.class_LPcountsDn[3] + self.dead_ctLPn[row][col]
            elif 18.0 <= self.dead_htLPn[row][col] < 23.0:
                self.class_LPcountsDn[4] = self.class_LPcountsDn[4] + self.dead_ctLPn[row][col]                       
            elif 23.0 <= self.dead_htLPn[row][col] < 27.0:
                self.class_LPcountsDn[5] = self.class_LPcountsDn[5] + self.dead_ctLPn[row][col]                  
            elif 27.0 <= self.dead_htLPn[row][col] < 30.0:
                self.class_LPcountsDn[6] = self.class_LPcountsDn[6] + self.dead_ctLPn[row][col]  
            elif 30.0 <= self.dead_htLPn[row][col] < 32.0:
                self.class_LPcountsDn[7] = self.class_LPcountsDn[7] + self.dead_ctLPn[row][col]
            elif 32.0 <= self.dead_htLPn[row][col] < 34.0:
                self.class_LPcountsDn[8] = self.class_LPcountsDn[8] + self.dead_ctLPn[row][col]  
            else: # >= 34 m ht
                self.class_LPcountsDn[9] = self.class_LPcountsDn[9] + self.dead_ctLPn[row][col]

        if self.dead_htLPf[row][col] > 0.0: # dead tree present, from fire

            if self.dead_htLPf[row][col] < 1.0:
                self.class_LPcountsDf[0] = self.class_LPcountsDf[0] + self.dead_ctLPf[row][col]
            elif 1.0 <= self.dead_htLPf[row][col] < 3.0:
                self.class_LPcountsDf[1] = self.class_LPcountsDf[1] + self.dead_ctLPf[row][col]
            elif 3.0 <= self.dead_htLPf[row][col] < 11.0:
                self.class_LPcountsDf[2] = self.class_LPcountsDf[2] + self.dead_ctLPf[row][col]                
            elif 11.0 <= self.dead_htLPf[row][col] < 18.0:
                self.class_LPcountsDf[3] = self.class_LPcountsDf[3] + self.dead_ctLPf[row][col]
            elif 18.0 <= self.dead_htLPf[row][col] < 23.0:
                self.class_LPcountsDf[4] = self.class_LPcountsDf[4] + self.dead_ctLPf[row][col]                       
            elif 23.0 <= self.dead_htLPf[row][col] < 27.0:
                self.class_LPcountsDf[5] = self.class_LPcountsDf[5] + self.dead_ctLPf[row][col]                  
            elif 27.0 <= self.dead_htLPf[row][col] < 30.0:
                self.class_LPcountsDf[6] = self.class_LPcountsDf[6] + self.dead_ctLPf[row][col]  
            elif 30.0 <= self.dead_htLPf[row][col] < 32.0:
                self.class_LPcountsDf[7] = self.class_LPcountsDf[7] + self.dead_ctLPf[row][col]
            elif 32.0 <= self.dead_htLPf[row][col] < 34.0:
                self.class_LPcountsDf[8] = self.class_LPcountsDf[8] + self.dead_ctLPf[row][col]  
            else: # >= 34 m ht
                self.class_LPcountsDf[9] = self.class_LPcountsDf[9] + self.dead_ctLPf[row][col]

        ''' Dead HW size classes '''

        if self.dead_htHW[row][col] > 0.0: # dead HW tree present, all 
            
            if self.dead_htHW[row][col] <= 1.5:
                self.class_HWcountsD[0] = self.class_HWcountsD[0] + self.dead_ctHW[row][col]                
            elif 1.5 < self.dead_htHW[row][col] <= 3.0:
                self.class_HWcountsD[1] = self.class_HWcountsD[1] + self.dead_ctHW[row][col]               
            elif 3.0 < self.dead_htHW[row][col] <= 5.0:
                self.class_HWcountsD[2] = self.class_HWcountsD[2] + self.dead_ctHW[row][col]               
            elif 5.0 < self.dead_htHW[row][col] <= 10.0:
                self.class_HWcountsD[3] = self.class_HWcountsD[3] + self.dead_ctHW[row][col]              
            elif 10.0 < self.dead_htHW[row][col] <= 15.0:
                self.class_HWcountsD[4] = self.class_HWcountsD[4] + self.dead_ctHW[row][col]            
            elif 15.0 < self.dead_htHW[row][col] <= 20.0:
                self.class_HWcountsD[5] = self.class_HWcountsD[5] + self.dead_ctHW[row][col]           
            else: # > 20 m ht
                self.class_HWcountsD[6] = self.class_HWcountsD[6] + self.dead_ctHW[row][col]
                
        if self.dead_htHWn[row][col] > 0.0: # dead HW tree present, natural mortality 
            
            if self.dead_htHWn[row][col] <= 1.5:
                self.class_HWcountsDn[0] = self.class_HWcountsDn[0] + self.dead_ctHWn[row][col]       
            elif 1.5 < self.dead_htHWn[row][col] <= 3.0:
                self.class_HWcountsDn[1] = self.class_HWcountsDn[1] + self.dead_ctHWn[row][col]          
            elif 3.0 < self.dead_htHWn[row][col] <= 5.0:
                self.class_HWcountsDn[2] = self.class_HWcountsDn[2] + self.dead_ctHWn[row][col]       
            elif 5.0 < self.dead_htHW[row][col] <= 10.0:
                self.class_HWcountsDn[3] = self.class_HWcountsDn[3] + self.dead_ctHWn[row][col]                
            elif 10.0 < self.dead_htHWn[row][col] <= 15.0:
                self.class_HWcountsDn[4] = self.class_HWcountsDn[4] + self.dead_ctHWn[row][col]              
            elif 15.0 < self.dead_htHWn[row][col] <= 20.0:
                self.class_HWcountsDn[5] = self.class_HWcountsDn[5] + self.dead_ctHWn[row][col]            
            else: # > 20 m ht
                self.class_HWcountsDn[6] = self.class_HWcountsDn[6] + self.dead_ctHWn[row][col]

        if self.dead_htHWf[row][col] > 0.0: # dead HW tree present, fire 
            
            if self.dead_htHWf[row][col] <= 1.5:
                self.class_HWcountsDf[0] = self.class_HWcountsDf[0] + self.dead_ctHWf[row][col]               
            elif 1.5 < self.dead_htHWf[row][col] <= 3.0:
                self.class_HWcountsDf[1] = self.class_HWcountsDf[1] + self.dead_ctHWf[row][col]               
            elif 3.0 < self.dead_htHWf[row][col] <= 5.0:
                self.class_HWcountsDf[2] = self.class_HWcountsDf[2] + self.dead_ctHWf[row][col]            
            elif 5.0 < self.dead_htHWf[row][col] <= 10.0:
                self.class_HWcountsDf[3] = self.class_HWcountsDf[3] + self.dead_ctHWf[row][col]           
            elif 10.0 < self.dead_htHWf[row][col] <= 15.0:
                self.class_HWcountsDf[4] = self.class_HWcountsDf[4] + self.dead_ctHWf[row][col]               
            elif 15.0 < self.dead_htHWf[row][col] <= 20.0:
                self.class_HWcountsDf[5] = self.class_HWcountsDf[5] + self.dead_ctHWf[row][col]                
            else: # > 20 m ht
                self.class_HWcountsDf[6] = self.class_HWcountsDf[6] + self.dead_ctHWf[row][col]
                
    def size_mort_rates(self):

        ''' Function to check mortality rates in size classes for comparison with MPM values;
        Tracking values for all simulation years; Should be performed after size_classLP()'''

        mort  = [0] * self.size_LP
        mortc = [0] * self.size_LP
        mortn = [0] * self.size_LP
        mortf = [0] * self.size_LP
        
        morth  = [0] * self.size_HW
        mortnh = [0] * self.size_HW
        mortfh = [0] * self.size_HW

        for i in range(len(mort)):
            
            if self.class_LPcountsA[i] > 0.0:

                mort[i] = self.class_LPcountsD[i]/self.class_LPcountsA[i] 
                mortc[i] = self.class_LPcountsDc[i]/self.class_LPcountsA[i] 
                mortn[i] = self.class_LPcountsDn[i]/self.class_LPcountsA[i] 
                mortf[i] = self.class_LPcountsDf[i]/self.class_LPcountsA[i]

            else: #elif self.class_LPcountsA[i] == 0.0:
                mort[i], mortc[i],mortn[i],mortf[i] = 0.0, 0.0, 0.0, 0.0

        for i in range(len(morth)):
            
            if self.class_HWcountsA[i] > 0.0:

                morth[i] = self.class_HWcountsD[i]/self.class_HWcountsA[i] 
                mortnh[i] = self.class_HWcountsDn[i]/self.class_HWcountsA[i] 
                mortfh[i] = self.class_HWcountsDf[i]/self.class_HWcountsA[i]

            else: #elif self.class_HWcountsA[i] == 0.0:
                
                morth[i], mortnh[i],mortfh[i] = 0.0, 0.0, 0.0
                
        self.tracking_mort_rates.append(mort)
        self.tracking_mortc_rates.append(mortc)
        self.tracking_mortn_rates.append(mortn)
        self.tracking_mortf_rates.append(mortf)

        self.tracking_morth_rates.append(morth)
        self.tracking_mortnh_rates.append(mortnh)
        self.tracking_mortfh_rates.append(mortfh) 
    
    def run(self,nyears): #  where all functions above are run, Note: CONTINUOUSLY check order of functions!
        
        #EJ introduced this parameters to track the mast and fire years
        self.fyears=np.zeros(nyears+1) 
        self.myears=np.zeros(nyears+1)
        #JE
        #self.beg_stand()
        
        if self.graph_yn == 'on':
            ion()      
            self.display_stand()
            if self.D.displayed != 'f_fire' and self.D.displayed != 'tree_histLP'\
               and self.D.displayed != 'tree_histHW' and self.D.displayed != 'dtree_histLP'\
               and self.D.displayed != 'dtree_histHW':
                colorbar()
                

        elif self.graph_yn_m == 'on':
            ion()
            
            n = len(self.D.displayed_multi)           
            self.multi_display(n, 0)

        else: pass

        max_litter = 0.0
        max_sum_ltter = 0.0
        max_htsum = 0.0
        max_seeds = 0.0
        
        #EJ modified
        #readfireprobfromfile=1
        #readmastprobfromfile=1
        if self.readfireprobfromfile:
           print ('Read from fire file .....')
           yfire = np.loadtxt('LLM_fire_years.out', unpack=True)
        if self.readmastprobfromfile:
           print ('Read from must file .....')
           ymast = np.loadtxt('LLM_mast_years.out', unpack=True)
        #JE

        fir = 0
        for tim in range(1, nyears+1):
           if self.verbose:
              print ('year ',tim)

           ntree = 0.0 # tree count, adults
           ntree_all = 0.0 # tree counts - all
           ntreeHW = 0.0
           ntreeHW_all = 0.0
           dtree = 0.0 # tree density
           dtreeHW = 0.0
           nt8 = 0.0 # tree count, young
           nt8HW = 0.0    
            
##            fir = fir + 1
##            if fir > 49: # fire at specific time interval
            #EJ modifed
           if self.readfireprobfromfile:
              if yfire[tim-1] ==1:
                 if self.verbose:
                    print ('fire occurs in year ', tim)
                 fire = 'yes'
                 self.fyears[tim]=1
                 self.tracking_burnyrs.append(1)
                 self.fire_ign()
                 if (self.fire_prob>0):
                     self.fire_spread_flag()
                 fir = 0
                 self.tsf_mort = self.time_since_fire[0][0] + 1 # holding time since fire for mort_fire()
              else:
                 fire = 'no'
                 self.tracking_burnyrs.append(0)
                 self.burn_area.append(0)
           else: 
            #JE           
              self.rand_calc()
              if self.rnd_value < self.fire_prob or tim == 1: # prob of fires occuring
              # fire always occurs the 1st year - easier to do time since fire & tsf_mort . . .
                 if self.verbose:
                    print ('fire occurs in year ', tim)
                 fire = 'yes'
                 self.fyears[tim]=1
                 self.tracking_burnyrs.append(1)
                 self.fire_ign()
                 if (self.fire_prob>0):
                    self.fire_spread_flag()
                 fir = 0
                 self.tsf_mort = self.time_since_fire[0][0] + 1 # holding time since fire for mort_fire()
   
              else:
                 fire = 'no'
                 self.tracking_burnyrs.append(0)
                 self.burn_area.append(0)
 
           for row in range(self.dim):
                for col in range(self.dim):
                    if self.dorm_tag[row][col] > 0.0:
                        self.dorm_tag[row][col] = self.dorm_tag[row][col] - 1
                    if self.fireflag[row][col] == 1.0:
                        self.time_since_fire[row][col] = 0.0 # reset time since fire to zero
                    else: # no fire, adding to time since fire
                        self.time_since_fire[row][col] = self.time_since_fire[row][col] + 1

            # no recruitment after tsf reaches 20 (or other critical_TSF)
           if self.time_since_fire[0,0] == self.critical_TSF:
                self.recruit = 'no'
                if self.verbose:
                   print ('no more recruitment! ')
            #EJ modified
           if self.readmastprobfromfile:
              if ymast[tim-1] ==1: 
                 self.myears[tim]=1 
                 self.mast_yr = 'yes'# if mast year
                 self.tracking_masting.append(1)
                 if self.verbose:
                    print ('mast year')
              else: self.tracking_masting.append(0)
           else:
            #JE                        
              self.rand_calc()        
              if self.rnd_value < self.mast_prob: # based on Boyer's data
                 self.myears[tim]=1 
                 self.mast_yr = 'yes'# if mast year
                 self.tracking_masting.append(1)
                 if self.verbose:
                    print ('mast year')
              else: self.tracking_masting.append(0)
                          
              self.seed_prodLP()            
            
              for row in range(self.dim):
                for col in range(self.dim):
                    self.size_class_A(row,col) # get beg. size classes for alive trees
                    self.seed_dispLP(row,col)
                    self.germLP()
                    self.reproHW(row,col)
                    self.HtSumLP(row, col)
                    self.HtSumHW(row,col)

              a = sum(self.tracking_HWsprouts)
              self.total_HWsprouts.append(a)
                    
              tracking = 0
              for row in range(self.dim):
                for col in range(self.dim):
                    if self.new_age[row][col] == 1.0: # seedling established in cell
                        tracking = tracking + 1 # + self.new_LPcount[row][col]
              self.tracking_seedlingsLP.append(tracking) # tracking no. of cells/yr with germ seedlings
                       
              if fire == 'yes': # if fires occurred this year
                 self.fire_spread_burn() # burns litter

                 for i in range(len(self.litter)):
                    for j in range(len(self.litter)):
                        self.litter[i][j] = self.litter[i][j] - self.burntlit[i][j] # = 0
                        self.litterHW[i][j] = self.litterHW[i][j] - self.burntlitHW[i][j] # = 0
                        self.litterWG[i][j] = self.litterWG[i][j] - self.burntlitWG[i][j] # = 0
                        # note: total_burntlit is used for fire intensity, adds all burntlit together
                    
              for row in range(self.dim):
                for col in range(self.dim):
                    self.litLP(row,col)
                    self.litHW(row,col) # HW
                    self.litWG(row,col) # wiregrass                
                    self.mort_comp(row,col) # order is important here for mortflag and dead_ctLP
                    self.mort_nat(row,col)
                    self.mort_natHW(row,col)
                    self.mort_fire(row,col) # uses burntlit
                    self.mort_fireHW(row,col) # HW
                    self.tot_biomass[row][col] = self.litter[row][col] + self.litterHW[row][col] + self.litterWG[row][col]
                    # tot_biomass (after fire) used for graphing ONLY, not fire intensity
                    
              if fire == 'yes':  
                try:                   
##                    print '# survived', self.surv1
##                    print '# alive', self.alive1
##                    print '# dead ', self.dead1
                    sa = self.surv1/self.alive1 # survival rate for time step, <2m
                    da = self.dead1/self.alive1
##                    print '  survival rate ', sa
##                    print '  mortality rate ', da
                    self.surv_smHW_fire.append(sa)
                    self.mort_smHW_fire.append(da)
                except: pass # if no HW< 2m

##            print 'dead from comp ', self.dead_comp
##            print 'dead from nat ', self.dead_nat
##            print 'dead from fire ', self.dead_fire
##            print ''
##            print 'HW dead from nat ', self.dead_natHW
##            print 'HW dead from fire ', self.dead_fireHW
##            print 'HW top-killed from fire ', self.TK_HW

              ''' Counting Trees for Graphs Section '''
            
              tot_HWdead = self.dead_natHW + self.dead_fireHW # total dead for year
              tot_LPdead = self.dead_nat + self.dead_fire + self.dead_comp

              self.HWdead.append(tot_HWdead)
              self.LPdead.append(tot_LPdead)
              self.tracking_deadnatHW.append(self.dead_natHW)
              self.tracking_deadfireHW.append(self.dead_fireHW)
              self.tracking_deadnatLP.append(self.dead_nat)
              self.tracking_deadfireLP.append(self.dead_fire)
              self.tracking_deadcompLP.append(self.dead_comp)
            
              for row in range(self.dim): # number alive
                for col in range(self.dim):
                    
                    self.LPgrowth(row, col)
                    self.HWgrowth(row,col)
                    self.size_class_D(row,col) # get resulting size classes for dead trees

                    ntree_all = ntree_all + self.new_LPcount[row][col]
                    ntreeHW_all = ntreeHW_all + self.new_HWcount[row][col]
                    
                    if self.new_ht[row,col] > self.graph_lg: # = subadult size, there is also a 'graph_sm = 5'

                        ntree = ntree + self.new_LPcount[row][col] # count no. trees > particular ht

                    else:

                        nt8 = nt8 + self.new_LPcount[row][col]

                    if self.new_htHW[row][col] > self.graph_lgHW:
                        
                        ntreeHW = ntreeHW + self.new_HWcount[row][col] # count HWs

                    else:

                        nt8HW = nt8HW + self.new_HWcount[row][col]
##            print ''  
##            print 'LP size class dead ', self.class_LPcountsD
##            print 'LP size class alive', self.class_LPcountsA
##            print ''
##            print ''  
##            print 'HW size class dead ', self.class_HWcountsD
##            print 'HW size class alive', self.class_HWcountsA
##            print ''
            
              self.size_mort_rates()
                     
              self.te = tim

              # tree no. and density
              self.tree_number.append(ntree)
              self.HWtree_number.append(ntreeHW)
              dtree = ntree / (25 * self.dim**2/10000.0) # calc. density (of larger trees) across given size landscape; 25 = size (sq. m) of cell
              dtreeHW = ntreeHW / (25 * self.dim**2/10000.0) # 25 sq. m/cell
              self.tree_density.append(dtree)
              self.HWtree_density.append(dtreeHW)
              self.tlt8.append(nt8)
              self.tlt8HW.append(nt8HW)

              # relative abundance graphs
              rel_ab1 = ntree_all/(ntree_all + ntreeHW_all)
              rel_ab2 = ntreeHW_all/(ntree_all + ntreeHW_all)
              self.relA_LP_all.append(rel_ab1)
              self.relA_HW_all.append(rel_ab2)            

              self.D.c = self.D.c + 1 # for colorbar in multi graph
              self.graph_set()
              if self.graph_yn == 'on':            
                self.display_stand(tim)
              elif self.graph_yn_m == 'on':
                self.multi_display(yr=tim)
              else: pass
            

              ''' Max Values '''
            
              mxline = []
              for i in range(len(self.old_ht)):
                mx = max(self.old_ht[i])
                mxline.append(mx)
              self.max_yrht.append(max(mxline))

              mxline2 = []
              for i in range(len(self.old_htHW)):
                mx = max(self.old_htHW[i])
                mxline2.append(mx)
              self.max_yrhtHW.append(max(mxline2))

              ''' copy data for next time step '''

              self.old_ht = copy.deepcopy(self.new_ht)
              self.old_htHW = copy.deepcopy(self.new_htHW) # HW
              self.old_age = copy.deepcopy(self.new_age)
              self.old_ageHW = copy.deepcopy(self.new_ageHW) # HW
              self.old_LPcount = copy.deepcopy(self.new_LPcount)
              self.old_HWcount = copy.deepcopy(self.new_HWcount) # HW
              self.old_litter = copy.deepcopy(self.litter) # for final map
              self.old_litterHW = copy.deepcopy(self.litterHW) # for final map
              self.old_litterWG = copy.deepcopy(self.litterWG) # for final map
              self.old_tot_biomass = copy.deepcopy(self.tot_biomass) # for final map
              self.old_seeds = copy.deepcopy(self.seeds)
              self.old_sum_htLP1 = copy.deepcopy(self.sum_htLP1)
              self.old_sum_htLP2 = copy.deepcopy(self.sum_htLP2)
              self.old_sum_htHW1 = copy.deepcopy(self.sum_htHW1)
              self.old_sum_htHW2 = copy.deepcopy(self.sum_htHW2)


              if tim<self.noHWyear:
                print (tim,'inside the noHW')
                self.old_HWcount = zeros((self.dim,self.dim), 'float')
                self.old_ageHW = zeros((self.dim,self.dim), 'float')
                self.old_htHW = zeros((self.dim,self.dim), 'float')
                self.litterHW = zeros((self.dim,self.dim), 'float')
                self.old_litterHW = zeros((self.dim,self.dim), 'float') # for final map

              if tim==self.noHWyear: # NOTE, works only for noHWyear>0 
                #start HW stands
                for row in range(self.dim):
                  for col in range(self.dim):
                    self.rand_calc()
                    if self.rnd_value < self.area_of_HW_resprouting: # random landscape of hardwoods, within ~10% of the area (cells)
                      self.old_htHW[row][col] = 0.5
                      self.old_ageHW[row][col] = 1.0 # difficult to determine age of young HW with frequent fires and resprouting
                      self.old_HWcount[row][col] = 1.0

              ''' reset values for next time step '''
            
              self.seeds = zeros((self.dim, self.dim),'float')
              self.germ_seeds = zeros((self.dim, self.dim),'float')
              self.Pg = zeros((self.dim, self.dim),'float')
              self.new_ht = zeros((self.dim, self.dim),'float')
              self.new_htHW = zeros((self.dim, self.dim),'float') # HW
              self.sum_htLP1 = zeros((self.dim,self.dim), 'float')
              self.sum_htLP2 = zeros((self.dim,self.dim), 'float')
              self.sum_htHW1 = zeros((self.dim,self.dim), 'float')
              self.sum_htHW2 = zeros((self.dim,self.dim), 'float')
              self.new_age = zeros((self.dim, self.dim),'float')
              self.new_ageHW = zeros((self.dim, self.dim),'float') # HW
              self.new_LPcount = zeros((self.dim, self.dim),'float') 
              self.new_HWcount = zeros((self.dim, self.dim),'float') # HW
              self.ef = zeros((self.dim,self.dim), 'float')
              self.burntlit = zeros((self.dim, self.dim),'float')
              self.burntlitHW = zeros((self.dim, self.dim),'float') # HW
              self.burntlitWG = zeros((self.dim, self.dim),'float')
              self.total_burntlit = zeros((self.dim,self.dim), 'float') # burnt litter in cell from HW, LLP, WG
              self.mortflag = zeros((self.dim,self.dim))
              self.mortflagHW = zeros((self.dim,self.dim)) # HW
              self.fireflag = zeros((self.dim,self.dim), 'float')
              self.dead_ctHW = zeros((self.dim,self.dim), 'float')# HW
              self.dead_ctHWn = zeros((self.dim,self.dim), 'float')
              self.dead_ctHWf = zeros((self.dim,self.dim), 'float')
              self.dead_ctLP = zeros((self.dim,self.dim), 'float')
              self.dead_ctLPc = zeros((self.dim,self.dim), 'float')
              self.dead_ctLPn = zeros((self.dim,self.dim), 'float')
              self.dead_ctLPf = zeros((self.dim,self.dim), 'float')
              self.dead_htHW = zeros((self.dim,self.dim), 'float')# HW
              self.dead_htHWn = zeros((self.dim,self.dim), 'float')
              self.dead_htHWf = zeros((self.dim,self.dim), 'float')
              self.dead_htLP = zeros((self.dim,self.dim), 'float')
              self.dead_htLPc = zeros((self.dim,self.dim), 'float')
              self.dead_htLPn = zeros((self.dim,self.dim), 'float')
              self.dead_htLPf = zeros((self.dim,self.dim), 'float')
              self.mast_yr = 'no'
              self.copius_yr = 'no'
              self.dead_comp = 0
              self.dead_nat = 0
              self.dead_natHW = 0 # HW
              self.dead_fire = 0
              self.dead_fireHW = 0 # HW
              self.TK_HW = 0 # HW
              self.tsf_mort = 0
              self.disperse = []

              self.class_LPcountsA = [0] * self.size_LP
              self.class_LPcountsD = [0] * self.size_LP
              self.class_LPcountsDc = [0] * self.size_LP
              self.class_LPcountsDn = [0] * self.size_LP
              self.class_LPcountsDf = [0] * self.size_LP
            
              self.class_HWcountsA = [0] * self.size_HW
              self.class_HWcountsD = [0] * self.size_HW
              self.class_HWcountsDn = [0] * self.size_HW
              self.class_HWcountsDf = [0] * self.size_HW
            
              #EJ modified
              hsi=self.HSI
              hsi.lb = self.lb[0]
              hsi.ub = self.ub[0]
              hsi.bench = self.bn[0]
              hsi.smooth = self.sm[0]
              hsi.nsize = int(self.xsize[0])
              hsi.msize = int(self.ysize[0])
              self.age_sc.append(hsi.score_subset(self.old_age))
              hsi.lb = self.lb[1]
              hsi.ub = self.ub[1]
              hsi.bench = self.bn[1]
              hsi.smooth = self.sm[1]
              hsi.nsize = int(self.xsize[1])
              hsi.msize = int(self.ysize[1])
              self.hw_sc.append(hsi.score_subset(self.old_ht)) 
              hsi.lb = self.lb[3]
              hsi.ub = self.ub[3]
              hsi.bench = self.bn[3]
              hsi.smooth = self.sm[3]
              hsi.nsize = int(self.xsize[3])
              hsi.msize = int(self.ysize[3])
              self.ageHW_sc.append(hsi.score_subset(self.old_ageHW))
              hsi.lb = self.lb[4]
              hsi.ub = self.ub[4]
              hsi.bench = self.bn[4]
              hsi.smooth = self.sm[4]
              hsi.nsize = int(self.xsize[4])
              hsi.msize = int(self.ysize[4])
              self.hwHW_sc.append(hsi.score_subset(self.old_htHW))
        
              hsi.nsize = 5
              hsi.msize = 5
        
              self.sq_sc.append(hsi.prob_of_use_sq(self.old_age, self.old_LPcount, self.old_ageHW, self.old_HWcount, self.tree_mature_age))

              #self.gt_sc.append(hsi.prob_of_use_gt(self.old_ht, self.old_LPcount, self.old_htHW, self.old_HWcount))
              self.gt_sc.append(hsi.prob_of_use_gt_new(self.old_ht, self.old_LPcount, self.old_htHW, self.old_HWcount))

              if self.verbose:
                 print ('age score:',round(self.age_sc[-1],2), self.lb[0], self.ub[0], self.bn[0])
                 print ('height score:' , round(self.hw_sc[-1],2), self.lb[1], self.ub[1], self.bn[1])
                 #print 'dbh score:' , round(self.dbh_sc[-1],2),lb[2], ub[2], bn[2]
                 print ('ageHW score:',round(self.ageHW_sc[-1],2), self.lb[3], self.ub[3], self.bn[3])
                 print ('height score:' , round(self.hwHW_sc[-1], 2), self.lb[4], self.ub[4], self.bn[4])
                 print ('squirrel score:',round(self.sq_sc[-1],2))
                 print ('GT score:',round(self.gt_sc[-1],2))
                 print ('----------------')
                 print ('')               

              self.hw_gt_age.append(count_trees(self.old_ageHW,self.old_HWcount,self.tree_mature_age))
              self.hw_tot_age.append(count_trees(self.old_ageHW,self.old_HWcount,0))
              self.lp_gt_age.append(count_trees(self.old_age,self.old_LPcount,self.tree_mature_age))
              self.lp_tot_age.append(count_trees(self.old_age,self.old_LPcount,0))
        #print 'HW trees > 30:',count_trees(self.old_ageHW,self.old_HWcount,30)
        #print 'HW total:',count_trees(self.old_ageHW,self.old_HWcount,0)
        #print 'LLP trees > 30:',count_trees(self.old_age,self.old_LPcount,30)
        #print 'LLP total:',count_trees(self.old_age,self.old_LPcount,0)
    
