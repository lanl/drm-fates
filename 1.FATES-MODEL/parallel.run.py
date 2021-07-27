#!/usr/bin/env python
"""
Computes a sample using templated run, example usage:

    mpirun -n 16 python parallel.run.py -c 'BCI.ICLM45ED.wolf.intel.Cb14cb81-F812a621.'  -r /lustre/scratch3/turquoise/cxu/ACME/cases  -f clm2.h0.2003-12.nc -n 1 -s 1 -t 1000 -l ./R/Missing.txt -g log
        
Phillip J. Wolfram
08/25/2017

Modified by Chonggang Xu
May 8th,2021

Modified by Rutuja Chitra-Tarak
July 24th, 2021
"""
import numpy as np
import shutil
import os
import subprocess
import logging

def run_case(casebase,runroot, finalfiletag, startitem,samplenum): # {{{

  casename = casebase + str(samplenum+startitem)
  logging.info(casename)
  rundir = runroot + "/"+casename + '/run'
  finalfile = casename + "." + finalfiletag 
  logging.info("file for check is " + rundir + '/' + finalfile)
#  test_array = np.zeros( (365, 5) )
#  for i in range(3):
#     test_array[:,2+i] = np.random.uniform(0,1,365)
#  test_array[:,0]=2000
#  test_array[:,1]= range(365)
#  test_array[:,1]= test_array[:,1]+1
#  header='year, doy, wetness1, wetness2, wetness3'
  
 # dir_path = os.path.dirname(os.path.realpath(__file__))
 # outdir=dir_path + '/OutputExtract'
 # if not os.path.exists(outdir):
 #   os.makedirs(outdir)
 # filename = outdir + '/HU.'+ str(samplenum+startitem)+'.csv'
 # np.savetxt(filename, test_array,comments='', header=header, delimiter=',') 
 
  if os.path.isdir(casename):
    if os.path.exists(rundir + '/' + finalfile):
      # case was already run succesfully!
      print 'Skipping case {}'.format(casename)
      logging.info('success')
      return
    else:
      # run the sample to get the results
      os.chdir(casename)
      print 'Running case {}'.format(casename)
      program = "./case.submit --no-batch"
      os.system(program)
      os.chdir('../')
      if os.path.exists(rundir + '/' + finalfile):
        logging.info('success')
      return 
  else:
    print('Case has not been built for {}'.format(casename))
  return
 

def subdivide_work(aproc, procsavail, nwork, overall=False):
  if procsavail > 1:
    nsubset = int(np.ceil(float(nwork)/float(procsavail)))
    idxs = aproc*nsubset
    idxe = (aproc+1)*nsubset
    samplerange = np.arange(nsubset*procsavail)[idxs:idxe]
  else:
    if overall:
      samplerange = np.arange(nwork)
    else:
      samplerange = np.array([aproc])

  return samplerange


def compute_results(casebase,runroot, finalfiletag, rank, startitem, totalitems, totalprocs,indxarr,overall=False):

    samplerange = subdivide_work(rank, totalprocs, totalitems,overall)
    for asample in samplerange:
      sample_num = -1
      if len(indxarr) > 0:
           if(asample<len(indxarr)):
               sample_num = indxarr[asample]
      else:
           sample_num = asample+1

      if sample_num < 0:
           print('Desired 1-indexed sample number of {} is skipped for processor of rank {}'.format(sample_num,rank))
           return
      
      run_case(casebase,runroot, finalfiletag, startitem,sample_num)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--casebase', dest='casebase',
                        help='base name of a case (casename=casebase${samplenumber+start_item})', metavar='STRING', required=True)
    parser.add_argument('-r', '--runroot', dest='runroot',
                        help='root folder to host run output for differnt cases.', metavar='FOLDER', required=True)
    parser.add_argument('-f', '--finalfiletag', dest='finalfiletag',
                        help='file tag to check if a run is finished.', metavar='STRING', required=True)			
    parser.add_argument('-n', '--sample_number', dest='sample_number', type=int,
                        help='Sample number to compute (1-indexed)', metavar='INT', required=False)
    parser.add_argument('-s', '--start_item', dest='start_item', type=int,
                        help='Start sample to be computed (default=1)', metavar='INT', required=False)			
    parser.add_argument('-t', '--total_items', dest='total_items', type=int,
                        help='Total number of samples to be computed', metavar='INT', required=True)
    parser.add_argument('-l', '--file', dest='file', type=str,
                        help='The file with the list of case numbers', metavar='FILE', required=False)
    parser.add_argument('-g', '--log', dest='log', type=str,
                        help='The log file', metavar='FILE', required=False)
    args = parser.parse_args()
    
    if args.start_item is None:
       args.start_item = 0
    else:
       args.start_item = args.start_item - 1
           
    indxarr = []
    if args.file is not None:
       txf = open(args.file,"r")
       lines = txf.readlines()
       indxarr = [0]*len(lines)
       args.start_item = 0
       for linei in range(0,len(lines)):
          indxarr[linei] = int(lines[linei].strip())          

    if args.log is None:
        args.log = 'log'
    
    logging.basicConfig(filename=args.log,filemode="w",level=logging.INFO,format='%(message)s')
    logging.info('Started')
    
    if args.sample_number is not None:
      args.sample_number = args.sample_number - 1
      rank = args.sample_number
      totalprocs = 1
      overall = False
    else:
      from mpi4py import MPI
      comm = MPI.COMM_WORLD
      rank = comm.Get_rank()
      totalprocs = comm.Get_size()
      overall = False
      if totalprocs < 2:
         overall = True
      logging.info('Computing using rank {} of {}'.format(rank, totalprocs))
      
    compute_results(args.casebase, args.runroot,args.finalfiletag, rank, args.start_item, args.total_items, totalprocs, indxarr,overall)

