import hviscore as hvi
import numpy as np
# Nosetests developed to test functions in hviscore.py
# Example: >> nosetests hvi_unit_tests.py
#          >> nosetests -v hvi_unit_tests.py (verbose on)

def test_norm_gaussian():
	# Note: first term in the Gaussian has to be float
	x=np.linspace(-50,50,1000)
	[y0, y]=hvi.norm_gaussian(1.0,x,1)
	[y1, y]=hvi.norm_gaussian(-1.0,x,1)
	assert y0 == y1
	[y0, y]=hvi.norm_gaussian(0.0,x,1)
	assert y0 == 1

def test_dbh2ba():
	assert hvi.dbh2ba(10) == 78.53981633974483

def test_gscore():
	# benchmarking number =11 , mean = 0, sigma = 2
	# equal spacing from mean, should be the same
    a = hvi.gscore(12, 11, 0, 2)
    b = hvi.gscore(10, 11, 0, 2)
    assert a == b
    
def test_age2dbh():
	# input is in years
	# output is in meters 
	assert hvi.age2dbh(0.0) == 0.0
	assert hvi.age2dbh(50) == 0.2827
	assert hvi.age2dbh(125) == 1.0402
	
def test_subset_matrix():
	a = np.matrix('1 2 3 ; 3 4 5; 6 7 8')
	b = np.matrix(' 2 3 ;  4 5')
	# 0 and 1 are indexes in the matrix
	c = hvi.subset_matrix(a, 0, 1, 2, 2)
	assert np.array_equal(c,b)
	
def test_score_relaxed():
	# Note, that test is not symmetric w.r.t. bn 
	a = np.matrix('1 2 3 ; 3 4 5; 6 7 8')
	assert hvi.score_relaxed(a, 2, 6, 6, 0) == 1
	assert hvi.score_relaxed(a, 2, 6, 3, 0) == 0.5
	assert hvi.score_relaxed(a, 2, 6, 12, 0) == 0.5
	
def test_gscore_relaxed():
	# Note, that test is not symmetric w.r.t. bn 
	a = np.matrix('1 2 3 ; 3 4 5; 6 7 8')
	assert hvi.gscore_relaxed(a, 2, 6, 6, 0) == 1
	b = hvi.gscore_relaxed(a, 2, 6, 8, 0)
	assert hvi.gscore_relaxed(a, 2, 6, 4, 0) == b
	
def test_score_subset():	
	x = np.ones([20,20])
	x1= np.cumsum(x)
	x2=np.reshape(x1, (20, 20))  
	
	#Work with all relaxeds score 
	assert hvi.score_subset(x2, 5, 10, 7, 5, 5, 1) == 0.0625
	
	#works with the gscore
	#assert hvi.score_subset(x2, 5, 10, 7, 5, 5, 1) == 0.06762248
	
def	test_score1d():
	# Note, that this test is valid if and only if it returns gscore(num, bn, 1)
	# The sm needs to be added to the score1d function
	x = np.ones([3,3])
	x1= np.cumsum(x)
	a = hvi.score1d(x1,3,6,2)
	assert hvi.score1d(x1,3,6,6) == a
	
	
	