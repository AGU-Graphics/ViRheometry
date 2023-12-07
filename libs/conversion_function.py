import numpy as np
import math
import copy

EPS = 0.0001

class CExtent:
	def __init__(self, extent):
		self.extent = copy.deepcopy(extent)
	def at(self, rate):
		return rate * (self.extent[1] - self.extent[0]) + self.extent[0]
	def get_rate(self, c):
		return (c - self.extent[0]) / (self.extent[1] - self.extent[0])

def f_mat_scalar_compact_with_inverse(M, HW, Theta, order):
	MHW = np.concatenate([M, HW])
	ret = Theta[0]
	prev = 1
	# the first ten prime numbers
	test_v = np.array([2, 3, 5, 7, 11, 13, 17, 19, 23, 29])

	MHW_with_inverse = np.concatenate([MHW, np.array([1.0 / MHW[0], 1.0 / MHW[1], 1.0 / (MHW[2] + EPS), 1.0 / MHW[3], 1.0 / MHW[4]])])
	for k in range(1, order+1):
		prod_t = test_v
		prod = MHW_with_inverse
		for i in range(k-1):
			prod_t = np.tensordot(prod_t, test_v, axes=0)
			prod = np.tensordot(prod, MHW_with_inverse, axes=0)

		prod = prod.reshape(-1)
		prod_t = prod_t.reshape(-1)
		prod_t_uni, prod_t_idx, prod_t_c = np.unique(prod_t, return_index = True, return_counts = True)
		prod_uni = np.multiply(prod[prod_t_idx], prod_t_c)
		num_terms = len(prod_uni)
		w_vec = Theta[prev:(prev+num_terms)]
		prev += num_terms
		ret += (np.dot(w_vec, prod_uni)) / math.factorial(k)
	return ret
