import numpy as np
from const import CGS, MKS
import copy
from conversion_function import *

class MaterialExtent:
    def __init__(self, extent_eta, extent_n, extent_sigmaY):
        self.extent_eta = copy.deepcopy(extent_eta)
        self.extent_n = copy.deepcopy(extent_n)
        self.extent_sigmaY = copy.deepcopy(extent_sigmaY)

CONST = {
    "cgs": CGS,
    "mks": MKS
}

const = CONST["mks"]

C_EXTENT = CExtent([5.0, 500.0])
LOG_S_EXTENT = CExtent([-6.0, 3.0])

def mat_hw_to_PL(eta_mks, n, sigmaY_mks, Hcm, Wcm):
    order = 2
    Theta_C = np.array([-1.37321779e+01,  8.14275110e-03, -1.09765321e+00, -8.00907812e-04,
        4.94833928e+00, -2.37988743e+00, -3.02972824e-03, -3.29177317e-01,
        3.79945757e-04,  2.93519945e-03, -2.55175315e-03, -2.57625369e-04,
       -1.09840507e-03,  1.04878754e+00, -1.27088316e-05,  3.57459552e-02,
       -3.86847818e-04,  1.06611942e+00, -1.72646031e-02,  4.11885913e-05,
        5.11073420e+00, -1.26161737e+00, -1.24167310e-04,  8.26656079e-03,
       -5.16803036e-08, -1.38268842e-05, -1.06641050e-05, -8.82504323e+01,
       -7.30494656e+00, -7.78391030e-03,  2.28734558e-05, -2.84809983e-05,
        2.59136974e-08,  5.31670213e-04,  2.13557773e+01,  4.34277901e-04,
       -1.92921774e-03, -6.04553328e-05,  4.17722294e+00,  2.49944076e-05,
        9.16347259e-03,  2.60463552e+01,  7.04927334e-05,  7.35523906e-06,
       -1.36230114e-05,  2.25393983e+01,  6.08146911e-07, -3.07427873e-01,
        1.29895721e-02, -9.46927385e-06, -1.90529838e-06, -4.84522719e-10,
        2.79472840e-02,  7.56647218e-02,  4.18255477e-08, -9.91701259e+00,
        6.77000023e-06,  1.72139168e-09,  1.83779067e-08,  4.39811697e-04,
        6.66829013e-08, -3.57121832e-04, -1.03725625e-04, -1.59476660e-09,
        2.57842478e-05,  5.92792423e-05])

    HW = np.array([Hcm * 0.01, Wcm * 0.01])
    C_pred = C_EXTENT.at(f_mat_scalar_compact_with_inverse([eta_mks, n, sigmaY_mks], HW, Theta_C, order))

    P0 = 2500.0
    L0 = C_pred / P0
    
    return np.array([P0, L0])



