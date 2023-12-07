import numpy as np
class Param:
    def __init__(self, eta, n, sigmaY):
        self.eta    = eta
        self.n      = n
        self.sigmaY = sigmaY

    def vectorize(self):
        return np.array([self.eta, self.n, self.sigmaY])

    def display_status(self):
        print("*------------ Param ------------*")
        print("eta = ", self.eta)
        print("n = ", self.n)
        print("sigmaY  = ", self.sigmaY)
