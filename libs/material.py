from param import Param
class Material:
    def __init__(self, name: str, rho: float, m: Param):
        self.name = name
        self.param = m
        self.rho = rho