import numpy as np
import functools
import csv
from param import Param
from mechanism import Mechanism, Mechanism3D
from const import CGS
const = CGS
extent_eta = const.extent_eta
extent_n = const.extent_n
extent_sigmaY = const.extent_sigmaY


class Optimizer:
    def __init__(self):
        self.param_data = []
        self.loss_data = []
        self.eta_values = []
        self.n_values = []
        self.sigmaY_values = []

    def optimize(self, cmaes, init_m, setups):
        fitfunc = functools.partial(self.fitfunc, setups=setups)
        m_breve_list = cmaes.optimize(fitfunc, [init_m.eta, init_m.n, init_m.sigmaY])
        m_breve = Param(m_breve_list[0], m_breve_list[1], m_breve_list[2])
        return m_breve

    def fitfunc(self, m_list, dim, setups):
        mechanism = Mechanism()
        m = Param(m_list[0],m_list[1],m_list[2])
        total_loss = mechanism.totalLoss(m, setups)
        for i in range(len(setups)):
            self.loss_data.append(total_loss)
            self.param_data.append(m)
            self.eta_values.append(m.eta)
            self.n_values.append(m.n)
            self.sigmaY_values.append(m.sigmaY)

        return total_loss
    

class Optimizer3D(Optimizer):
    def __init__(self, outdir, setting_sim_sh_path, ref_dir_paths, sim_dir_paths_list, \
        skin_app_path, gl_app_path, rho):
        super().__init__()
        self.setting_sim_sh_path = setting_sim_sh_path
        self.ref_dir_paths = ref_dir_paths
        self.sim_dir_paths_list = sim_dir_paths_list
        self.skin_app_path = skin_app_path
        self.gl_app_path = gl_app_path
        self.eval_count = 0
        self.loss_data_file = outdir+"/loss_data.csv"
        self.param_data_file = outdir+"/param_data.csv"
        self.best_loss_data_file = outdir + "/best_loss_data.csv"
        self.outdir = outdir
        self.rho = rho

    def fitfunc(self, m_list, dim, setups):
        self.eval_count += 1
        mechanism3D = Mechanism3D()
        m = Param(m_list[0],m_list[1],m_list[2])
        total_loss = mechanism3D.totalLoss(m, setups, self.setting_sim_sh_path,\
            self.ref_dir_paths, self.sim_dir_paths_list[-1], self.skin_app_path, self.gl_app_path, self.eval_count,self.rho)
        for i in range(len(setups)):
            self.loss_data.append(total_loss)
            self.saveLossData(total_loss)
            self.param_data.append(m)
            self.saveParamData(m)
            self.saveBestLossData()
        return total_loss

    def saveBestLossData(self):
        with open(self.best_loss_data_file, 'a') as f:
            writer = csv.writer(f)
            writer.writerow([np.amin(self.loss_data)])        

    def saveLossData(self, loss):
        with open(self.loss_data_file, 'a') as f:
            writer = csv.writer(f)
            writer.writerow([loss])

    def saveParamData(self, param):
        with open(self.param_data_file, 'a') as f:
            writer = csv.writer(f)
            writer.writerow([param.eta, param.n, param.sigmaY])

