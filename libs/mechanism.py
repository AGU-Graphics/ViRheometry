import math
import numpy as np
import copy
import shutil
from setup import Setup
import copy
import numpy as np
import sys
import os
from const import CGS
const = CGS
MAX_ETA, MAX_N, MAX_SIGMA_Y, MIN_ETA, MIN_N, MIN_SIGMA_Y = const.MAX_ETA, const.MAX_N, const.MAX_SIGMA_Y, const.MIN_ETA, const.MIN_N, const.MIN_SIGMA_Y
WEIGHTS_STEPS = const.WEIGHTS_STEPS
import time
from compare_loss import mat_hw_to_PL

extent_eta = const.extent_eta
extent_n = const.extent_n
extent_sigmaY = const.extent_sigmaY

MIN_W = const.MIN_W
MAX_W = const.MAX_W

MIN_H = const.MIN_H
MAX_H = const.MAX_H

class InvalidParamException(Exception):
    def __init__(self, param):
        self.param = param
    
    def __str__(self):
        return (
            f"[{self.param.eta, self.param.n, self.param.sigmaY}] is not valid on Loss(or Hesse). Check the range of params."
        )

class MaterialExtent:
    def __init__(self, extent_eta, extent_n, extent_sigmaY):
        self.extent_eta = copy.deepcopy(extent_eta)
        self.extent_n = copy.deepcopy(extent_n)
        self.extent_sigmaY = copy.deepcopy(extent_sigmaY)

class SimpleMechanism:
    def __init__(self, const):
        self.const = const

    def isSetupValid(self, P, L, m):
        if P*L - m.sigmaY <= 0.0:
            return False
        else:
            return True

    def isSetupValid2(self, P, L, m):
        eta, n, sigmaY = m.eta, m.n, m.sigmaY
        l = sigmaY / P
        W = P * (L - l) / eta
        if math.isnan(W) or W <= 0.0:
            return False
        else:
            return True
            
    def singleLoss(self, m, m_star, P, L): 
        if not self.isSetupValid(P,L,m) or not self.isSetupValid2(P,L,m):
            raise InvalidParamException(m)

        eta, n, sigmaY = m.eta, m.n, m.sigmaY
        eta_star, n_star, sigmaY_star = m_star.eta, m_star.n, m_star.sigmaY

        U1 = n / (n + 1.0)
        U2 = math.pow(P / eta, 1.0 / n)
        U3 = math.pow(L - sigmaY / P, (n + 1.0) / n)

        U1_star = n_star / (n_star + 1.0)
        U2_star = math.pow(P / eta_star, 1.0 / n_star)
        U3_star = math.pow(L - sigmaY_star / P, (n_star + 1.0) / n_star)

        D_inner_denom = P * L - sigmaY
        D_pwr = (n + 1.0) / n

        D_star_inner_denom = P*L - sigmaY_star
        D_star_pwr = (n_star + 1.0) / n_star

        num_intervals = 100
        dy = L / num_intervals
        Loss = 0.0
        for i in range(num_intervals):
            y = L * (i+0.5) / num_intervals

            D_inner_numer = max(0.0, P * y - sigmaY)
            D = 1.0 - math.pow(D_inner_numer / D_inner_denom, D_pwr)
            u = D * U1 * U2 * U3

            D_star_inner_numer = max(0.0, P * y - sigmaY_star)
            D_star = 1.0 - math.pow(D_star_inner_numer / D_star_inner_denom, D_star_pwr)
            u_star = D_star * U1_star * U2_star * U3_star

            Loss += (u - u_star) * (u - u_star) * dy

        return Loss

    def singleHessian(self, m, P, L):
        if not self.isSetupValid(P,L,m) or not self.isSetupValid2(P,L,m):
            raise InvalidParamException(m)
        eta, n, sigmaY = m.eta, m.n, m.sigmaY

        l = sigmaY / P
        W = P * (L - l) / eta

        S        = math.pow(W, (n+1.0)/n) * eta / (P * (n+1.0))
        A_eta    = -S / eta
        A_n      = S * (1.0/(n+1.0) - math.log(W) / n)
        B_n      = S / n
        A_sigmaY = -math.pow(W, 1.0/n) / P

        C1 = (1.0+n)*(1.0+n) / ((1.0+2.0*n) * (2.0+3.0*n))
        C2 = n*n*(3.0+5.0*n)*(1.0+n) / ((1.0+2.0*n)*(1.0+2.0*n)*(2.0+3.0*n)*(2.0+3.0*n))
        C3 = (2.0+3.0*n) / (2.0*(1.0+n)*(1.0+2.0*n))
        C4 = n*n*n/((2.0+3.0*n)*(2.0+3.0*n)*(2.0+3.0*n))
        C5 = n*n*(3.0+4.0*n)/(4.0*(1.0+n)*(1.0+n)*(1.0+2.0*n)*(1.0+2.0*n))
        C6 = 1.0 / ((1.0+n)*(2.0+n))

        H_eta_eta       = 2.0 * A_eta * A_eta * l + 4.0 * A_eta * A_eta * (L-l) * C1
        H_eta_n         = 2.0 * A_eta * A_n * l + 4.0 * A_eta * A_n * (L-l) * C1 - 2.0 * A_eta * B_n * (L-l) * C2
        H_eta_sigmaY    = 2.0 * A_eta * A_sigmaY * l + 2.0 * A_eta * A_sigmaY * (L-l) * C3
        H_n_n           = 2.0 * A_n * A_n * l + 4.0 * B_n * B_n * (L-l) * C4 - 4.0 * A_n * B_n * (L-l) * C2 + 4.0 * A_n * A_n * (L-l) * C1
        H_n_sigmaY      = 2.0 * A_n * A_sigmaY * l + 2.0 * A_n * A_sigmaY * (L-l) * C3 - 2.0 * B_n * A_sigmaY * (L-l) * C5
        H_sigmaY_sigmaY = 2.0 * A_sigmaY * A_sigmaY * l + 4.0 * A_sigmaY * A_sigmaY * (L-l) * C6

        return np.matrix([[H_eta_eta, H_eta_n, H_eta_sigmaY],[H_eta_n, H_n_n, H_n_sigmaY],[H_eta_sigmaY, H_n_sigmaY, H_sigmaY_sigmaY]])

    def rescaleHessian(self, Hessian):
        dm_dmtilde      = np.zeros((3,3))
        dm_dmtilde[0,0] = self.const.extent_eta[1] - self.const.extent_eta[0]
        dm_dmtilde[1,1] = self.const.extent_n[1] - self.const.extent_n[0]
        dm_dmtilde[2,2] = self.const.extent_sigmaY[1] - self.const.extent_sigmaY[0]
        return dm_dmtilde @ Hessian @ dm_dmtilde.transpose()

    def totalLossPL(self,m,m_star,P,L):
        loss = 0.0
        try:
            loss = self.singleLoss(m, m_star, P, L)
        except Exception as e:
            print(str(e))  
            return float('NaN')
        return loss
    
class Mechanism:

    def isSetupValid(self, P, L, m):
        if P*L - m.sigmaY <= 0.0:
            return False
        else:
            return True

    def isSetupValid2(self, P, L, m):
        eta, n, sigmaY = m.eta, m.n, m.sigmaY
        l = sigmaY / P
        W = P * (L - l) / eta
        if math.isnan(W) or W <= 0.0:
            return False
        else:
            return True
        
    def singleLoss(self, m, m_star, P, L): 
        if not self.isSetupValid(P,L,m) or not self.isSetupValid2(P,L,m):
            raise InvalidParamException(m)

        eta, n, sigmaY = m.eta, m.n, m.sigmaY
        eta_star, n_star, sigmaY_star = m_star.eta, m_star.n, m_star.sigmaY

        U1 = n / (n + 1.0)
        U2 = math.pow(P / eta, 1.0 / n)
        U3 = math.pow(L - sigmaY / P, (n + 1.0) / n)

        U1_star = n_star / (n_star + 1.0)
        U2_star = math.pow(P / eta_star, 1.0 / n_star)
        U3_star = math.pow(L - sigmaY_star / P, (n_star + 1.0) / n_star)

        D_inner_denom = P * L - sigmaY
        D_pwr = (n + 1.0) / n

        D_star_inner_denom = P*L - sigmaY_star
        D_star_pwr = (n_star + 1.0) / n_star

        num_intervals = 100
        dy = L / num_intervals
        Loss = 0.0
        for i in range(num_intervals):
            y = L * (i+0.5) / num_intervals

            D_inner_numer = max(0.0, P * y - sigmaY)
            D = 1.0 - math.pow(D_inner_numer / D_inner_denom, D_pwr)
            u = D * U1 * U2 * U3

            D_star_inner_numer = max(0.0, P * y - sigmaY_star)
            D_star = 1.0 - math.pow(D_star_inner_numer / D_star_inner_denom, D_star_pwr)
            u_star = D_star * U1_star * U2_star * U3_star

            Loss += (u - u_star) * (u - u_star) * dy

        return Loss

    def singleHessian(self, m, P, L):
        if not self.isSetupValid(P,L,m) or not self.isSetupValid2(P,L,m):
            raise InvalidParamException(m)
        eta, n, sigmaY = m.eta, m.n, m.sigmaY

        l = sigmaY / P
        W = P * (L - l) / eta

        S        = math.pow(W, (n+1.0)/n) * eta / (P * (n+1.0))
        A_eta    = -S / eta
        A_n      = S * (1.0/(n+1.0) - math.log(W) / n)
        B_n      = S / n
        A_sigmaY = -math.pow(W, 1.0/n) / P

        C1 = (1.0+n)*(1.0+n) / ((1.0+2.0*n) * (2.0+3.0*n))
        C2 = n*n*(3.0+5.0*n)*(1.0+n) / ((1.0+2.0*n)*(1.0+2.0*n)*(2.0+3.0*n)*(2.0+3.0*n))
        C3 = (2.0+3.0*n) / (2.0*(1.0+n)*(1.0+2.0*n))
        C4 = n*n*n/((2.0+3.0*n)*(2.0+3.0*n)*(2.0+3.0*n))
        C5 = n*n*(3.0+4.0*n)/(4.0*(1.0+n)*(1.0+n)*(1.0+2.0*n)*(1.0+2.0*n))
        C6 = 1.0 / ((1.0+n)*(2.0+n))

        H_eta_eta       = 2.0 * A_eta * A_eta * l + 4.0 * A_eta * A_eta * (L-l) * C1
        H_eta_n         = 2.0 * A_eta * A_n * l + 4.0 * A_eta * A_n * (L-l) * C1 - 2.0 * A_eta * B_n * (L-l) * C2
        H_eta_sigmaY    = 2.0 * A_eta * A_sigmaY * l + 2.0 * A_eta * A_sigmaY * (L-l) * C3
        H_n_n           = 2.0 * A_n * A_n * l + 4.0 * B_n * B_n * (L-l) * C4 - 4.0 * A_n * B_n * (L-l) * C2 + 4.0 * A_n * A_n * (L-l) * C1
        H_n_sigmaY      = 2.0 * A_n * A_sigmaY * l + 2.0 * A_n * A_sigmaY * (L-l) * C3 - 2.0 * B_n * A_sigmaY * (L-l) * C5
        H_sigmaY_sigmaY = 2.0 * A_sigmaY * A_sigmaY * l + 4.0 * A_sigmaY * A_sigmaY * (L-l) * C6

        return np.matrix([[H_eta_eta, H_eta_n, H_eta_sigmaY],[H_eta_n, H_n_n, H_n_sigmaY],[H_eta_sigmaY, H_n_sigmaY, H_sigmaY_sigmaY]])

    def rescaleHessian(self, Hessian):
        dm_dmtilde      = np.zeros((3,3))
        dm_dmtilde[0,0] = extent_eta[1] - extent_eta[0]
        dm_dmtilde[1,1] = extent_n[1] - extent_n[0]
        dm_dmtilde[2,2] = extent_sigmaY[1] - extent_sigmaY[0]
        return dm_dmtilde @ Hessian @ dm_dmtilde.transpose()

    def computeRescaledNormal(self, H, mat_ext):
        s, Q = np.linalg.eig(H)
        n = Q[:,np.argmax(s)]
        n = np.ravel(n)

        A = np.zeros((3, 3))
        A[0][0] = (mat_ext.extent_eta[1] - mat_ext.extent_eta[0])
        A[1][1] = (mat_ext.extent_n[1] - mat_ext.extent_n[0])
        A[2][2] = (mat_ext.extent_sigmaY[1] - mat_ext.extent_sigmaY[0])

        n_tilde = A @ n
        n_tilde = n_tilde / np.linalg.norm(n_tilde)

        return n_tilde

    def searchNewSetup_orthognality_for_forth_setup(self, m, S):
        mat_ext = MaterialExtent([const.MIN_ETA, const.MAX_ETA], [const.MIN_N, const.MAX_N], [const.MIN_SIGMA_Y, const.MAX_SIGMA_Y])
        s1 = S[0]
        P0, L0 = mat_hw_to_PL(m.eta * 0.1, m.n, m.sigmaY * 0.1, s1.H, s1.W)
        P, L = P0 / 10.0, L0 * 100.0
        H1 = self.singleHessian(m, P, L)
        q1 = self.computeRescaledNormal(H1, mat_ext)

        s2 = S[1]
        P0, L0 = mat_hw_to_PL(m.eta * 0.1, m.n, m.sigmaY * 0.1, s2.H, s2.W)
        P, L = P0 / 10.0, L0 * 100.0
        H2 = self.singleHessian(m, P, L)
        q2 = self.computeRescaledNormal(H2, mat_ext)

        s3 = S[2]
        print("s3: ")
        s3.display_status()
        P0, L0 = mat_hw_to_PL(m.eta * 0.1, m.n, m.sigmaY * 0.1, s3.H, s3.W)
        P, L = P0 / 10.0, L0 * 100.0
        H3 = self.singleHessian(m, P, L)
        q3 = self.computeRescaledNormal(H3, mat_ext)

        print("q1: ", q1)
        print("q2: ", q2)
        print("q3: ", q3)

        min_proj = sys.float_info.max
        best_setup = Setup(100,100,1.0)
        Hx10_list = np.arange(int(MIN_H*10), int(MAX_H*10)+1)
        Wx10_list = np.arange(int(MIN_W*10), int(MAX_W*10)+1)
        candidate_of_setups = []
        for Hx10 in Hx10_list:
            for Wx10 in Wx10_list:
                candidate_of_setups.append(Setup(Hx10 * 0.1, Wx10 * 0.1, 1.0))

        for s in candidate_of_setups:
            Hcm = s.H
            Wcm = s.W
            P0, L0 = mat_hw_to_PL(m.eta * 0.1, m.n, m.sigmaY * 0.1, Hcm, Wcm)
            P, L = P0 / 10.0, L0 * 100.0
            if self.isSetupValid(P,L,m):
                H_dash = np.zeros((3,3))
                H_dash = self.singleHessian(m, P, L)
                q_dash = self.computeRescaledNormal(H_dash, mat_ext)
                proj_1 = np.dot(q1, q_dash)
                proj_2 = np.dot(q2, q_dash)
                proj_3 = np.dot(q3, q_dash)
                print("proj 1 > 2 > 3: ", proj_1, " > ", proj_2, " > ", proj_3)
                proj = np.sqrt(proj_1*proj_1 + proj_2*proj_2 + proj_3*proj_3)
                if min_proj > proj:
                    min_proj = proj
                    best_setup = Setup(Hcm, Wcm, 1.0)
            else:
                continue
        best_setup.display_status()
        s1.weight = 1/4
        s2.weight = 1/4
        s3.weight = 1/4
        best_setup.weight = 1/4
        new_setups = [s1, s2, s3, best_setup]
        for s in new_setups:
            s.display_status()
        return new_setups

    def searchNewSetup_orthognality_for_third_setup(self, m, S):
        mat_ext = MaterialExtent([const.MIN_ETA, const.MAX_ETA], [const.MIN_N, const.MAX_N], [const.MIN_SIGMA_Y, const.MAX_SIGMA_Y])
        s1 = S[0]
        P0, L0 = mat_hw_to_PL(m.eta * 0.1, m.n, m.sigmaY * 0.1, s1.H, s1.W)
        P, L = P0 / 10.0, L0 * 100.0
        H1 = self.singleHessian(m, P, L)
        q1 = self.computeRescaledNormal(H1, mat_ext)

        s2 = S[1]
        P0, L0 = mat_hw_to_PL(m.eta * 0.1, m.n, m.sigmaY * 0.1, s2.H, s2.W)
        P, L = P0 / 10.0, L0 * 100.0
        H2 = self.singleHessian(m, P, L)
        q2 = self.computeRescaledNormal(H2, mat_ext)

        min_proj = sys.float_info.max
        best_setup = Setup(100,100,1.0)
        Hx10_list = np.arange(int(MIN_H*10), int(MAX_H*10)+1)
        Wx10_list = np.arange(int(MIN_W*10), int(MAX_W*10)+1)
        candidate_of_setups = []
        for Hx10 in Hx10_list:
            for Wx10 in Wx10_list:
                candidate_of_setups.append(Setup(Hx10 * 0.1, Wx10 * 0.1, 1.0))
        for s in candidate_of_setups:
            Hcm = s.H
            Wcm = s.W
            P0, L0 = mat_hw_to_PL(m.eta * 0.1, m.n, m.sigmaY * 0.1, Hcm, Wcm)
            P, L = P0 / 10.0, L0 * 100.0
            if self.isSetupValid(P,L,m):
                H_dash = np.zeros((3,3))
                H_dash = self.singleHessian(m, P, L)
                q_dash = self.computeRescaledNormal(H_dash, mat_ext)
                proj_1 = np.dot(q1, q_dash)
                proj_2 = np.dot(q2, q_dash)
                proj = np.sqrt(proj_1*proj_1 + proj_2*proj_2)
                if min_proj > proj:
                    min_proj = proj
                    best_setup = Setup(Hcm, Wcm, 1.0)
            else:
                continue
        best_setup.display_status()
        s1.weight = 1/3
        s2.weight = 1/3
        best_setup.weight = 1/3
        new_setups = [s1, s2, best_setup]
        for s in new_setups:
            s.display_status()
        return new_setups


    def searchNewSetup_orthognality_for_second_setup(self, m, S):
        mat_ext = MaterialExtent([const.MIN_ETA, const.MAX_ETA], [const.MIN_N, const.MAX_N], [const.MIN_SIGMA_Y, const.MAX_SIGMA_Y])
        s1 = S[0]
        P0, L0 = mat_hw_to_PL(m.eta * 0.1, m.n, m.sigmaY * 0.1, s1.H, s1.W)
        P, L = P0 / 10.0, L0 * 100.0
        H1 = self.singleHessian(m, P, L)
        q1 = self.computeRescaledNormal(H1, mat_ext)

        min_proj = sys.float_info.max
        best_setup = Setup(100,100,1.0)

        Hx10_list = np.arange(int(MIN_H*10), int(MAX_H*10)+1)
        Wx10_list = np.arange(int(MIN_W*10), int(MAX_W*10)+1)
        candidate_of_setups = []
        for Hx10 in Hx10_list:
            for Wx10 in Wx10_list:
                candidate_of_setups.append(Setup(Hx10 * 0.1, Wx10 * 0.1, 1.0))
                
        for s in candidate_of_setups:
            Hcm = s.H
            Wcm = s.W
            P0, L0 = mat_hw_to_PL(m.eta * 0.1, m.n, m.sigmaY * 0.1, Hcm, Wcm)
            P, L = P0 / 10.0, L0 * 100.0

            if self.isSetupValid(P,L,m):
                H_dash = np.zeros((3,3))

                H_dash = self.singleHessian(m, P, L)
                q_dash = self.computeRescaledNormal(H_dash, mat_ext)
                proj = np.abs(np.dot(q1, q_dash))
                if min_proj > proj:
                    min_proj = proj
                    best_setup = Setup(Hcm, Wcm, 1.0)
            else:
                continue
        
        best_setup.display_status()
        s1.weight = 0.5
        best_setup.weight = 0.5
        new_setups = [s1, best_setup]

        return new_setups

class Mechanism3D(Mechanism):

    def totalLoss(self, m, setups, setting_sim_sh_path, ref_dir_paths, sim_dir_paths, skin_app_path, gl_app_path, eval_count, rho):
        loss = 0.0
        for i,s in enumerate(setups):
            print("iter: ", i)
            m.display_status()
            s.display_status()
            loss += s.weight * self.singleLoss(m, s.H, s.W, setting_sim_sh_path,\
                 ref_dir_paths[i], sim_dir_paths[i], skin_app_path, gl_app_path, eval_count, rho)
        return loss

    def singleLoss(self, m, H, W, shellScriptPath, ref_dir_path, sim_dir_path, skin_app_path, gl_app_path, eval_count, rho): 
        max_frames = 8

        f = open(sim_dir_path+"/diffvalues.dat", 'w')
        f.close()

        makeXml = 'sh "' + shellScriptPath + '" "' + sim_dir_path + '/to_process' + '" ' + str(eval_count) + ' '\
        + str(m.eta) + ' ' + str(m.n) + ' ' + str(m.sigmaY) + ' '\
        + str(H) + ' ' + str(W) + ' "'\
        + skin_app_path + '" "' + gl_app_path + '" "'\
        + ref_dir_path + '" "' + sim_dir_path + '" ' + str(rho)

        os.system(makeXml)

        while True:
            with open(sim_dir_path+"/diffvalues.dat") as f:
                diff_values = np.loadtxt(f)
                if diff_values.size >= max_frames:
                    shutil.copyfile(sim_dir_path+"/snapcal_08.png", sim_dir_path+"/backup_sim_endframe/endframe_"+str(eval_count)+".png")
                    return np.mean(diff_values)
                else:
                    time.sleep(3)
