import sys
sys.path.append("../libs")
import os
import argparse
import random
from plotly.subplots import make_subplots
import numpy as np
import random
import csv
import xmlParser
import time
import const
import copy
from optimizer import Optimizer3D
from cmaes import CMAES_PWQB_LS
from cmaes import CMAES
from mechanism import Mechanism, Mechanism3D
from setup import Setup
from param import Param
from const import CGS
const = CGS
MAX_ETA, MAX_N, MAX_SIGMA_Y, MIN_ETA, MIN_N, MIN_SIGMA_Y = const.MAX_ETA, const.MAX_N, const.MAX_SIGMA_Y, const.MIN_ETA, const.MIN_N, const.MIN_SIGMA_Y
MAX_LOSS_THR, MIN_H, MAX_H, MIN_W, MAX_W = const.MAX_LOSS_THR, const.MIN_H, const.MAX_H, const.MIN_W, const.MAX_W

import shutil

extent_eta = const.extent_eta
extent_n = const.extent_n
extent_sigmaY = const.extent_sigmaY
max_frame = 8

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--first_dir')
args = parser.parse_args()

class Main:
    def __init__(self, xmlData):
        self.root_folder = os.path.abspath(xmlData.pathData.root_dir_path)
        self.sh_dir_path = os.path.abspath(xmlData.pathData.shell_script_dir_path)
        self.particle_skinner_path = os.path.abspath(xmlData.pathData.particle_skinner_path)
        self.mpm_path = os.path.abspath(xmlData.pathData.mpm_path)
        self.sim_render_path = os.path.abspath(xmlData.pathData.GL_render_path)
        self.rho = xmlData.setupData.rho
        
        self.setups = []
        H = xmlData.setupData.H
        W = xmlData.setupData.W
        self.setups.append(Setup(H,W,1.0))
        self.setups[0].display_status()

        self.sim_dir_name = "sim"
        self.ref_dir_name = "ref"
        self.ref_dir_paths = []
        self.sim_dir_paths_list = []
        self.makeOptimizationDir()
        self.taichiActivation()

    def taichiActivation(self):
        script ='sh "'+ self.sh_dir_path + '/run_taichi_3d.sh" "' + self.root_folder + '" "' + self.mpm_path + '" ' + str(self.rho)
        os.system(script)

    def makeOptimizationDir(self):
        if not os.path.isdir(self.root_folder):
            os.mkdir(self.root_folder)
            os.mkdir(self.root_folder+"/"+self.sim_dir_name)
            os.mkdir(self.root_folder+"/"+self.ref_dir_name)
        else:
            print("this directory has already exist.")
            while True:
                yorN = input("overwrite and start the optimization? 'yes' or 'no' [y/N]: ").lower()
                if yorN in ['y', 'ye', 'yes']:
                    os.system("rm -r " + self.root_folder)
                    os.mkdir(self.root_folder)
                    os.mkdir(self.root_folder+"/"+self.sim_dir_name)
                    os.mkdir(self.root_folder+"/"+self.ref_dir_name)
                    break
                elif yorN in ['n', 'no', ""]:
                    exit(1)
            

    def makeReferenceDir(self, n_setup):
        dir = self.root_folder+"/"+self.ref_dir_name+"/"+str(n_setup)+"setup"
        os.mkdir(dir)  
        os.mkdir(dir+"/to_process")
        os.mkdir(dir+"/processing")
        os.mkdir(dir+"/processed")
        self.ref_dir_paths.append(dir)

    def makeSimulationDir(self, num_setups):
        dir = self.root_folder+"/"+self.sim_dir_name+"/"+str(num_setups)+"setups"
        os.mkdir(dir)
        process_dirs = []
        for n in range(num_setups):
            process_dir = dir+"/"+str(n+1)+"setup"
            process_dirs.append(process_dir)
            os.mkdir(process_dir)
            os.mkdir(process_dir+"/to_process")
            os.mkdir(process_dir+"/processing")
            os.mkdir(process_dir+"/processed")
            os.mkdir(process_dir+"/backup_sim_endframe")
        self.sim_dir_paths_list.append(process_dirs)

    def waitForFile(self, file_path):
        while True:
            if os.path.isfile(file_path):
                break
            else:
                print("waiting file : "+file_path)
                time.sleep(1)

    def trueParamSimulation(self):
        self.makeReferenceDir(len(self.setups))

        from_dir = args.first_dir
        fname_list = [filename for filename in os.listdir(from_dir) if not filename.startswith('.')]

        for fname in fname_list:
            print(fname)
            if os.path.isfile(from_dir + "/" + fname):
                from_file = from_dir + "/" + fname
                to_file = self.ref_dir_paths[-1] + "/" + fname
                print("from: ",from_file, "to: ", to_file)
                shutil.copy(from_file, self.ref_dir_paths[-1])

        for frame in range(max_frame):
            frame_number = frame + 1
            self.waitForFile(self.ref_dir_paths[-1]+"/config_"+ str(frame_number).zfill(2)  + ".png")

    def exportSecondSetupXML(self, m_breve):
        material_name = os.path.basename(self.root_folder).split(",")[0].split("_")[0]
        s2 = self.setups[1]
        Hcm_2 = str(round(s2.H, 2))
        Wcm_2 = str(round(s2.W, 2))
        arg_1 = "../data/ref_" + material_name + "_" + Hcm_2 + "_" + Wcm_2 + "_2/"
        if not os.path.exists(arg_1):
            os.mkdir(arg_1)
        arg_2 = "settings.xml"
        arg_3 = "../data/" + material_name + "_2"
        arg_4 = Hcm_2
        arg_5 = Wcm_2
        arg_6 = str(self.rho)
        arg_7 = str(m_breve.eta)
        arg_8 = str(m_breve.n)
        arg_9 = str(m_breve.sigmaY)

        command = "sh out_xml.sh"
        command += " "
        command += arg_1
        command += " "
        command += arg_2
        command += " "
        command += arg_3
        command += " "
        command += arg_4
        command += " "
        command += arg_5
        command += " "
        command += arg_6
        command += " "
        command += arg_7
        command += " "
        command += arg_8
        command += " "
        command += arg_9

        os.system(command)
        explanation = "Before running the optimization for the 2nd experiment, please put the following files into the specified folder \"" + arg_1 + "\": \n"
        explanation += "  -  the binary frames (config01.png to config08.png) from the experiment with the setting W = " + arg_5 + ", and H = " + arg_4 + "\n"
        explanation += "  -  the camera parameter setting (camera_params.xml) from the calibration"
        print(explanation)

    def optimization_1(self):
        import codecs
        print("------double setup optiizatoin 1-------")
        self.trueParamSimulation()
        self.makeSimulationDir(1)
        
        lbounds = [MIN_ETA, MIN_N, MIN_SIGMA_Y]
        ubounds = [MAX_ETA, MAX_N, MAX_SIGMA_Y]

        cmaes = CMAES_PWQB_LS(7,10,1.0,["MAXITER"],100,"example1.txt", lbounds, ubounds)

        opt3d = Optimizer3D(self.root_folder, self.sh_dir_path+"/setting_simulation_3d.sh",self.ref_dir_paths,\
            self.sim_dir_paths_list, self.particle_skinner_path, self.sim_render_path,self.rho)

        init_m = Param( (lbounds[0] + ubounds[0]) / 2.0, (lbounds[1] + ubounds[1]) / 2.0, (lbounds[2] + ubounds[2]) / 2.0 )
        m_breve = opt3d.optimize(cmaes,init_m,self.setups)
        m_breve_1 = opt3d.param_data[np.argmin(opt3d.loss_data)]

        mechanism = Mechanism3D()
        self.setups = mechanism.searchNewSetup_orthognality_for_second_setup(m_breve, self.setups)
        self.exportSecondSetupXML(m_breve)

def main():
    xmlData = xmlParser.OptXMLData(args.first_dir + "/settings.xml")
    optimizer = Main(xmlData)
    optimizer.optimization_1()

if __name__ == '__main__':
    main()
