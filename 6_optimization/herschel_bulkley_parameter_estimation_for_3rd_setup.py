import sys
sys.path.append("../libs")
import os
import argparse
import numpy as np
import xmlParser
import time
import const
import copy
from optimizer import Optimizer3D
from cmaes import CMAES_PWQB_LS
from mechanism import Mechanism, Mechanism3D
from setup import Setup
from param import Param
from const import CGS
const = CGS
MAX_ETA, MAX_N, MAX_SIGMA_Y, MIN_ETA, MIN_N, MIN_SIGMA_Y = const.MAX_ETA, const.MAX_N, const.MAX_SIGMA_Y, const.MIN_ETA, const.MIN_N, const.MIN_SIGMA_Y
MAX_LOSS_THR, MIN_H, MAX_H, MIN_W, MAX_W = const.MAX_LOSS_THR, const.MIN_H, const.MAX_H, const.MIN_W, const.MAX_W
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--first_dir')
parser.add_argument('-s', '--second_dir')
parser.add_argument('-t', '--third_dir')
args = parser.parse_args()
from_dirs = [args.first_dir, args.second_dir, args.third_dir]

extent_eta = const.extent_eta
extent_n = const.extent_n
extent_sigmaY = const.extent_sigmaY
max_frame = 8

class Main:
    def __init__(self, xmlDatas):
        self.root_folder = os.path.abspath(xmlDatas[-1].pathData.root_dir_path)
        self.sh_dir_path = os.path.abspath(xmlDatas[-1].pathData.shell_script_dir_path)
        self.particle_skinner_path = os.path.abspath(xmlDatas[-1].pathData.particle_skinner_path)
        self.mpm_path = os.path.abspath(xmlDatas[-1].pathData.mpm_path)
        self.sim_render_path = os.path.abspath(xmlDatas[-1].pathData.GL_render_path)
        self.rho = xmlDatas[-1].setupData.rho
        eta = xmlDatas[-1].initialMaterialData.eta
        n = xmlDatas[-1].initialMaterialData.n
        sigmaY = xmlDatas[-1].initialMaterialData.sigmaY
        self.initial_material = Param(eta, n, sigmaY)
        
        self.setups = []
        H = xmlDatas[0].setupData.H
        W = xmlDatas[0].setupData.W
        self.setups.append(Setup(H,W,0.5))

        H2 = xmlDatas[1].setupData.H
        W2 = xmlDatas[1].setupData.W
        self.setups.append(Setup(H2,W2,0.5))
    
        H3 = xmlDatas[1].setupData.H
        W3 = xmlDatas[1].setupData.W
        self.setups.append(Setup(H3, W3, 0.5))


        self.setups[0].display_status()
        self.setups[1].display_status()

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
                elif yorN in ['n', 'no', '']:
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

    def trueParamSimulation(self, i):
        self.makeReferenceDir(i)

        from_dir = from_dirs[i - 1]
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

    def exportXML(self, m_breve):
        material_name = os.path.basename(self.root_folder).split(",")[0].split("_")[0]
        s = self.setups[-1]
        Hcm = str(round(s.H, 2))
        Wcm = str(round(s.W, 2))
        arg_1 = "../data/ref_" + material_name + "_" + Hcm + "_" + Wcm + "_4/"
        if not os.path.exists(arg_1):
            os.mkdir(arg_1)
        arg_2 = "settings.xml"
        arg_3 = "../data/" + material_name + "_4"
        arg_4 = Hcm
        arg_5 = Wcm
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
        explanation = "Please place the experiment frames from config01.png to config08.png based on the initial setup (W, H) = (" + Hcm + ", " + Wcm + ") ,\n along with the calibration file camera_params.xml, into the specified folder -> " + arg_1 
        print(explanation)

    def optimization_3(self):

        m_breve = self.initial_material
        print("initial material: ")
        m_breve.display_status()

        self.trueParamSimulation(1)
        self.makeSimulationDir(1)

        self.trueParamSimulation(2)
        self.makeSimulationDir(2)

        self.trueParamSimulation(3)
        self.makeSimulationDir(3)

        lbounds = [MIN_ETA, MIN_N, MIN_SIGMA_Y]
        ubounds = [MAX_ETA, MAX_N, MAX_SIGMA_Y]

        cmaes = CMAES_PWQB_LS(7,10,(2.0/3) ** 2,["MAXITER"],100, "example3.txt", lbounds, ubounds)
        opt3d = Optimizer3D(self.root_folder, self.sh_dir_path+"/setting_simulation_3d.sh",self.ref_dir_paths,\
            self.sim_dir_paths_list, self.particle_skinner_path, self.sim_render_path,self.rho)
        m_breve = opt3d.optimize(cmaes, m_breve, self.setups)
        print("m_breve3 : ", m_breve.eta,", ", m_breve.n,", ", m_breve.sigmaY)
        for s in self.setups:
            print("setup: ", "H=", s.H, "W=", s.W, "weight=", s.weight)
        mechanism = Mechanism3D()
        self.setups = mechanism.searchNewSetup_orthognality_for_quadruple_setup(m_breve, self.setups)
        self.exportXML(m_breve)

def main():
    xmlData1 = xmlParser.OptXMLData(args.first_dir + "/settings.xml")
    xmlData2 = xmlParser.OptXMLData(args.second_dir + "/settings.xml")
    xmlData3 = xmlParser.OptXMLData(args.third_dir + "/settings.xml")
    xmlDatas = [xmlData1, xmlData2, xmlData3]
    optimizer = Main(xmlDatas)
    optimizer.optimization_3() 

if __name__ == '__main__':
    main()
