import argparse
import os
import random
from datetime import datetime
random.seed(datetime.now().timestamp())

H = random.randint(20, 70)
W = random.randint(20, 70)
Hcm = str(round(H * 0.1, 2))
Wcm = str(round(W * 0.1, 2))
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--density')
parser.add_argument('-m', '--material_name')

args = parser.parse_args()

arg_1 = "../data/ref_" + args.material_name + "_" + Hcm + "_" + Wcm + "_1/"
if not os.path.exists(arg_1):
    os.mkdir(arg_1)
arg_2 = "settings.xml"
arg_3 = "../data/" + args.material_name + "_1"
arg_4 = str(round(H * 0.1, 2))
arg_5 = str(round(W * 0.1, 2))
arg_6 = args.density

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

os.system(command)
explanation = "Before running the optimization for the 1st experiment, please put the following files into the specified folder \"" + arg_1 + "\": \n"
explanation += "  -  the binary frames (config01.png to config08.png) from the experiment with the setting W = " + arg_5 + ", and H = " + arg_4 + "\n"
explanation += "  -  the camera parameter setting (camera_params.xml) from the calibration"
print(explanation)
