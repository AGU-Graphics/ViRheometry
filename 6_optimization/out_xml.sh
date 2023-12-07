#!/bin/bash

OUT_DIR=${1}
OUT_FNAME=${2}
ROOT_DIR_PATH=${3}

DAM_H0=${4}
DAM_W0=${5}
RHO=${6}

ETA=${7}
N=${8}
SIGMA_Y=${9}

cat <<EOF > ${OUT_DIR}"/"${OUT_FNAME}
<?xml version="1.0"?>
<Optimizer>
  <path
    root_dir_path="${ROOT_DIR_PATH}"
    GL_render_path="../libs/3D/GLRender3d/build/GLRender3d"
    mpm_path="../libs/3D/MPM3d/AGTaichiMPM.py"
    particle_skinner_path="../libs/ParticleSkinner3DTaichi.py"
    shell_script_dir_path="../libs/3D/shellScript3d"
    GL_emulation_render_path="../libs/3D/GLEmulationRender3d/build/GLEmulationRender3d"
  />

  <setup
    RHO="${RHO}"
    H="${DAM_H0}"
    W="${DAM_W0}"
  />

  <initial_material
    eta="${ETA}"
    n="${N}"
    sigmaY="${SIGMA_Y}"
  />
  <cuboid min="-0.150000 -0.150000 -0.150000" max="${DAM_W0} ${DAM_H0} 4.150000" density="${RHO}" cell_samples_per_dim="2" vel="0.0 0.0 0.0" omega="0.0 0.0 0.0" />
  <static_box min="-100.000000 -1.000000 -100.000000" max="100.000000 0.000000 100.000000" boundary_behavior="sticking"/>
  <static_box min="-1.000000 0.000000 0.000000" max="0.000000 20.000000 4.000000" boundary_behavior="sticking"/>
  <static_box min="-1.000000 0.000000 -0.300000" max="${DAM_W0} 20.000000 0.000000" boundary_behavior="sticking"/>
  <static_box min="-1.000000 0.000000 4.000000" max="${DAM_W0} 20.000000 4.300000" boundary_behavior="sticking"/>

</Optimizer>
EOF

