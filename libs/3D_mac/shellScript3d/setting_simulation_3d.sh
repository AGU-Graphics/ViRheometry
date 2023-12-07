#!/bin/bash
OUTPATH=${1}
OUTDIR=${1}"/setting_"${2}".xml"

ETA=${3}
HPOWER=${4}
SIGMA_Y=${5}

DAM_H0=${6}
DAM_W0=${7}

SKINNER_PATH=${8}
GLRENDER_PATH=${9}

REF_DIR_PATH=${10}
SIM_DIR_PATH=${11}
DENSITY=${12}

#PLANE_MAX_X=`echo "scale=5; ${SIMBOX_X} - ${SIM_BUFFA}" | bc -l`

echo "OutDir:"${OUTDIR}

echo "eta="${ETA}
echo "herschel_bulkley_power="${HPOWER}
echo "yield_stress="${SIGMA_Y}

echo "DAM_H0="${DAM_H0}
echo "DAM_W0="${DAM_W0}

cat <<EOF > ${OUTDIR}
<?xml version="1.0"?>
<AGTaichiMPM3D>
    
  <path
    ref_dir_path="${REF_DIR_PATH}"
    sim_dir_path="${SIM_DIR_PATH}"
  />

  <emulation
    emulation="0"
    emu_dir_path=""
  />

  <particle_skinner
    path="${SKINNER_PATH}"
    grid_space="0.063"
    file_type="obj"
  />
  
  <GLRender path="${GLRENDER_PATH}">
    <camera 
      eyepos="8.45696138 14.32447243 25.97375739"  
      quat="0.96642661 -0.25693407 0.00158011 -0.00141509" 
      window_size="960 540" 
      fov="36.17541"
    />
  </GLRender>
 
  <integrator
    herschel_bulkley_power = "${HPOWER}"
    eta = "${ETA}"
    yield_stress = "${SIGMA_Y}"
  />
  
<cuboid min="-0.150000 -0.150000 -0.150000" max="${DAM_W0} ${DAM_H0} 4.300000" density="${DENSITY}" cell_samples_per_dim="2" vel="0.0 0.0 0.0" omega="0.0 0.0 0.0" />
<static_box min="-100.000000 -1.000000 -100.000000" max="100.000000 0.000000 100.000000" boundary_behavior="sticking"/>
<static_box min="-1.000000 0.000000 0.000000" max="0.000000 20.000000 4.150000" boundary_behavior="sticking"/>
<static_box min="-1.000000 0.000000 -0.300000" max="${DAM_W0} 20.000000 0.000000" boundary_behavior="sticking"/>
<static_box min="-1.000000 0.000000 4.150000" max="${DAM_W0} 20.000000 4.450000" boundary_behavior="sticking"/>
  
</AGTaichiMPM3D>
EOF

