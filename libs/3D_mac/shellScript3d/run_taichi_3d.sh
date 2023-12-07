outpath=${1}
outfile=${1}"/setting.xml"
taichi_script=${2}
tmp=/tmp/$$exec_taichi

### Simulation Settings ###
dt="0.000075"
bulk_modulus="100000.0"
shear_modulus="10000.0"
flip_pic_alpha="0.95"

fps="24"
max_frames="8"
emulation="0" # -f

#grid settings
grid_min="-1.0 -1.0 -10.0"
grid_max="30.0 8.0 14.0"
cell_width="0.126"


#cuboid settings
cuboid_min="-0.15 -0.15 -0.15"
cuboid_max="7.000000 7.000000 4.300000"
density=${3}
cell_samples_per_dim="2"
vel="0.0 0.0 0.0"

#gravity setting
gravity="-981.0"

{
cat <<EOF > ${outfile}
<?xml version="1.0"?>
<AGTaichiMPM3D>

  <integrator
    dt = "${dt}"
    bulk_modulus = "${bulk_modulus}"
    shear_modulus = "${shear_modulus}"
    flip_pic_alpha = "${flip_pic_alpha}"
    max_frames = "${max_frames}"
    fps = "${fps}"
    optimization_emulation = "${emulation}"
  />
  <grid min="${grid_min}" max="${grid_max}" cell_width="${cell_width}"/>
  <cuboid min="${cuboid_min}" max="${cuboid_max}" density="${density}" cell_samples_per_dim="${cell_samples_per_dim}" vel="${vel}"/>
  <near_earth_gravity f="0.0 ${gravity} 0.0"/>
</AGTaichiMPM3D>
EOF

} &&{
cat << EOF > ${tmp}
#!/bin/bash
python3 "${taichi_script}" "${outpath}" "${outfile}"
EOF
chmod +x ${tmp}

} &&{
open ${tmp}

} &&{
echo "taichi is activated."

} ||{
echo "Activating taichi is Failed."
}
