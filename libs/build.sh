#!/bin/zsh
env=${1}
ENV1="ubuntu"
ENV2="mac"

if [ ${env} = ${ENV2} ]; then
    rm -r ./libs/3D
    cp -r libs/3D_mac ./libs/3D
elif [ ${env} = ${ENV1} ]; then
    rm -r ./libs/3D
    cp -r libs/3D_ubuntu ./libs/3D
fi

cp -r ./libs/MPM3d ./libs/3D/MPM3d

mkdir ./libs/3D/GLRender3d/build
cd ./libs/3D/GLRender3d/build
cmake ..
make
cd ../../../../

mkdir ./libs/3D/ParticleSkinner3DTaichi/cpp_marching_cubes/build
cd ./libs/3D/ParticleSkinner3DTaichi/cpp_marching_cubes/build
cmake ..
make
cd ../../../../../

self_dir=$(cd $(dirname $0); pwd)

ln -s $self_dir/3D/ParticleSkinner3DTaichi/ParticleSkinner3DTaichi.py $self_dir/ParticleSkinner3DTaichi.py
ln -s $self_dir/3D/ParticleSkinner3DTaichi/cpp_marching_cubes/build/cpp_marching_cubes $self_dir/cpp_marching_cubes