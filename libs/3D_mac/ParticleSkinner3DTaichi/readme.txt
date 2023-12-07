How to use ParticleSkinner3DTaichi

* Preparation
1. compile cpp_marching_cubes
2. At the working directory, set a symbolic link to ParticleSkinner3DTaichi.py のショートカットをおく
3. At the same working directory, set a symbolic link to ./cpp_marching_cubes 
4. Put mpm particle data (e.g., config_01.dat) into the working directory
5. run ParticleSkinner3DTaichi for each particle data via ParticleSkinner3DTaichi.py config_01.dat config_01_phi.dat 

* Note
1. The rasterize procedure only does the computation at a narrow band around the particles. The outside will have a positive sign, while the inner region will have a negative sign
2. redistance operation is omitted for performance
3. The smoothing is only performed near the 0 level set
4. marching cubes is run via the c++ code