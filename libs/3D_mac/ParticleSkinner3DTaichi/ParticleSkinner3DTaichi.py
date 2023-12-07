# This code is based on the code https://github.com/haimasree/particleskinner by Haimasree Bhattacharya and Adam Bargteil.
# We have done a few modifications for the purpose of optimizing Herschel-Bulkley parameters of complex fluids:
# - The redistance operation is omitted for performance
# - The smoothing is only performed near the 0 level set
# These modifications may not be valid for general skinning purposes.

import ctypes
from ctypes import *
import numpy as np
import sys
import taichi as ti
import time as tm
import os
import subprocess
import math

# ti.init(arch=ti.cuda, default_fp=real)
real = ti.f32
realnp = np.float32
mem = 3
ti.init(default_fp=real, arch=ti.cpu, device_memory_GB=mem)
LARGE_FLOAT = 1.0e10
LARGE_FLOAT_THR = 1.0e9

@ti.data_oriented
class PartickeSkinner3DTaichi:
    def __init__(self, h, config_fn, phi_fn, obj_fn, marching_cube_path):
        self.py_h = h
        self.py_config_fn = config_fn
        self.py_phi_fn = phi_fn
        self.py_obj_fn = obj_fn
        self.py_marching_cube_path = marching_cube_path

        f = open(self.py_config_fn, "rb")
        self.py_num_points = int.from_bytes(f.read(4), 'little')
        print("num_points: ", self.py_num_points)

        self.py_positions = np.frombuffer(f.read(3*self.py_num_points*4), dtype=realnp).reshape([self.py_num_points, 3], order='C')
        self.py_radii = np.frombuffer(f.read(self.py_num_points*4), dtype=realnp)
        self.py_velocities = np.frombuffer(f.read(3*self.py_num_points*4), dtype=realnp).reshape([self.py_num_points, 3], order='C')

        self.py_r_min = 0.5 * math.sqrt(3.0) * h
        self.py_r_max = 4.0 * self.py_r_min
        self.py_r_init = 0.5 * (self.py_r_min + self.py_r_max)
        self.py_dtLaplace = 0.1 * h * h
        dtBiharmonicGain = 1.0
        self.py_dtBiharmonic = 0.01*dtBiharmonicGain*h*h*h*h
        self.py_max_stretch = 4.0

        print("particle data of the first 5 particles: ")
        for i in range(5):
            print("idx ", i, ": p =", self.py_positions[i, :], ", r =", self.py_radii[i], ", vel =", self.py_velocities[i, :])

        f.close()

        self.py_bb_min = self.py_positions.min(axis=0)
        self.py_bb_max = self.py_positions.max(axis=0)
        print("bb: ", self.py_bb_min, " - ", self.py_bb_max)

        #self.py_bb_min -= 5*self.py_h + 1.0*self.py_r_max + 10*self.py_h
        self.py_bb_min -= 5*self.py_h + 1.0*self.py_r_max + 15*self.py_h
        #self.py_bb_max += 5*self.py_h + 1.0*self.py_r_max + 10*self.py_h
        self.py_bb_max += 5*self.py_h + 1.0*self.py_r_max + 15*self.py_h

        n_min = np.floor(self.py_bb_min/self.py_h).astype(int)
        n_max = np.ceil(self.py_bb_max/self.py_h).astype(int)

        self.py_n = n_max - n_min
        self.py_bb_min = n_min.astype(float) * self.py_h
        self.py_bb_max = n_max.astype(float) * self.py_h
        print("resized bb: ", self.py_bb_min, " - ", self.py_bb_max)
        print("bb resolution: ", self.py_n)

        self.ti_num_points = ti.field(dtype=ti.i32, shape=())
        self.ti_num_points[None] = self.py_num_points
        self.ti_positions = ti.Vector.field(3, dtype=real, shape=self.ti_num_points[None])
        self.ti_positions.from_numpy(self.py_positions)
        self.ti_radii = ti.field(dtype=real, shape=self.ti_num_points[None])
        self.ti_radii.from_numpy(self.py_radii)
        self.ti_velocities = ti.Vector.field(3, dtype=real, shape=self.ti_num_points[None])
        self.ti_velocities.from_numpy(self.py_velocities)

        test_pos = self.ti_positions.to_numpy()
        print("test_pos[0, 0]: ", test_pos[0, 0])

        self.ti_phis = [ti.field(dtype=real, shape=self.py_n), ti.field(dtype=real, shape=self.py_n)]
        self.py_phi_current = 0
        self.ti_phi_min = ti.field(dtype=real, shape=self.py_n)
        self.ti_phi_max = ti.field(dtype=real, shape=self.py_n)
        self.ti_accepted = ti.field(dtype=int, shape=self.py_n)
        self.ti_laplacian = ti.field(dtype=real, shape=self.py_n)
        self.ti_biharmonic = ti.field(dtype=real, shape=self.py_n)

    @ti.kernel
    def rasterize(self):
        ti_bb_min = ti.Vector(self.py_bb_min)

        for I in ti.grouped(self.ti_phis[self.py_phi_current]):
            self.ti_phis[self.py_phi_current][I] = LARGE_FLOAT

        #print("num_pts: ", self.ti_num_points[None])

        #width = int(ti.ceil(self.py_r_max / self.py_h) + 1)
        width = int(ti.ceil(self.py_r_max / self.py_h) + 8)
        for p in range(self.ti_num_points[None]):
        #for p in range(5):
            bin = ti.cast(ti.floor((self.ti_positions[p] - ti_bb_min) / self.py_h), ti.i32);
            base = bin - width

            #print("pos: ", self.ti_positions[p])
            #print("bin: ", bin)
            #print("base: ", base)

            for i, j, k in ti.ndrange(width*2, width*2, width*2):
                offset = ti.Vector([i, j, k])
                gp = ti_bb_min + (base + offset).cast(real) * self.py_h
                d = (gp - self.ti_positions[p]).norm()
                ti.atomic_min(self.ti_phis[self.py_phi_current][base+offset], d)
                # print("gp: ", gp)
                # print("d: ", d)

        for I in ti.grouped(self.ti_phis[self.py_phi_current]):
            if self.ti_phis[self.py_phi_current][I] > LARGE_FLOAT_THR:
                self.ti_phi_min[I] = LARGE_FLOAT
                self.ti_phi_max[I] = LARGE_FLOAT
            else:
                d = self.ti_phis[self.py_phi_current][I]
                self.ti_phis[self.py_phi_current][I] = d - self.py_r_init
                self.ti_phi_min[I] = d - self.py_r_min
                self.ti_phi_max[I] = d - self.py_r_max


    @ti.func
    def cdX(self, i: int, j: int, k: int, phi: ti.template()):
        return (phi[i+1,j,k] - phi[i-1,j,k]) / (2*self.py_h)

    @ti.func
    def cdY(self, i: int, j: int, k: int, phi: ti.template()):
        return (phi[i,j+1,k] - phi[i,j-1,k]) / (2*self.py_h)

    @ti.func
    def cdZ(self, i: int, j: int, k: int, phi: ti.template()):
        return (phi[i,j,k+1] - phi[i,j,k-1]) / (2*self.py_h)

    @ti.kernel
    def computeLaplacian(self, phi: ti.template()):
        divisor = self.py_h * self.py_h
        updateBand = 4 * self.py_h
        for _i, _j, _k in ti.ndrange(phi.shape[0]-2, phi.shape[1]-2, phi.shape[2]-2):
            i = _i + 1
            j = _j + 1
            k = _k + 1
            if ti.abs(phi[i, j, k]) <= updateBand:
                self.ti_laplacian[i, j, k] = (phi[i+1, j, k] + phi[i-1, j, k] + phi[i, j+1, k] + phi[i, j-1, k] + phi[i, j, k+1] + phi[i, j, k-1] - 6 * phi[i, j, k]) / divisor

    @ti.kernel
    def stepLaplacian(self, dt: real, phi_from: ti.template(), phi_to: ti.template()):
        #change = 0.0
        updateBand = 4 * self.py_h

        for _i, _j, _k in ti.ndrange(phi_from.shape[0]-4, phi_from.shape[1]-4, phi_from.shape[2]-4):
            i = _i + 2
            j = _j + 2
            k = _k + 2
            if ti.abs(phi_from[i, j, k]) <= updateBand:
                phix = self.cdX(i, j, k, phi_from)
                phiy = self.cdY(i, j, k, phi_from)
                phiz = self.cdZ(i, j, k, phi_from)
                val = self.ti_laplacian[i, j, k]
                gradMag = ti.sqrt(phix*phix + phiy*phiy + phiz*phiz)
                updatedPhi = phi_from[i, j, k] + val * dt * gradMag;
                updatedPhi = ti.min(updatedPhi, self.ti_phi_min[i, j, k])
                updatedPhi = ti.max(updatedPhi, self.ti_phi_max[i, j, k])
                phi_to[i, j, k] = updatedPhi
                #change += ti.abs(val)

        #print("change in this step: ", change)

    def doLaplacianSmoothing(self, iter, dt):
        for i in range(iter):
            self.computeLaplacian(self.ti_phis[self.py_phi_current])
            self.stepLaplacian(dt, self.ti_phis[self.py_phi_current], self.ti_phis[1-self.py_phi_current])
            self.py_phi_current = 1 - self.py_phi_current

    @ti.kernel
    def computeBiharmonic(self, phi: ti.template()):
        divider = self.py_h * self.py_h
        updateBand = 3 * self.py_h

        for _i, _j, _k in ti.ndrange(phi.shape[0]-4, phi.shape[1]-4, phi.shape[2]-4):
            i = _i + 2
            j = _j + 2
            k = _k + 2

            if ti.abs(phi[i, j, k]) <= updateBand:
                self.ti_biharmonic[i, j, k] = (self.ti_laplacian[i+1, j, k] + self.ti_laplacian[i-1, j, k] + self.ti_laplacian[i, j+1, k] + self.ti_laplacian[i, j-1, k] + self.ti_laplacian[i, j, k+1] + self.ti_laplacian[i, j, k-1] - 6 * self.ti_laplacian[i, j, k]) / divider

    @ti.kernel
    def stepBiharmonic(self, dt: real, phi_from: ti.template(), phi_to: ti.template()):
        #change = 0.0
        updateBand = 3 * self.py_h

        for _i, _j, _k in ti.ndrange(phi_from.shape[0]-4, phi_from.shape[1]-4, phi_from.shape[2]-4):
            i = _i + 2
            j = _j + 2
            k = _k + 2

            if ti.abs(phi_from[i, j, k]) <= updateBand:
                phix = self.cdX(i, j, k, phi_from)
                phiy = self.cdY(i, j, k, phi_from)
                phiz = self.cdZ(i, j, k, phi_from)
                gradMag = ti.sqrt(phix * phix + phiy * phiy + phiz * phiz)
                val = self.ti_biharmonic[i, j, k]
                updatedPhi = phi_from[i, j, k] - val * dt * gradMag
                updatedPhi = ti.min(updatedPhi, self.ti_phi_min[i, j, k])
                updatedPhi = ti.max(updatedPhi, self.ti_phi_max[i, j, k])
                phi_to[i, j, k] = updatedPhi
                #change += ti.abs(val)

        #print("change in this step: ", change)


    def doBiharmonicSmoothing(self, iter, dt):
        for i in range(iter):
            self.computeLaplacian(self.ti_phis[self.py_phi_current])
            self.computeBiharmonic(self.ti_phis[self.py_phi_current])
            self.stepBiharmonic(dt, self.ti_phis[self.py_phi_current], self.ti_phis[1-self.py_phi_current])
            self.py_phi_current = 1 - self.py_phi_current

    def process(self):
        self.py_phi_current = 0
        print("rasterize...")
        self.rasterize()
        #print("redistance...")
        #self.redistance( self.ti_phi_min, self.ti_phi_temp, self.ti_accepted )
        # self.sweep_and_set( self.ti_phi_temp, self.ti_accepted, self.ti_phi_min )
        #self.redistance( self.ti_phi_max, self.ti_phi_temp, self.ti_accepted )
        # self.sweep_and_set( self.ti_phi_temp, self.ti_accepted, self.ti_phi_max )
        #self.redistance( self.ti_phi, self.ti_phi_temp, self.ti_accepted )
        # self.sweep_and_set( self.ti_phi_temp, self.ti_accepted, self.ti_phi )

        print("doLaplacianSmoothing...")
        self.doLaplacianSmoothing(15, self.py_dtLaplace)
        print("doBiharmonicSmoothing...")
        self.doBiharmonicSmoothing(500, self.py_dtBiharmonic)
        print("process done...")

    def postprocess(self):
        bb_min = np.array(self.py_bb_min, dtype='float32')
        b_bb_min = bb_min.tobytes()
        h = np.array(self.py_h, dtype='float32')
        b_h = h.tobytes()
        resolution = np.array(self.py_n, dtype='int32')
        b_res = resolution.tobytes()
        phi = self.ti_phis[self.py_phi_current].to_numpy().astype('float32')
        b_phi = phi.tobytes('F')

        print("phi.shape: ", phi.shape)

        print("phi: ", np.min(phi), " - ", np.max(phi))

        f = open(self.py_phi_fn, "wb")
        f.write(b_bb_min)
        f.write(b_h)
        f.write(b_res)
        f.write(b_phi)
        f.close()

        subprocess.call([self.py_marching_cube_path, self.py_phi_fn, self.py_obj_fn])


if len(sys.argv) < 6:
    print("usage: python ParticleSkinner3DTaichi.py <[input] h> <[input] mpm dat file> <[output] phi file> <[output] obj file> <marching cube path>")
    exit(-1)

particleSkinner3D = PartickeSkinner3DTaichi(float(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
print("processing...")
particleSkinner3D.process()
particleSkinner3D.postprocess()
print("done.")
