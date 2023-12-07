import ctypes
import numpy as np
import subprocess
from subprocess import PIPE
import sys
import initXmlParser
import xmlParser
import setUpXmlParser
import setUpXmlParser2
import taichi as ti
from ctypes import *
import glob, os, ntpath, shutil
import time as tm
import os
from enum import Enum
from pathlib import Path
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN, ROUND_FLOOR

real = ti.f32
realnp = np.float32
ti.init(default_fp=real, arch=ti.gpu)

class Mode(Enum):
    IDLE = 1
    PROCESSING = 2

count = 0
mode = Mode.IDLE


@ti.data_oriented
class AGTaichiMPM:
    def __init__(self, xmlData):
        self.ti_hb_n = ti.field(dtype=real, shape=())
        self.ti_hb_eta = ti.field(dtype=real, shape=())
        self.ti_hb_sigmaY = ti.field(dtype=real, shape=())

        self.ti_num_boxes = ti.field(dtype=int, shape=())
        self.ti_num_boxes[None] = 4
        self.ti_static_box_min = ti.Vector.field(3, dtype=real, shape=self.ti_num_boxes[None])
        self.ti_static_box_max = ti.Vector.field(3, dtype=real, shape=self.ti_num_boxes[None])
        self.ti_static_box_type = ti.field(ti.i32, shape=self.ti_num_boxes[None])

        self.py_camera_position = ""
        self.py_camera_quat = ""
        self.py_camera_window = ""
        self.py_camera_fov = ""

        self.py_ref_dir_path = ""
        self.py_sim_dir_path = ""

        self.py_emulation = 0
        self.py_emu_dir_path = ""

        self.py_particleSkinnerApp = ""
        self.py_skin_size = 0.05
        self.py_file_type = "obj"
        self.py_GLRenderApp = ""

        self.py_fps = xmlData.integratorData.fps
        self.py_optimization_emulation = xmlData.integratorData.optimization_emulation
        # material parameters
        self.py_kappa     = xmlData.integratorData.bulk_modulus
        self.py_mu        = xmlData.integratorData.shear_modulus

        # flip-pic alpha
        self.py_alpha = xmlData.integratorData.flip_pic_alpha

        # temporal/spatial resolution
        self.py_dt = xmlData.integratorData.dt
        self.py_dx = xmlData.gridData.cell_width
        self.py_invdx = 1.0 / self.py_dx

        # near earth gravity
        self.ti_g = ti.Vector(xmlData.nearEarthGravityData.g)

        # iteration count
        self.ti_iteration = ti.field(dtype=int, shape=())
        self.ti_iteration[None] = 0

        # max time
        self.py_max_frames = xmlData.integratorData.max_frames
        self.py_num_saved_frames = 1

        # configuring grid by using the specified grid center and cell width as is
        # min and max will be recomputed because the specified grid size may not agree with the specified cell width

        # compute grid center and tentative grid width
        grid_center = (xmlData.gridData.max + xmlData.gridData.min) * 0.5
        grid_width = xmlData.gridData.max - xmlData.gridData.min
        self.py_cell_count = np.ceil(grid_width / self.py_dx).astype(int)

        # recompute grid width, min and max
        grid_width = self.py_cell_count.astype(realnp) * self.py_dx
        self.ti_grid_min = ti.Vector(grid_center - 0.5 * grid_width)
        self.ti_grid_max = ti.Vector(grid_center + 0.5 * grid_width)

        # allocating fields for grid mass and velocity (momentum)

        self.ti_grid_m = ti.field(dtype=real, shape=self.py_cell_count)
        self.ti_grid_x = ti.Vector.field(3, dtype=real, shape=self.py_cell_count)
        self.ti_grid_v = ti.Vector.field(3, dtype=real, shape=self.py_cell_count)
        self.ti_grid_a = ti.Vector.field(3, dtype=real, shape=self.py_cell_count)
        # for debug
        self.ti_grid_pos = ti.Vector.field(3, dtype=real, shape=np.prod(self.py_cell_count))
        self.ti_grid_color = ti.field(ti.i32, shape=np.prod(self.py_cell_count))

        #self.ti_particle_init_min = ti.Vector(xmlData.cuboidData.min)
        self.ti_particle_init_min = ti.Vector.field(3, dtype=real, shape=1)
        self.ti_particle_init_min.from_numpy(xmlData.cuboidData.min.astype(realnp).reshape(1,3))
        self.py_particle_init_cell_samples_per_dim = xmlData.cuboidData.cell_samples_per_dim
        self.ti_particle_init_vel = ti.Vector(xmlData.cuboidData.vel)

        self.py_particle_hl = 0.5 * self.py_dx / xmlData.cuboidData.cell_samples_per_dim
        print('py_particle_hl: ', self.py_particle_hl)

        self.py_particle_volume = (self.py_dx / xmlData.cuboidData.cell_samples_per_dim)**3
        self.py_particle_mass = xmlData.cuboidData.density * self.py_particle_volume

        # in
        # itialize max number of particles
        cuboid_width = xmlData.cuboidData.max - xmlData.cuboidData.min
        self.ti_particle_ndcount = ti.field(dtype=int, shape=3)
        self.ti_particle_ndcount.from_numpy(np.ceil(cuboid_width * xmlData.cuboidData.cell_samples_per_dim / self.py_dx).astype(np.int32))
        print('cuboid_width: ', cuboid_width)
        print('xmlData.cuboidData.cell_samples_per_dim: ', xmlData.cuboidData.cell_samples_per_dim)
        print('py_dx: ', self.py_dx)

        self.ti_particle_count = ti.field(dtype=int, shape=())
        self.ti_particle_count[None] = np.prod(self.ti_particle_ndcount.to_numpy())
        print('ti_particle_count: ', self.ti_particle_count[None])
        self.ti_particle_is_inner_of_box = ti.field(int, shape=self.ti_particle_count[None])
        self.ti_particle_x = ti.Vector.field(3, dtype=real, shape=self.ti_particle_count[None])
        self.ti_particle_v = ti.Vector.field(3, dtype=real, shape=self.ti_particle_count[None])
        self.ti_particle_be = ti.Matrix.field(3, 3, dtype=real, shape=self.ti_particle_count[None])
        self.ti_particle_C = ti.Matrix.field(3, 3, dtype=real, shape=self.ti_particle_count[None])
        # for debug
        self.ti_particle_color_f = ti.field(real, shape=self.ti_particle_count[None])
        self.ti_particle_color = ti.field(ti.i32, shape=self.ti_particle_count[None])

    @ti.kernel
    def changeSetUpDataKernel_Box0(self, box_0_min_x: ti.f32, box_0_min_y: ti.f32, box_0_min_z: ti.f32, box_0_max_x: ti.f32, box_0_max_y: ti.f32, box_0_max_z: ti.f32, isSticky: ti.i32):
        self.ti_static_box_min[0][0] = box_0_min_x
        self.ti_static_box_min[0][1] = box_0_min_y
        self.ti_static_box_min[0][2] = box_0_min_z
        self.ti_static_box_max[0][0] = box_0_max_x
        self.ti_static_box_max[0][1] = box_0_max_y
        self.ti_static_box_max[0][2] = box_0_max_z
        self.ti_static_box_type[0] = isSticky

    @ti.kernel
    def changeSetUpDataKernel_Box1(self, box_1_min_x: ti.f32, box_1_min_y: ti.f32, box_1_min_z: ti.f32, box_1_max_x: ti.f32, box_1_max_y: ti.f32, box_1_max_z: ti.f32, isSticky: ti.i32):
        self.ti_static_box_min[1][0] = box_1_min_x
        self.ti_static_box_min[1][1] = box_1_min_y
        self.ti_static_box_min[1][2] = box_1_min_z
        self.ti_static_box_max[1][0] = box_1_max_x
        self.ti_static_box_max[1][1] = box_1_max_y
        self.ti_static_box_max[1][2] = box_1_max_z
        self.ti_static_box_type[1] = isSticky

    @ti.kernel
    def changeSetUpDataKernel_Box2(self, box_2_min_x: ti.f32, box_2_min_y: ti.f32, box_2_min_z: ti.f32, box_2_max_x: ti.f32, box_2_max_y: ti.f32, box_2_max_z: ti.f32, isSticky: ti.i32):
        self.ti_static_box_min[2][0] = box_2_min_x
        self.ti_static_box_min[2][1] = box_2_min_y
        self.ti_static_box_min[2][2] = box_2_min_z
        self.ti_static_box_max[2][0] = box_2_max_x
        self.ti_static_box_max[2][1] = box_2_max_y
        self.ti_static_box_max[2][2] = box_2_max_z
        self.ti_static_box_type[2] = isSticky

    @ti.kernel
    def changeSetUpDataKernel_Box3(self, box_3_min_x: ti.f32, box_3_min_y: ti.f32, box_3_min_z: ti.f32, box_3_max_x: ti.f32, box_3_max_y: ti.f32, box_3_max_z: ti.f32, isSticky: ti.i32):
        self.ti_static_box_min[3][0] = box_3_min_x
        self.ti_static_box_min[3][1] = box_3_min_y
        self.ti_static_box_min[3][2] = box_3_min_z
        self.ti_static_box_max[3][0] = box_3_max_x
        self.ti_static_box_max[3][1] = box_3_max_y
        self.ti_static_box_max[3][2] = box_3_max_z
        self.ti_static_box_type[3] = isSticky

    @ti.kernel
    def changeCuboid(self, init_min_x: ti.f32, init_min_y: ti.f32, init_min_z: ti.f32, nd_count_x: ti.i32, nd_count_y: ti.i32, nd_count_z: ti.i32):
        self.ti_particle_init_min[0][0] = init_min_x
        self.ti_particle_init_min[0][1] = init_min_y
        self.ti_particle_init_min[0][2] = init_min_z

        self.ti_particle_ndcount[0] = nd_count_x
        self.ti_particle_ndcount[1] = nd_count_y
        self.ti_particle_ndcount[2] = nd_count_z
        self.ti_particle_count[None] = nd_count_x * nd_count_y * nd_count_z

        print('changeCuboid - ti_particle_count[None]: ', self.ti_particle_count[None])

    @ti.kernel
    def changeHBParamKernel(self, hb_n: ti.f32, hb_eta: ti.f32, hb_sigma_Y: ti.f32):
        self.ti_hb_n[None]      = hb_n
        self.ti_hb_eta[None]    = hb_eta
        self.ti_hb_sigmaY[None] = hb_sigma_Y

    def changeSetUpData(self, xmlData_camera, xmlData_cuboid):
        xmlData_camera.show()
        eyepos = xmlData_camera.cameraData.eyepos
        quat = xmlData_camera.cameraData.quat
        window_size = xmlData_camera.cameraData.window_size
        fov = xmlData_camera.cameraData.fov
        self.py_camera_position = str(eyepos[0]) + "," + str(eyepos[1])  + "," + str(eyepos[2])
        self.py_camera_quat = str(quat[0]) + "," + str(quat[1]) + "," + str(quat[2]) + "," +str(quat[3])
        self.py_camera_window = str(window_size[0]) + "," + str(window_size[1])
        self.py_camera_fov = str(fov)

        xmlData_cuboid.show()
        cuboid_width = xmlData_cuboid.cuboidData.max - xmlData_cuboid.cuboidData.min
        _nd_count = np.ceil(cuboid_width * xmlData_cuboid.cuboidData.cell_samples_per_dim / self.py_dx).astype(np.int32)
        print('changeSetUpData - cuboid_width: ', cuboid_width)
        print('changeSetUpData - cell_samples_per_dim: ', xmlData_cuboid.cuboidData.cell_samples_per_dim)
        print('changeSetUpData - py_dx: ', self.py_dx)
        self.changeCuboid(xmlData_cuboid.cuboidData.min[0], xmlData_cuboid.cuboidData.min[1], xmlData_cuboid.cuboidData.min[2], _nd_count[0].item(), _nd_count[1].item(), _nd_count[2].item())

        self.ti_num_boxes[None] = len(xmlData.staticBoxList)

        self.changeSetUpDataKernel_Box0(xmlData_cuboid.staticBoxList[0].min[0], xmlData_cuboid.staticBoxList[0].min[1], xmlData_cuboid.staticBoxList[0].min[2], xmlData_cuboid.staticBoxList[0].max[0], xmlData_cuboid.staticBoxList[0].max[1], xmlData_cuboid.staticBoxList[0].max[2], int(xmlData_cuboid.staticBoxList[0].isSticky))
        self.changeSetUpDataKernel_Box1(xmlData_cuboid.staticBoxList[1].min[0], xmlData_cuboid.staticBoxList[1].min[1], xmlData_cuboid.staticBoxList[1].min[2], xmlData_cuboid.staticBoxList[1].max[0], xmlData_cuboid.staticBoxList[1].max[1], xmlData_cuboid.staticBoxList[1].max[2], int(xmlData_cuboid.staticBoxList[1].isSticky))
        self.changeSetUpDataKernel_Box2(xmlData_cuboid.staticBoxList[2].min[0], xmlData_cuboid.staticBoxList[2].min[1], xmlData_cuboid.staticBoxList[2].min[2], xmlData_cuboid.staticBoxList[2].max[0], xmlData_cuboid.staticBoxList[2].max[1], xmlData_cuboid.staticBoxList[2].max[2], int(xmlData_cuboid.staticBoxList[2].isSticky))
        self.changeSetUpDataKernel_Box3(xmlData_cuboid.staticBoxList[3].min[0], xmlData_cuboid.staticBoxList[3].min[1], xmlData_cuboid.staticBoxList[3].min[2], xmlData_cuboid.staticBoxList[3].max[0], xmlData_cuboid.staticBoxList[3].max[1], xmlData_cuboid.staticBoxList[3].max[2], int(xmlData_cuboid.staticBoxList[3].isSticky))

    def changeXMLData(self, xmlData):
        self.py_ref_dir_path = xmlData.pathData.ref_dir_path
        self.py_sim_dir_path = xmlData.pathData.sim_dir_path

        self.py_emulation = xmlData.emulationData.emulation
        self.py_emu_dir_path = xmlData.emulationData.emu_dir_path

        if not self.py_optimization_emulation:
            while True:
                if os.path.isfile(self.py_ref_dir_path+ "/camera_params.xml") and os.path.isfile(self.py_ref_dir_path+ "/settings.xml") :
                    try:
                        setUpData1 = setUpXmlParser.SetUpXMLData(self.py_ref_dir_path+ "/camera_params.xml")
                        setUpData2 = setUpXmlParser2.SetUpXMLData(self.py_ref_dir_path+ "/settings.xml")
                        self.changeSetUpData(setUpData1, setUpData2)
                    except Exception as e:
                        print(e)
                        print("setup.camera_params.xml or settings.xml is something wrong.")
                    else:
                        break
                else:
                    print("cannot find camera_params.xml or settings.xml")
                tm.sleep(1.0)
        else:
            eyepos = xmlData.GLRenderData.cameraData.eyepos
            quat = xmlData.GLRenderData.cameraData.quat
            window_size = xmlData.GLRenderData.cameraData.window_size
            fov = xmlData.GLRenderData.cameraData.fov
            self.py_camera_position = str(eyepos[0]) + "," + str(eyepos[1])  + "," + str(eyepos[2])
            self.py_camera_quat = str(quat[0]) + "," + str(quat[1]) + "," + str(quat[2]) + "," +str(quat[3])
            self.py_camera_window = str(window_size[0]) + "," + str(window_size[1])
            self.py_camera_fov = str(fov)
            cuboid_width = xmlData.cuboidData.max - xmlData.cuboidData.min
            _nd_count = np.ceil(cuboid_width * xmlData.cuboidData.cell_samples_per_dim / self.py_dx).astype(np.int32)
            self.changeCuboid(xmlData.cuboidData.min[0], xmlData.cuboidData.min[1], xmlData.cuboidData.min[2], _nd_count[0].item(), _nd_count[1].item(), _nd_count[2].item())
        
            self.ti_num_boxes[None] = len(xmlData.staticBoxList)
        
            self.changeSetUpDataKernel_Box0(xmlData.staticBoxList[0].min[0], xmlData.staticBoxList[0].min[1], xmlData.staticBoxList[0].min[2], xmlData.staticBoxList[0].max[0], xmlData.staticBoxList[0].max[1], xmlData.staticBoxList[0].max[2], int(xmlData.staticBoxList[0].isSticky))
            self.changeSetUpDataKernel_Box1(xmlData.staticBoxList[1].min[0], xmlData.staticBoxList[1].min[1], xmlData.staticBoxList[1].min[2], xmlData.staticBoxList[1].max[0], xmlData.staticBoxList[1].max[1], xmlData.staticBoxList[1].max[2], int(xmlData.staticBoxList[1].isSticky))
            self.changeSetUpDataKernel_Box2(xmlData.staticBoxList[2].min[0], xmlData.staticBoxList[2].min[1], xmlData.staticBoxList[2].min[2], xmlData.staticBoxList[2].max[0], xmlData.staticBoxList[2].max[1], xmlData.staticBoxList[2].max[2], int(xmlData.staticBoxList[2].isSticky))
            self.changeSetUpDataKernel_Box3(xmlData.staticBoxList[3].min[0], xmlData.staticBoxList[3].min[1], xmlData.staticBoxList[3].min[2], xmlData.staticBoxList[3].max[0], xmlData.staticBoxList[3].max[1], xmlData.staticBoxList[3].max[2], int(xmlData.staticBoxList[3].isSticky))


        self.py_particleSkinnerApp = xmlData.particleSkinnerData.path
        self.py_skin_size = xmlData.particleSkinnerData.grid_space
        self.py_file_type = xmlData.particleSkinnerData.file_type

        self.py_GLRenderApp = xmlData.GLRenderData.path

        self.changeHBParamKernel(xmlData.integratorData.herschel_bulkley_power, xmlData.integratorData.eta, xmlData.integratorData.yield_stress)

    @ti.kernel
    def initialize(self):
        self.ti_iteration[None] = 0
        # clear grid values
        for I in ti.grouped(self.ti_grid_m):
            self.ti_grid_m[I] = 0.0
            self.ti_grid_v[I] = ti.Vector.zero(real, 3)
            self.ti_grid_a[I] = ti.Vector.zero(real, 3)
            self.ti_grid_x[I] = self.ti_grid_min + I * self.py_dx

        # compute grid point locations (for debug)
        for i in range(self.py_cell_count[0]*self.py_cell_count[1]):
            gi = i % self.py_cell_count[0]
            gj = (i // self.py_cell_count[0]) % self.py_cell_count[1]
            gk = i // (self.py_cell_count[0] * self.py_cell_count[1])
            I = ti.Vector([gi, gj, gk])
            self.ti_grid_pos[i] = self.ti_grid_min + I.cast(real) * self.py_dx

        # initialize particles
        for i in range(self.ti_particle_count[None]):

            pi = i % self.ti_particle_ndcount[0]
            pj = (i // self.ti_particle_ndcount[0]) % self.ti_particle_ndcount[1]
            pk = i // (self.ti_particle_ndcount[0] * self.ti_particle_ndcount[1])
            # r = ti.Vector([ti.random() - 0.5, ti.random() - 0.5, ti.random() - 0.5])
            r = ti.Vector([ 0.5, 0.5, 0.5])

            _I = ti.Vector([pi, pj, pk]).cast(real) + r
            self.ti_particle_x[i] = self.ti_particle_init_min[0] + (self.py_dx / self.py_particle_init_cell_samples_per_dim) * _I
            self.ti_particle_v[i] = self.ti_particle_init_vel
            self.ti_particle_be[i] = ti.Matrix.identity(real, 3)
            self.ti_particle_C[i] = ti.Matrix.zero(real,3, 3)
            self.ti_particle_is_inner_of_box[i] = 0

    # uGIMP basis functions
    @staticmethod
    @ti.func
    def linearIntegral(xp, hl, xi, w):
        diff = ti.abs(xp - xi)
        ret = 0.0
        if diff >= w + hl:
            ret = 0.0
        elif diff >= w - hl:
            ret = ((w + hl - diff) ** 2) / (2.0 * w)
        elif diff >= hl:
            ret = 2.0 * hl * (1.0 - diff / w)
        else:
            ret = 2.0 * hl - (hl * hl + diff * diff) / w
        return ret

    @staticmethod
    @ti.func
    def linearIntegralGrad(xp, hl, xi, w):
        diff = ti.abs(xp - xi)
        sgn = 1.0 if xp - xi >= 0.0 else -1.0
        ret = 0.0
        if diff >= w + hl:
            ret = 0.0
        elif diff >= w - hl:
            ret = -sgn * (w + hl - diff) / w
        elif diff >= hl:
            ret = -sgn * 2.0 * hl / w
        else:
            ret = 2.0 * (xi - xp) / w
        return ret

    @staticmethod
    @ti.func
    def uGIMPStencil():
        return ti.ndrange(3, 3, 3)

    @ti.func
    def uGIMPBase(self, particle_pos):
        return ((particle_pos - self.py_particle_hl - self.ti_grid_min) * self.py_invdx).cast(int)

    @ti.func
    def uGIMPWeightAndGrad(self, particle_pos, grid_pos):
        wx = self.linearIntegral(particle_pos[0], self.py_particle_hl, grid_pos[0], self.py_dx)
        wy = self.linearIntegral(particle_pos[1], self.py_particle_hl, grid_pos[1], self.py_dx)
        wz = self.linearIntegral(particle_pos[2], self.py_particle_hl, grid_pos[2], self.py_dx)
        weight = wx * wy * wz / self.py_particle_volume
        weight_grad = ti.Vector([wy * wz * self.linearIntegralGrad(particle_pos[0], self.py_particle_hl, grid_pos[0], self.py_dx), wx * wz * self.linearIntegralGrad(particle_pos[1], self.py_particle_hl, grid_pos[1], self.py_dx), wx * wy * self.linearIntegralGrad(particle_pos[2], self.py_particle_hl, grid_pos[2], self.py_dx)]) / self.py_particle_volume
        return weight, weight_grad

    @staticmethod
    @ti.func
    def bar_3d(A):
        return A / ti.pow(A.determinant(), 1.0/3.0)

    @staticmethod
    @ti.func
    def dev_3d(A):
        return A - (1.0/3.0) * A.trace() * ti.Matrix.identity(float, 3)

    @staticmethod
    @ti.func
    def hb_eval_3d(x, sigma_len_pre, mu_div_J, hb_sigma_y, hb_n, hb_eta, trace_be_bar, dt):
        return x - sigma_len_pre + ti.sqrt(2.0) * dt * mu_div_J * trace_be_bar * ti.pow( ( x / ti.sqrt(2.0) - hb_sigma_y ) / hb_eta, 1.0 / hb_n ) / 3.0

    @staticmethod
    @ti.func
    def hb_eval_deriv_3d(x, sigma_len_pre, mu_div_J, hb_sigma_y, hb_n, hb_eta, trace_be_bar, dt):
        return 1.0 + dt * mu_div_J * trace_be_bar * ti.pow( ( x / ti.sqrt(2.0) - hb_sigma_y ) / hb_eta, 1.0 / hb_n - 1.0 ) / (3.0 * hb_n * hb_eta)

    @ti.func
    def scalar_hb_solve_3d(self, sigma_len_pre, mu_div_J, hb_sigma_y, hb_n, hb_eta, trace_be_bar, dt):
        x = sigma_len_pre

        #while True:
        for i in range(14):
            fx = self.hb_eval_3d(x, sigma_len_pre, mu_div_J, hb_sigma_y, hb_n, hb_eta, trace_be_bar, dt)
            dfx = self.hb_eval_deriv_3d(x, sigma_len_pre, mu_div_J, hb_sigma_y, hb_n, hb_eta, trace_be_bar, dt)
            dx = - fx / dfx

            for j in range(20):
                x_new = x + dx
                if ( x_new / ti.sqrt(2.0) - hb_sigma_y ) >= 0:
                    x = x_new
                    break
                dx = dx / 2.0    

            if ti.abs(dx) < 1.0e-6:
                break

        return x

    @ti.func
    def isnan(self, x):
        return not (x < 0 or 0 < x or x == 0)

    @ti.kernel
    def step(self):
        self.ti_iteration[None] += 1

        # clear grid data
        for I in ti.grouped(self.ti_grid_m):
            self.ti_grid_m[I] = 0.0
            self.ti_grid_v[I] = ti.Vector.zero(real, 3)
            self.ti_grid_a[I] = ti.Vector.zero(real, 3)

        # particle status update and p2g
        for p in range(self.ti_particle_count[None]):
            base = self.uGIMPBase(self.ti_particle_x[p])
            stencil = self.uGIMPStencil()
            
            # compute particle stress
            J = ti.sqrt(self.ti_particle_be[p].determinant())
            be_bar = self.ti_particle_be[p] * pow(J, -2.0/3.0)
            dev_be_bar = be_bar - be_bar.trace() * ti.Matrix.identity(real, 3) / 3.0
            tau = self.py_kappa * 0.5 * (J+1.0) * (J-1.0) * ti.Matrix.identity(real, 3) + self.py_mu * dev_be_bar

            # p2g
            for i, j, k in ti.static(stencil):
                offset = ti.Vector([i, j, k])
                # grid point position
                gp = self.ti_grid_min + (base + offset).cast(real) * self.py_dx

                # compute weight and weight grad
                weight, weight_grad = self.uGIMPWeightAndGrad(self.ti_particle_x[p], gp)

                #internal force   
                f_internal = - self.py_particle_volume * tau @ weight_grad

                # accumulate grid velocity, acceleration and mass
                self.ti_grid_v[base + offset] += weight * self.py_particle_mass * ( self.ti_particle_v[p] + self.ti_particle_C[p] @ ( gp - self.ti_particle_x[p] ) )
                self.ti_grid_a[base + offset] += f_internal
                self.ti_grid_m[base + offset] += weight * self.py_particle_mass

        # grid update
        for I in ti.grouped(self.ti_grid_m):
            if self.ti_grid_m[I] > 0:
                old_momentum = self.ti_grid_v[I]
                new_momentum = old_momentum + self.py_dt * ( self.ti_grid_a[I] + self.ti_grid_m[I] * self.ti_g )

                # boundary conditions
                for s in range(self.ti_num_boxes[None]):
                    if self.ti_static_box_min[s][0] <= self.ti_grid_x[I][0] <= self.ti_static_box_max[s][0]:
                        if self.ti_static_box_min[s][1] <= self.ti_grid_x[I][1] <= self.ti_static_box_max[s][1]:
                            if self.ti_static_box_min[s][2] <= self.ti_grid_x[I][2] <= self.ti_static_box_max[s][2]:
                                new_momentum = ti.Vector.zero(real, 3)

                self.ti_grid_v[I] = new_momentum / self.ti_grid_m[I]
                self.ti_grid_a[I] = ( new_momentum - old_momentum ) / ( self.ti_grid_m[I] * self.py_dt )

        # g2p and update deformation status
        for p in range(self.ti_particle_count[None]):
            base = self.uGIMPBase(self.ti_particle_x[p])
            stencil = self.uGIMPStencil()

            v_pic = ti.Vector.zero(real, 3)
            grid_a = ti.Vector.zero(real, 3)
            vel_grad = ti.Matrix.zero(real, 3, 3)

            # compute velocity gradient and particle velocity
            for i, j, k in ti.static(stencil):
                offset = ti.Vector([i, j, k])
                # grid point position
                gp = self.ti_grid_min + (base + offset).cast(real) * self.py_dx

                # compute weight and weight grad
                weight, weight_grad = self.uGIMPWeightAndGrad(self.ti_particle_x[p], gp)

                vel_grad += self.ti_grid_v[base + offset].outer_product(weight_grad)
                v_pic += weight * self.ti_grid_v[base + offset]
                grid_a += weight * self.ti_grid_a[base + offset]

            self.ti_particle_v[p] = v_pic
            self.ti_particle_C[p] = vel_grad

            # elastic prediction
            f = ti.Matrix.identity(float, 3) + self.py_dt * vel_grad
            f_bar = self.bar_3d(f)
            be_bar = self.bar_3d(self.ti_particle_be[p])
            be_bar_pre = f_bar @ be_bar @ f_bar.transpose()

            be = f @ self.ti_particle_be[p] @ f.transpose()
            det_be = be.determinant()
            J = ti.sqrt(det_be)

            sigma_s_pre = self.py_mu * self.dev_3d(be_bar_pre) / J
            sigma_s_pre_len = sigma_s_pre.norm()

            scalar_sigma_pre = sigma_s_pre_len / ti.sqrt(2.0)

            # plastic correction
            if scalar_sigma_pre - self.ti_hb_sigmaY[None] > 0.0:
                sigma_s_pre_hat = sigma_s_pre / sigma_s_pre_len
                sigma_s_len = self.scalar_hb_solve_3d(sigma_s_pre_len, self.py_mu / J, self.ti_hb_sigmaY[None], self.ti_hb_n[None], self.ti_hb_eta[None], be_bar.trace(), self.py_dt)

                be_bar = (be_bar.trace() / 3.0) * ti.Matrix.identity(float, 3) + sigma_s_len * J * sigma_s_pre_hat / self.py_mu
                det_be_bar = be_bar.determinant()
                be = be_bar * ti.pow(det_be, 1.0/3.0) / ti.pow(det_be_bar, 1.0 / 3.0)


            self.ti_particle_be[p] = be
               
            # boundary conditions
            for s in range(self.ti_num_boxes[None]):
                if self.ti_static_box_min[s][0] <= self.ti_particle_x[p][0]<= self.ti_static_box_max[s][0]:
                    if self.ti_static_box_min[s][1] <= self.ti_particle_x[p][1] <= self.ti_static_box_max[s][1]:
                        if self.ti_static_box_min[s][2] <= self.ti_particle_x[p][2] <= self.ti_static_box_max[s][2]:
                            self.ti_particle_v[p] = ti.Vector.zero(real, 3)
                            self.ti_particle_C[p] = ti.Matrix.zero(real, 3, 3)
                            self.ti_particle_is_inner_of_box[p] = 1
                        else:
                            self.ti_particle_is_inner_of_box[p] = 0

            # advect
            self.ti_particle_x[p] += self.py_dt * self.ti_particle_v[p]

def system(cmd):
        proc = subprocess.run(cmd, shell=True, encoding='utf-8', stdout=PIPE, stderr=PIPE, text=True)
        print(proc.stdout)
        print(proc.stderr)

def nohup(cmd):
        proc = subprocess.Popen(cmd, shell=True, encoding='utf-8', stdout=PIPE, stderr=PIPE, text=True)


#@ti.data_oriented
class fileOperation:
    static_to_processing_dir = '/to_process/'
    static_processing_dir = '/processing/'
    static_processed_dir = '/processed/'
    static_file_wild_card = '*.xml'

    def __init__(self, root_dir_path):
        self.py_saved_iteration = 0
        self.py_filename = ''
        self.py_save_count = 1

        self.py_root_dir_path = str(root_dir_path)
        self.py_file_processing = ''

    def initialize(self, filename):
        self.py_saved_iteration = 0
        self.py_filename = filename
        self.py_save_count = 1

    def fetchFileToProcess(self):
        files = glob.glob(self.py_root_dir_path + "/**" + self.static_to_processing_dir + self.static_file_wild_card, recursive=True)
        if len(files) > 0:
            self.initialize(files[0])
        else:
            self.initialize('')

    def getFileName(self):
        return self.py_filename

    def moveFileToProcessing(self):
        path = Path(self.py_filename)
        path_processing = str(path.parent.parent) + self.static_processing_dir
        self.py_file_processing = path_processing + ntpath.basename(self.py_filename)
        shutil.move(self.py_filename, path_processing)

    def moveFileToProcessed(self):
        path = Path(self.py_filename)
        path_processed = str(path.parent.parent) + self.static_processed_dir
        shutil.move(self.py_file_processing, path_processed)

    def copy(self,src,dest):
        shutil.copy(src, dest)

    def saveState(self , agTaichiMPM):
        saveStateFilePath = os.path.dirname(self.py_filename) + '/config_' + str(self.py_save_count).zfill(2) + ".dat"
        saveStateIntermediateFilePath = os.path.dirname(self.py_filename) + '/config_' + str(self.py_save_count).zfill(2) + "_phi" + ".dat"
        outObjFilePath = os.path.dirname(self.py_filename) + '/config_' + str(self.py_save_count).zfill(2) + ".obj"

        if os.path.exists(saveStateFilePath):
            print("remove: " + saveStateFilePath)
            os.remove(saveStateFilePath)
        if os.path.exists(saveStateIntermediateFilePath):
            print("remove: " + saveStateIntermediateFilePath)
            os.remove(saveStateIntermediateFilePath)
        if os.path.exists(outObjFilePath):
            print("remove: " + outObjFilePath)
            os.remove(outObjFilePath)
        marching_cube_path = os.path.dirname(agTaichiMPM.py_particleSkinnerApp) + "/cpp_marching_cubes"
       
        print('[AGTaichiMPM] saving state to ' + saveStateFilePath)
        f = open(saveStateFilePath, 'wb')
        particle_is_inner_of_box_id = np.where(agTaichiMPM.ti_particle_is_inner_of_box.to_numpy()[0:agTaichiMPM.ti_particle_count[None]].astype(np.int32) == 1)
        f.write(ctypes.c_int32(agTaichiMPM.ti_particle_count[None] -  particle_is_inner_of_box_id[0].size))
        #output x
        p_x = agTaichiMPM.ti_particle_x.to_numpy()[0:agTaichiMPM.ti_particle_count[None]].astype(np.float32)
        np.delete(p_x, particle_is_inner_of_box_id,axis=0).flatten().tofile(f)
        #output radius
        np.delete((np.ones(agTaichiMPM.ti_particle_count[None], np.float32) * agTaichiMPM.py_particle_hl).astype(np.float32), particle_is_inner_of_box_id,axis=0).flatten().tofile(f)
        #output velocity
        np.delete(agTaichiMPM.ti_particle_v.to_numpy()[0:agTaichiMPM.ti_particle_count[None]].astype(np.float32), particle_is_inner_of_box_id,axis=0).flatten().tofile(f)
        #output id
        np.delete(np.ones(agTaichiMPM.ti_particle_count[None], ctypes.c_int32), particle_is_inner_of_box_id,axis=0).flatten().tofile(f)
        f.close()
        #exec particle skinner
        
        cmd = 'python3 "' + agTaichiMPM.py_particleSkinnerApp + '" ' + str(agTaichiMPM.py_skin_size) + ' "' + saveStateFilePath + '" "' + saveStateIntermediateFilePath + '" "' + outObjFilePath + '" "' + marching_cube_path + '"'

        if agTaichiMPM.py_emulation:
            cmd += '; "' + agTaichiMPM.py_GLRenderApp + '" -a "' + agTaichiMPM.py_emu_dir_path \
            + '" -b "' + outObjFilePath + '" -c ' + str(self.py_save_count).zfill(2) \
            + ' -d ' + agTaichiMPM.py_camera_position + ' -e ' + agTaichiMPM.py_camera_quat \
            + ' -f ' + agTaichiMPM.py_camera_fov + ' -g ' + agTaichiMPM.py_camera_window
            print(cmd)
            nohup(cmd)
        else:
            refPngFilePath = agTaichiMPM.py_ref_dir_path + '/config_' + str(self.py_save_count).zfill(2) + ".png"
            outdiff_dir = agTaichiMPM.py_sim_dir_path
            cmd += '; "' + agTaichiMPM.py_GLRenderApp + '" -a "' + outdiff_dir\
            + '" -b "' + refPngFilePath + '" -c "' + outObjFilePath\
            + '" -d ' + str(self.py_save_count).zfill(2)\
            + ' -e ' + agTaichiMPM.py_camera_position + ' -f ' + agTaichiMPM.py_camera_quat\
            + ' -g ' + agTaichiMPM.py_camera_fov + ' -h ' + agTaichiMPM.py_camera_window
            print(cmd)
            nohup(cmd)
        self.py_save_count += 1
        
        isNan = np.isnan(p_x.flatten()).any(axis=0)
        if isNan == True:
            sys.exit("Raised Nan Error.")


def T(a):
    phi, theta = np.radians(28), np.radians(32)

    a = a - 0.5
    x, y, z = a[:, 0], a[:, 1], a[:, 2]
    c, s = np.cos(phi), np.sin(phi)
    C, S = np.cos(theta), np.sin(theta)
    x, z = x * c + z * s, z * c - x * s
    u, v = x, y * C + z * S
    return np.array([u, v]).swapaxes(0, 1) + 0.5

gui = ti.GUI("AGTaichiMPM")

if len(sys.argv) <= 2:
    print("usage: python AGTaichiMPM.py <path to simulation> <initial setting xml>")
    exit(-1)

fileOp = fileOperation(sys.argv[1])
init_xmlData = initXmlParser.MPMXMLData(sys.argv[2])
agTaichiMPM = AGTaichiMPM(init_xmlData)
cc = 0

while True:
    if mode == Mode.IDLE:
        print('idle... [' + str(count) + ']')
        fileOp.fetchFileToProcess()
        fn = fileOp.getFileName()
        files = glob.glob(fileOp.py_root_dir_path + "/**/*.txt", recursive=True)
        for f in files:
            if os.path.basename(f) == "exit.txt":
                sys.exit()
        if  fn != '':
            print(fileOp.py_root_dir_path)

            for _ in range(5):
                try:
                    xmlData = xmlParser.MPMXMLData(fn)
                except Exception as e:
                    print(e)
                    tm.sleep(1.0)
                else:
                    break
            else:
                print("XML parse error occur.")
                exit()

            xmlData.show()
            agTaichiMPM.changeXMLData(xmlData)
            agTaichiMPM.initialize()
            # cleaning up idle mode
            count = 0
            # preparing for switching to processing mode
            print('found :' + fn)
            fileOp.moveFileToProcessing()
            mode = Mode.PROCESSING
            print('processing phase will start')
            continue
        count += 1
        tm.sleep(0.5)

    elif mode == Mode.PROCESSING:
        time = agTaichiMPM.ti_iteration[None] * agTaichiMPM.py_dt
        print("stepping : " + str(time))

        if agTaichiMPM.py_num_saved_frames > agTaichiMPM.py_max_frames:
            agTaichiMPM.py_num_saved_frames = 1
            fileOp.moveFileToProcessed()
            gui.circles(T(agTaichiMPM.ti_particle_x.to_numpy() / 20 + 0.3), radius=2, color=0xFFFFFF)
            path = Path(fileOp.py_filename)
            path_processed = str(path.parent.parent) + fileOp.static_processed_dir
            gui.show(path_processed + "/ti_endframe" + str(cc) + ".png")
            #if not emulation, output lossvalues.dat
            cc+=1
            mode = Mode.IDLE
            print('process done. turning to idle mode')
        else:
            for i in range(100):
                agTaichiMPM.step()
                time = agTaichiMPM.ti_iteration[None] * agTaichiMPM.py_dt

                if time * agTaichiMPM.py_fps >= agTaichiMPM.py_num_saved_frames:
                    fileOp.py_saved_iteration = agTaichiMPM.ti_iteration[None]
                    fileOp.saveState(agTaichiMPM)
                    agTaichiMPM.py_num_saved_frames += 1
            gui.circles(T(agTaichiMPM.ti_particle_x.to_numpy() / 20 + 0.3), radius=2, color=0xFFFFFF)
            gui.show()
