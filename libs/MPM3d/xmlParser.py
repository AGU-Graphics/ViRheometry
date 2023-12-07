import numpy as np
import xml.etree.ElementTree as ET


class MPMXMLStaticBoxData:
    def __init__(self, dim):
        self.min = np.zeros(dim)
        self.max = np.zeros(dim)
        self.isSticky = False

    def setFromXMLNode(self, node):
        self.min = np.array([float(e) for e in node.attrib['min'].split(' ')])
        self.max = np.array([float(e) for e in node.attrib['max'].split(' ')])
        if node.attrib['boundary_behavior'] == 'sticking':
            self.isSticky = True
        else:
            self.isSticky = False

    def show(self):
        print('*** Static box ***')
        print('  min: ' + str(self.min))
        print('  max: ' + str(self.max))
        print('  isSticky: ' + str(self.isSticky))


class MPMXMLEmulationData:
    def __init__(self):
        self.emulation = 0
        self.emu_dir_path = ""

    def setFromXMLNode(self, node):
        self.emulation = float(node.attrib['emulation'])
        self.emu_dir_path = str(node.attrib['emu_dir_path'])


    def show(self):
        print('*** emulation ***')
        print('  emulation      : ' + str(self.emulation))
        print('  emu_dir_path      : ' + str(self.emu_dir_path))


class MPMXMLPathData:
    def __init__(self):
        self.sim_dir_path = ""
        self.ref_dir_path = ""

    def setFromXMLNode(self, node):
        self.ref_dir_path = str(node.attrib['ref_dir_path'])
        self.sim_dir_path = str(node.attrib['sim_dir_path'])

    def show(self):
        print('*** Path ***')
        print('  ref_dir_path      : ' + str(self.ref_dir_path))
        print('  sim_dir_path      : ' + str(self.sim_dir_path))

class MPMXMLparticleSkinnerData:
    def __init__(self):
        self.path = ""
        self.grid_space = 0.0
        self.file_type = ""

    def setFromXMLNode(self, node):
        self.path = str(node.attrib['path'])
        self.grid_space = float(node.attrib['grid_space'])
        self.file_type = str(node.attrib['file_type'])


    def show(self):
        print('*** Particle Skinner Arguments ***')
        print('  path      : ' + str(self.path))
        print('  grid_space: ' + str(self.grid_space))
        print('  file_type : ' + str(self.file_type))


class MPMXMLGLRenderData:

    class MPMXMLCameraData:
        def __init__(self):
            self.eyepos = np.zeros(3)
            self.quat = np.zeros(4)
            self.window_size = np.zeros(2)
            self.fov = 0.0

        def setFromXMLNode(self, node):
            self.eyepos = np.array([float(e) for e in node.attrib['eyepos'].split(' ')])
            self.quat = np.array([float(e) for e in node.attrib['quat'].split(' ')])
            self.window_size = np.array([float(e) for e in node.attrib['window_size'].split(' ')])
            self.fov = float(node.attrib['fov'])

        def show(self):
            print('  camera :')
            print('    eyepos      : ' + str(self.eyepos))
            print('    quat: ' + str(self.quat))
            print('    window_size : ' + str(self.window_size))
            print('    fov : ' + str(self.fov))

    def __init__(self):
        self.path = ""
        self.cameraData = self.MPMXMLCameraData()

    def setFromXMLNode(self, node):
        self.path = str(node.attrib['path'])
        self.cameraData.setFromXMLNode(node.find('camera'))

    def show(self):
        print('*** GLRender Arguments ***')
        print('  path      : ' + str(self.path))
        self.cameraData.show()

class MPMXMLIntegratorData:
    def __init__(self):
        self.herschel_bulkley_power = 1.0 
        self.eta = 0.0 
        self.yield_stress = 0.0 


    def setFromXMLNode(self, node):
        self.herschel_bulkley_power = float(node.attrib['herschel_bulkley_power'])
        self.eta = float(node.attrib['eta'])
        self.yield_stress = float(node.attrib['yield_stress'])

    def show(self):
        print('*** Integrator ***')
        print('  herschel_bulkley_power: ' + str(self.herschel_bulkley_power))
        print('  eta: ' + str(self.eta))
        print('  yield_stress: ' + str(self.yield_stress))

class MPMXMLCuboidData:
    def __init__(self, dim):
        self.min = np.zeros(dim)
        self.max = np.zeros(dim)
        self.density = 1.0
        self.cell_samples_per_dim = 2
        self.vel = np.zeros(dim)

    def setFromXMLNode(self, node):
        self.min = np.array([float(e) for e in node.attrib['min'].split(' ')])
        self.max = np.array([float(e) for e in node.attrib['max'].split(' ')])
        self.density = float(node.attrib['density'])
        self.cell_samples_per_dim = int(node.attrib['cell_samples_per_dim'])
        self.vel = np.array([float(e) for e in node.attrib['vel'].split(' ')])

    def show(self):
        print('*** Cuboid ***')
        print('  min: ' + str(self.min))
        print('  max: ' + str(self.max))
        print('  density: ' + str(self.density))
        print('  cell_samples_per_dim: ' + str(self.cell_samples_per_dim))
        print('  vel: ' + str(self.vel))

class MPMXMLStaticPlaneData:
    def __init__(self, dim):
        self.x = np.zeros(dim)
        self.n = np.zeros(dim)
        self.isSticky = False

    def setFromXMLNode(self, node):
        self.x = np.array([float(e) for e in node.attrib['x'].split(' ')])
        self.n = np.array([float(e) for e in node.attrib['n'].split(' ')])
        if node.attrib['boundary_behavior'] == 'sticking':
            self.isSticky = True
        else:
            self.isSticky = False

    def show(self):
        print('*** Static plane ***')
        print('  x: ' + str(self.x))
        print('  n: ' + str(self.n))
        print('  isSticky: ' + str(self.isSticky))

class MPMXMLNearEarthGravityData:
    def __init__(self, dim):
        self.g = np.zeros(dim)

    def setFromXMLNode(self, node):
        self.g = np.array([float(e) for e in node.attrib['f'].split(' ')])

    def show(self):
        print('*** Near earth gravity ***')
        print('  g: ' + str(self.g))

class MPMXMLData:
    def __init__(self, file_name):
        print('[AGTaichiMPM3D] Parsing xml file: ' + str(file_name))
        tree = ET.parse(file_name)
        root = tree.getroot()
        if root.tag != 'AGTaichiMPM3D':
            print('[AGTaichiMPM3D] Could not find root note AGTaichiMPM3D. Exiting...')
            exit(-1)

        self.staticBoxList = []
        for static_box in root.findall('static_box'):
            staticBoxData = MPMXMLStaticBoxData(3)
            staticBoxData.setFromXMLNode(static_box)
            self.staticBoxList.append(staticBoxData)

        emulation = root.find('emulation')
        self.emulationData = MPMXMLEmulationData()
        self.emulationData.setFromXMLNode(emulation)
            
        path = root.find('path')
        self.pathData = MPMXMLPathData()
        self.pathData.setFromXMLNode(path)

        ps = root.find('particle_skinner')
        self.particleSkinnerData = MPMXMLparticleSkinnerData()
        self.particleSkinnerData.setFromXMLNode(ps)

        GLRender = root.find('GLRender')
        self.GLRenderData = MPMXMLGLRenderData()
        self.GLRenderData.setFromXMLNode(GLRender)
        
        integrator = root.find('integrator')
        self.integratorData = MPMXMLIntegratorData()
        self.integratorData.setFromXMLNode(integrator)

        cuboid = root.find('cuboid')
        self.cuboidData = MPMXMLCuboidData(3)
        self.cuboidData.setFromXMLNode(cuboid)


    def show(self):
        print('[AGTaichiMPM3D] XML Data:')
        self.emulationData.show()
        self.pathData.show()
        self.particleSkinnerData.show()
        self.GLRenderData.show()
        self.integratorData.show()
        self.cuboidData.show()
        for sb in self.staticBoxList:
            sb.show()

