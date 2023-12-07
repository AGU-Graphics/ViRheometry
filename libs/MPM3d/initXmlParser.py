import numpy as np
import xml.etree.ElementTree as ET

class MPMXMLIntegratorData:
    def __init__(self):
        self.dt = 0.0
        self.bulk_modulus = 1.0 
        self.shear_modulus = 1.0 
        self.flip_pic_alpha = 0.95

        self.fps = 50
        self.max_frames = 4
        self.optimization_emulation = 0

    def setFromXMLNode(self, node):
        self.dt = float(node.attrib['dt'])
        self.bulk_modulus = float(node.attrib['bulk_modulus'])
        self.shear_modulus = float(node.attrib['shear_modulus'])
        self.flip_pic_alpha = float(node.attrib['flip_pic_alpha'])

        self.fps = int(node.attrib['fps'])
        self.max_frames = int(node.attrib['max_frames'])
        self.optimization_emulation = int(node.attrib['optimization_emulation'])

    def show(self):
        print('*** Integrator ***')
        print('  dt: ' + str(self.dt))
        print('  bulk_modulus: ' + str(self.bulk_modulus))
        print('  shear_modulus: ' + str(self.shear_modulus))
        print('  flip_pic_alpha: ' + str(self.flip_pic_alpha))

        print('  fps: ' + str(self.fps))
        print('  max_frames: ' + str(self.max_frames))
        print('  optimization_emulation: ' + str(self.optimization_emulation))


class MPMXMLGridData:
    def __init__(self, dim):
        self.min = np.zeros(dim)
        self.max = np.zeros(dim)
        self.cell_width = 1.0

    def setFromXMLNode(self, node):
        self.min = np.array([float(e) for e in node.attrib['min'].split(' ')])
        self.max = np.array([float(e) for e in node.attrib['max'].split(' ')])
        self.cell_width = float(node.attrib['cell_width'])

    def show(self):
        print('*** Grid ***')
        print('  min: ' + str(self.min))
        print('  max: ' + str(self.max))
        print('  cell_width: ' + str(self.cell_width))


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
        
        integrator = root.find('integrator')
        self.integratorData = MPMXMLIntegratorData()
        self.integratorData.setFromXMLNode(integrator)

        grid = root.find('grid')
        self.gridData = MPMXMLGridData(3)
        self.gridData.setFromXMLNode(grid)

        cuboid = root.find('cuboid')
        self.cuboidData = MPMXMLCuboidData(3)
        self.cuboidData.setFromXMLNode(cuboid)

        

        nearEarthGravity = root.find('near_earth_gravity')
        self.nearEarthGravityData = MPMXMLNearEarthGravityData(3)
        self.nearEarthGravityData.setFromXMLNode(nearEarthGravity)

    def show(self):
        print('[AGTaichiMPM3D] XML Data:')
        self.integratorData.show()
        self.gridData.show()
        self.cuboidData.show()
        self.nearEarthGravityData.show()

