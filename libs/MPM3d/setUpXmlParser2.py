import numpy as np
import xml.etree.ElementTree as ET


class SetUpXMLStaticBoxData:
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

class SetUpXMLCuboidData:
    def __init__(self, dim):
        self.min = np.zeros(dim)
        self.max = np.zeros(dim)
        self.density = 1.0
        self.cell_samples_per_dim = 2
        self.vel = np.zeros(dim)

    def setFromXMLNode(self, node):
        self.min = np.array([float(e) for e in node.attrib['min'].split(' ')])
        self.max = np.array([float(e) for e in node.attrib['max'].split(' ')])

    def show(self):
        print('*** Cuboid ***')
        print('  min: ' + str(self.min))
        print('  max: ' + str(self.max))

class SetUpXMLData:
    def __init__(self, file_name):
        print('[Cuboid, Static box] Parsing xml file: ' + str(file_name))
        tree = ET.parse(file_name)
        root = tree.getroot()
        if root.tag != 'Optimizer':
            print('[Optimizer] Could not find root note AGTaichi SetUp3D. Exiting...')
            exit(-1)
        
        self.staticBoxList = []
        for static_box in root.findall('static_box'):
            staticBoxData =  SetUpXMLStaticBoxData(3)
            staticBoxData.setFromXMLNode(static_box)
            self.staticBoxList.append(staticBoxData)

        cuboid = root.find('cuboid')
        self.cuboidData =  SetUpXMLCuboidData(3)
        self.cuboidData.setFromXMLNode(cuboid)

    def show(self):
        print('[AGTaichi Optimizer] XML Data:')
        self.cuboidData.show()
        for sb in self.staticBoxList:
            sb.show()

