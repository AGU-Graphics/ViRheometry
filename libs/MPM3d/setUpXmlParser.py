import numpy as np
import xml.etree.ElementTree as ET


class SetUpXMLCameraData:
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


class SetUpXMLData:
    def __init__(self, file_name):
        print('[3D Camera] Parsing xml file: ' + str(file_name))
        tree = ET.parse(file_name)
        root = tree.getroot()
        if root.tag != 'setup3D':
            print('[AGTaichi SetUp3D] Could not find root note AGTaichi SetUp3D. Exiting...')
            exit(-1)
        
        camera = root.find('camera')
        self.cameraData =  SetUpXMLCameraData()
        self.cameraData.setFromXMLNode(camera)


    def show(self):
        print('[AGTaichi SetUp3D] XML Data:')
        self.cameraData.show()
