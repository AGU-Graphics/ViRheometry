import xml.etree.ElementTree as ET

class OptXMLPathData:
    def __init__(self):
        self.root_dir_path = ""
        self.GL_emulation_render_path = ""
        self.GL_render_path = ""
        self.particle_skinner_path = ""
        self.mpm_path = ""
        self.shell_script_dir_path = ""
        
    def setFromXMLNode(self, node):
        self.root_dir_path = str(node.attrib['root_dir_path'])
        self.GL_emulation_render_path = str(node.attrib['GL_emulation_render_path'])
        self.GL_render_path = str(node.attrib['GL_render_path'])
        self.particle_skinner_path = str(node.attrib['particle_skinner_path'])
        self.mpm_path = str(node.attrib['mpm_path'])
        self.shell_script_dir_path = str(node.attrib['shell_script_dir_path'])

    def show(self):
        print('*** Path ***')
        print('  root_dir_path      : ' + str(self.root_dir_path))
        print('  GL_emulation_render_path      : ' + str(self.GL_emulation_render_path))
        print('  GL_render_path      : ' + str(self.GL_render_path))
        print('  particle_skinner_path      : ' + str(self.particle_skinner_path))
        print('  mpm_path      : ' + str(self.mpm_path))
        print('  shell_script_dir_path      : ' + str(self.shell_script_dir_path))


class OptXMLInitialMaterialData:
    def __init__(self):
        self.eta = 0.0 
        self.n = 0.0 
        self.sigmaY = 0.0 

    def setFromXMLNode(self, node):
        self.eta = float(node.attrib['eta'])
        self.n = float(node.attrib['n'])
        self.sigmaY = float(node.attrib['sigmaY'])

    def show(self):
        print('*** initial material ***')
        print('  eta: ' + str(self.eta))
        print('  n: ' + str(self.n))
        print('  sigmaY: ' + str(self.sigmaY))

class OptXMLSetupData:
    def __init__(self):
        self.H = 0.0
        self.W = 0.0
        self.rho = 1.0

    def setFromXMLNode(self, node):
        self.H = float(node.attrib['H'])
        self.W = float(node.attrib['W'])
        self.rho = float(node.attrib['RHO'])

    def show(self):
        print('*** setup ***')
        print('  H: ' + str(self.H))
        print('  W: ' + str(self.W))
        print('  rho: ' + str(self.rho))

class OptXMLData:
    def __init__(self, file_name):
        print('[Optimizer] Parsing xml file: ' + str(file_name))
        tree = ET.parse(file_name)
        root = tree.getroot()
        if root.tag != 'Optimizer':
            print('[Optimizer] Could not find root note Optimizer. Exiting...')
            exit(-1)
        
            
        path = root.find('path')
        self.pathData = OptXMLPathData()
        self.pathData.setFromXMLNode(path)

        init_setup = root.find('setup')
        self.setupData = OptXMLSetupData()
        self.setupData.setFromXMLNode(init_setup)

        initial_material = root.find('initial_material')
        if not initial_material == None:
            self.initialMaterialData = OptXMLInitialMaterialData()
            self.initialMaterialData.setFromXMLNode(initial_material)

    def show(self):
        print('[Optimizier] XML Data:')
        self.pathData.show()
        self.setupData.show()
        self.initialMaterialData.show()

        for s in self.setupList:
            s.show()


  
