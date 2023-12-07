import bpy
import math
import bpy
import os
from bpy.props import StringProperty
from bpy_extras.io_utils import ImportHelper
from bpy.types import Operator

EXPORT_CONTENTS = ""

def Cam(offset, cam_arr):
    for objname in cam_arr:
        bpy.data.objects[objname].select_set(True) # 2.8+
    print("<!-- export camera -->")
    for obj in bpy.context.selected_objects:
        EP=[]
        QA=[]
        WS=[]
        FV=[]
        
        EP.append(obj.location[0])
        EP.append(obj.location[2])
        EP.append(-obj.location[1])

        theta = 2.0 * math.acos(obj.rotation_quaternion[0]) - math.pi * 0.5
        _axis_x = obj.rotation_quaternion[1]
        _axis_y = obj.rotation_quaternion[2]
        _axis_z = obj.rotation_quaternion[3]
        axis_length = math.sqrt(_axis_x*_axis_x+_axis_y*_axis_y+_axis_z*_axis_z)
        
        axis_x = _axis_x / axis_length
        axis_y = _axis_z / axis_length
        axis_z = -_axis_y / axis_length
        
        q0 = math.cos(theta*0.5)
        q1 = math.sin(theta*0.5) * axis_x
        q2 = math.sin(theta*0.5) * axis_y
        q3 = math.sin(theta*0.5) * axis_z

        QA.append(q0)
        QA.append(q1)
        QA.append(q2)
        QA.append(q3)
        
        fovy = math.atan(20.25 / (2.0 * obj.data.lens)) * 2.0 * 180.0 / math.pi

        WS.append(960)
        WS.append(540)
        
        global EXPORT_CONTENTS
        EXPORT_CONTENTS += '<camera eyepos="{:.8f} {:.8f} {:.8f}"  quat="{:.8f} {:.8f} {:.8f} {:.8f}" window_size="{:.0f} {:.0f}" fov="{:.5f}"/>\n'.format(EP[0] + offset[0],EP[1] + offset[1],EP[2] + offset[2],QA[0],QA[1],QA[2],QA[3],WS[0],WS[1],fovy)
        
def Obj(overlap_thickness, side_wall_thickness, obj_arr):
    for objname in obj_arr:
        bpy.data.objects[objname].select_set(True) # 2.8+
    
    for obj in bpy.context.selected_objects:
        X=[]
        Y=[]
        Z=[]
        Scale = obj.scale

        obj.data
        for i in obj.data.vertices:
            X.append(i.co.x)
            Y.append(i.co.y)
            Z.append(i.co.z)
        X.sort()
        Y.sort()
        Z.sort()

        # transform to opengl coordinates
        coord_min = [X[0]*Scale[0], Z[0]*Scale[2], -Y[7]*Scale[1]]
        coord_max = [X[7]*Scale[0], Z[7]*Scale[2], -Y[0]*Scale[1]]
        
        offset_x = coord_max[0] - coord_min[0]
        offset_y = 0.0
        offset_z = coord_max[2] - coord_min[2]

        return [offset_x, offset_y, offset_z]


class FileWindow(Operator, ImportHelper):

    bl_idname = "file.file_window"
    bl_label = "Open the file window"
    
    directory: StringProperty(
        name="Directory Path", 
        default="",            
        maxlen=1024,           
        subtype='FILE_PATH',  
        description="",        
    )

    def execute(self, context):
        dir = self.filepath
        global EXPORT_CONTENTS
        with open(dir + "/camera_params.xml", "w") as f:
            f.write(EXPORT_CONTENTS)
            f.close()
        return {'FINISHED'}


def register():
    bpy.utils.register_class(FileWindow)


def unregister():
    bpy.utils.unregister_class(FileWindow)


if __name__ == "__main__":
    register()
    bpy.ops.file.file_window('INVOKE_DEFAULT')
    OBJNAME=[]
    ObjName=[]
    CamName=[]

    for obj in bpy.context.selected_objects:
        OBJNAME.append(obj.name_full)
        #print(obj.data.name_full)
        
    for objname in OBJNAME:
        if objname.startswith("Obj"):
            ObjName.append(objname)
            
    for objname in OBJNAME:
        if objname.startswith("Cam"):
            CamName.append(objname)
                
    overlap_thickness = 0.15
    side_wall_thickness = 0.3
    coord_offset = [0.0, 0.0, 0.0]
    EXPORT_CONTENTS += "<?xml version=\"1.0\"?>\n"
    EXPORT_CONTENTS += "<setup3D>\n"

    bpy.ops.object.select_all(action='DESELECT')
    offset = Obj(overlap_thickness, side_wall_thickness, ObjName)

    bpy.ops.object.select_all(action='DESELECT')
    Cam(offset, CamName)
    EXPORT_CONTENTS += "</setup3D>\n"