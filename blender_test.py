import bpy
import bmesh

bpy.data.objects['Sphere'].select = True
bpy.context.scene.objects.active = bpy.data.objects['Sphere']

bpy.ops.object.mode_set(mode = 'EDIT')
bpy.ops.mesh.select_all(action = 'SELECT')
bpy.ops.uv.unwrap()

# bpy.ops.object.mode_set(mode='EDIT')

bpy.ops.wm.save_as_mainfile(filepath=bpy.data.filepath)
