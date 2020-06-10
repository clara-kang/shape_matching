import bpy

scene = bpy.context.scene

for obj in scene.objects:
    if obj.name.startswith("v_"):
        obj.select = True
