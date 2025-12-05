import bpy

ob = bpy.context.active_object
bpy.ops.object.origin_set(type="GEOMETRY_ORIGIN")

box = []
for i in ob.dimensions:
    box.append(i*2)

f = open('F9bK3qZ2L7xJ1D0wH8sV4P5aR.tsi', 'w')

f.write('version 1.1\n')
f.write('box\t{}\t{}\t{}\n'.format(box[0], box[1], box[2]))
f.write('vertex\t{}\n'.format(len(ob.data.vertices)))

for v in ob.data.vertices:
    f.write('{}\t{:.6f}\t{:.6f}\t{:.6f}\n'.format(v.index, v.co.x+box[0]/2, v.co.y+box[1]/2, v.co.z+box[2]/2))

f.write('triangle\t{}\n'.format(len(ob.data.polygons)))
for t in ob.data.polygons:
    f.write('{}\t'.format(t.index))
    for v in t.vertices:
        f.write('{}\t'.format(v))
    f.write('\n')

f.write('inclusion\t0\n')
