
import subprocess
import os



def convert_vtk(reader,file):
    mesh = reader.GetOutput()
    points = mesh.GetPoints()
    n_points = points.GetNumberOfPoints()
    bounds=mesh.GetBounds()
    dx = (bounds[1] - bounds[0])
    dy = (bounds[3] - bounds[2])
    dz = (bounds[5] - bounds[4])

    VTK_TRIANGLE = 5
    triangles = []
    for cell_id in range(mesh.GetNumberOfCells()):
        cell = mesh.GetCell(cell_id)
        if cell.GetCellType() == VTK_TRIANGLE:
            point_ids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
            triangles.append((cell_id, *point_ids))

    with open(file+".tsi","w",encoding="UTF8") as f:
        f.write('version 2.0\n')
        f.write('box\t{}\t{}\t{}\n'.format(dx*2,dy*2,dz*2))
        f.write('vertex\t{}\n'.format(n_points))
        for i in range(n_points):
            x,y,z=points.GetPoint(i)
            f.write('{}\t{:.6f}\t{:.6f}\t{:.6f}\n'.format(i, x+dx, y+dy, z+dz))
        

        f.write('triangle\t{}\n'.format(len(triangles)))
        for t in triangles:
            f.write('{}\t'.format(t[0]))
            for k in [1,2,3]:
                f.write('{}\t'.format(t[k]))
            f.write('\n')

        f.write('inclusion\t0\n')


def adjust_mesh_file(args):
    file = args[args.index("-TSfile")+1]
    _,file_extension = os.path.splitext(file)
    if file_extension == ".tsi":
        return args
    elif file_extension == ".vtk" or file_extension == ".vtu":
        import vtk
        if file_extension == ".vtk":
            reader = vtk.vtkUnstructuredGridReader()
        else:
            reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(file)
        reader.Update()
        convert_vtk(reader,file)
        args[args.index("-TSfile")+1]=file+".tsi"

    return args