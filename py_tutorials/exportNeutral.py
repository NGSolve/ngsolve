import sys

def Export (mesh, filename):
    """ export Netgen mesh to neutral format """
    
    print ("export mesh in neutral format to file = ", filename)

    f = open (filename, 'w')

    points = mesh.Points()
    print (len(points), file=f)
    for p in points:
        print (p.p[0], p.p[1], p.p[2], file=f)


    volels = mesh.Elements3D();
    print (len(volels), file=f)
    for el in volels:
        print (el.index, end="   ", file=f)
        for j in el.vertices:
            print (j.nr, end=" ", file=f)
        print(file=f)




