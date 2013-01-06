#!/usr/bin/env python
"""
gmsh2bms script
===============

Usage: python gmsh2bms.py filename

Created on Wed Aug 15 16:16:41 2012
by fwagner@gfz-potsdam.de
"""

import numpy as np
import pygimli as g

# script-version only
import sys
import time

def readGmsh(fname, verbose=False):
    """
        Read Gmsh ASCII file and return instance of GIMLI::Mesh class.

        Parameters
        ----------
        fname : string
            Filename of the file to read (*.msh). The file must conform
            to the `MSH ASCII file version 2 <http://geuz.org/gmsh/doc/
            texinfo/gmsh.html#MSH-ASCII-file-format>`_ format.
        verbose : boolean
            Be verbose during import.

        Notes
        -----
        Physical groups specified in Gmsh are interpreted as follows:

        - Points with the physical number 99 are interpreted as sensors.
        - Physical Lines and Surfaces define boundaries in 2D and 3D, respectively.
            - Physical Number 1: homogenous Neumann condition
            - Physical Number 2: mixed boundary condition
            - Physical Number 3: homogeneous Dirichlet condition
            - Physical Number 4: Dirichlet condition
        - Physical Surfaces and Volumes define regions in 2D and 3D, respectively.
            - Physical Number 1: No inversion region
            - Physical Number >= 2: Inversion region
    """

    inNodes, inElements, ncount, ecount = 0, 0, 0, 0
    fid = open(fname)
    if verbose: print 'Reading %s... \n' % fname

    for line in fid:

        if line[0] == '$':
            if line.find('Nodes') > 0: inNodes = 1
            if line.find('EndNodes') > 0: inNodes = 0
            if line.find('Elements') > 0: inElements = 1
            if line.find('EndElements') > 0: inElements = 0

        else:
            if inNodes == 1:
                if len(line.split()) == 1:
                    nodes = np.zeros((int(line), 3))
                    if verbose: print '  Nodes: %s' % int(line)
                else:
                    nodes[ncount, :] = np.array(line.split(), 'float')[1:]
                    ncount += 1

            elif inElements == 1:
                if len(line.split()) == 1:
                    if verbose: print '  Entries: %s' % int(line)
                    points, lines, triangles, tets = [], [], [], []

                else:
                    entry = map(int, line.split())[1:]

                    if entry[0] == 15:
                        points.append((entry[-2], entry[-3]))
                    elif entry[0] == 1:
                        lines.append((entry[-2], entry[-1], entry[2]))
                    elif entry[0] == 2:
                        triangles.append((entry[-3], entry[-2],
                                          entry[-1], entry[2]))
                    elif entry[0] == 4:
                        tets.append((entry[-4], entry[-3], entry[-2],
                                     entry[-1], entry[2]))
    fid.close()

    lines = np.asarray(lines)
    triangles = np.asarray(triangles)
    tets = np.asarray(tets)

    if verbose:
        print '    Points: %s' % len(points)
        print '    Lines: %s' % len(lines)
        print '    Triangles: %s' % len(triangles)
        print '    Tetrahedra: %s \n' % len(tets)
        print 'Creating mesh object... \n'

    # check dimension
    if len(tets) == 0:
        dim, bounds, cells = 2, lines, triangles
        zero_dim = np.abs(nodes.sum(0)).argmin()  # identify zero dimension
    else:
        dim, bounds, cells = 3, triangles, tets
    if verbose: print '  Dimension: %s-D' % dim

    # creating instance of GIMLI::Mesh class
    mesh = g.Mesh(dim)

    # replacing boundary markers (gmsh does not allow negative physical regions)
    bound_marker = (g.MARKER_BOUND_HOMOGEN_NEUMANN, g.MARKER_BOUND_MIXED,
                    g.MARKER_BOUND_HOMOGEN_DIRICHLET, g.MARKER_BOUND_DIRICHLET)
    for i in range(4):
        bounds[:,dim][bounds[:,dim] == i+1] = bound_marker[i]

    if verbose:
        bound_types = np.unique(bounds[:,dim])
        regions = np.unique(cells[:,dim+1])
        print '  Regions: %s ' % len(regions) + str(tuple(regions))
        print '  Boundary types: %s ' % len(bound_types) + str(tuple(bound_types))

    for node in nodes:
        if dim == 2: mesh.createNode(node[0], node[3-zero_dim], 0)
        else: mesh.createNode(node)

    for cell in cells:
        if dim == 2:
            mesh.createTriangle(mesh.node(int(cell[0]-1)), mesh.node(int(cell[1]-1)),
                                mesh.node(int(cell[2]-1)), marker=int(cell[3]))
        else:
            mesh.createTetrahedron(mesh.node(int(cell[0]-1)), mesh.node(int(cell[1]-1)),
                                   mesh.node(int(cell[2]-1)), mesh.node(int(cell[3]-1)),
                                   marker=int(cell[4]))

    mesh.createNeighbourInfos()

    for bound in bounds:
        if dim == 2:
            mesh.createEdge(mesh.node(int(bound[0]-1)), mesh.node(int(bound[1]-1)), marker=int(bound[2]))
        else:
            mesh.createTriangleFace(mesh.node(int(bound[0]-1)), mesh.node(int(bound[1]-1)),
                                    mesh.node(int(bound[2]-1)), marker=int(bound[3]))

    # assign marker to corresponding nodes (sensors, reference nodes, etc.)
    if len(points) > 0:
        for point in points:
            mesh.node(point[0]-1).setMarker(-point[1])

    if verbose:
        points = np.asarray(points)
        node_types = np.unique(points[:,1])
        print '  Marked nodes: %s ' % len(points) + str(tuple(node_types))
        print '\nDone. \n'
        print '  ' + str(mesh)

    return mesh

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print "Reads Gmsh ASCII file (*.msh) and returns binary mesh (*.bms) for BERT + VTK file. \n"
        print "Usage: python %s filename" % sys.argv[0]
        print readGmsh.func_doc[-688:]
    else:
        t_0 = time.time()
        mesh = readGmsh(sys.argv[1], verbose=True)
        mesh.save(sys.argv[1][:-4])
        mesh.exportVTK(sys.argv[1][:-4])
        elapsed = time.time() - t_0
        m, s = divmod(elapsed, 60)
        h, m = divmod(m, 60)
        print "\nElapsed run time: %02dh %02dmin %02dsec" % (h, m, s)
