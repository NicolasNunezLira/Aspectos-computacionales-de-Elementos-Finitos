"""
Programa que convierte una malla, en formato XML,
generada por dolfin-convert y que inicialmente
estaba hecha con un archivo msh generado por GMSH

python convert.py --name nombre_de_la_malla

ejemplo:

python convert.py --name stenosis

crea la malla stenosis.h5 sabiendo que la malla
fue creada con gmsh

"""

from dolfin import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--name', type=str, default="toto")
a = parser.parse_args()

mesh       = Mesh(a.name+".xml");
domains    = MeshFunction('size_t',mesh,a.name+"_physical_region.xml");
bc_markers = MeshFunction('size_t',mesh,a.name+"_facet_region.xml");
hdf        = HDF5File(mesh.mpi_comm(), a.name+".h5", "w")

hdf.write(mesh, "/mesh")
hdf.write(domains, "/domains")
hdf.write(bc_markers, "/bc_markers")
