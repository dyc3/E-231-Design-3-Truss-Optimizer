#!/usr/bin/env python3
from anastruct import SystemElements
import math

width = 15
height = 4

ss = SystemElements(EA=15000, EI=5000)
ss.add_truss_element(location=[[0, 0], [0, height]])
ss.add_truss_element(location=[[0, 0], [width, 0]])
ss.add_truss_element(location=[[0, 4], [width, height]])
ss.add_truss_element(location=[[width, height], [width, 0]])
ss.add_truss_element(location=[[0, height], [width/4, 0]])
ss.add_truss_element(location=[[width/4*2, height], [width/4, 0]])
ss.add_truss_element(location=[[width/4*2, height], [width/4*3, 0]])
ss.add_truss_element(location=[[width, height], [width/4*3, 0]])

ss.add_support_fixed(node_id=ss.find_node_id(vertex=[0, 0]))
ss.add_support_fixed(node_id=ss.find_node_id(vertex=[15, 0]))

ss.point_load(Fy=-300, node_id=ss.find_node_id(vertex=[width/2, height]))

ss.solve()

# Get visual results.
ss.show_structure()
# ss.show_reaction_force()
# ss.show_axial_force()
# ss.show_shear_force()
# ss.show_bending_moment()
# ss.show_displacement()
