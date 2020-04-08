#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import math
from anastruct import SystemElements
import random
import copy
import numpy as np

# trusses must span 15 inches, and there must be a connection at the top center of the truss
# member length must not exceed 72 inches, as 2 lengths of 36 inches

# extra credit: +2
# holds 320 lbs, less than 60 in

# extra credit: +2
# 6:1 strength to weight ratio

MIN_WIDTH = 15
MAX_HEIGHT = 4

def dist(a, b):
	return math.sqrt((b[0] - a[0])**2 + (b[1] - a[1])**2)

def midpoint(a, b):
    return [(a[0]+b[0])/2, (a[1]+b[1])/2]

def valmap(value, istart, istop, ostart, ostop):
	return ostart + (ostop - ostart) * ((value - istart) / (istop - istart))

def lbsToN(lbs):
	return lbs * 4.4482216282509

class Truss:
	def __init__(self):
		# format:
		# list of tuples of x and y coordinates
		self.nodes = []
		self.members = []

	def draw(self):
		lines = list([list(self.nodes[idx] for idx in member) for member in self.members])
		for line in lines:
			plt.plot([line[0][0], line[1][0]], [line[0][1], line[1][1]])
		plt.show()

	def is_valid(self):
		left_most = min(self.nodes, key=lambda p: p[0])
		right_most = max(self.nodes, key=lambda p: p[0])
		top_most = max(self.nodes, key=lambda p: p[1])
		bottom_most = min(self.nodes, key=lambda p: p[1])
		total_width = right_most[0] - left_most[0]
		total_height = top_most[1] - bottom_most[1]
		lines = list([list(self.nodes[idx] for idx in member) for member in self.members])
		total_length = sum([dist(*line) for line in lines])
		return total_width >= MIN_WIDTH and total_height <= MAX_HEIGHT and total_length <= 72 and len(self.members) >=  2 * len(self.nodes) - 3

	def calculate_member_forces(self):
		top_most_idx = self.nodes.index(max(self.nodes, key=lambda p: p[1])) # load will be placed on this node
		pass

test_truss = Truss()
test_truss.nodes = [
	(0, 0),
	(15/2, 0),
	(15, 0),
	(15/2, 4),
]
test_truss.members = [
	[0, 1],
	[1, 2],
	[1, 3],
	[0, 3],
	[2, 3],
]

print("is_valid", test_truss.is_valid())
# test_truss.draw()

def generate_truss():
	"""
	Randomly generate a valid truss
	"""
	ss = SystemElements(EA=15000, EI=5000)
	width = MIN_WIDTH
	height = MAX_HEIGHT
	subdivide_mode = random.choice(["triangle_subdivide", "radial_subdivide", "pillar_subdivide"])
	if subdivide_mode == "triangle_subdivide":
		subdivides = random.randint(1, 2)
		triangles = [
			[
				[[0, 0], [width, 0]],
				[[width, 0], [width/2, height]],
				[[width/2, height], [0, 0]],
			],
		]
		for _ in range(subdivides):
			new_triangles = []
			for triangle in triangles:
				mids = [midpoint(*line) for line in triangle]
				new_triangles += [
					[
						[triangle[0][0], mids[0]],
						[mids[0], mids[2]],
						[mids[2], triangle[0][0]],
					],
					[
						[mids[2], mids[1]],
						[mids[1], triangle[2][0]],
						[triangle[2][0], mids[2]],
					],
					[
						[mids[0], triangle[1][0]],
						[triangle[1][0], mids[1]],
						[mids[1], mids[0]],
					],
					[
						[mids[2], mids[0]],
						[mids[0], mids[1]],
						[mids[1], mids[2]],
					],
				]
			triangles = new_triangles
		# TODO: deduplicate lines
		for triangle in triangles:
			for line in triangle:
				ss.add_truss_element(location=line)
	elif subdivide_mode == "radial_subdivide":
		subdivides = random.randint(1, 4)
		step_size = width / 2 / subdivides
		bottom_midpoint = midpoint([0, 0], [width, 0])
		lines = []
		for x in np.arange(0, width + 0.1, step_size):
			lines += [
				[bottom_midpoint, [x, valmap(x, 0, width / 2, 0, height) if x <= width / 2 else valmap(x, width / 2, width, height, 0)]],
			]
		lines[-1][1][1] = 0 # HACK: set last y value to 0
		top_points = [p[1] for p in lines]
		top_lines = []
		for i in range(1, len(top_points)):
			top_lines += [
				[top_points[i - 1], top_points[i]]
			]
		lines += top_lines
		for line in lines:
			ss.add_truss_element(location=line)
	elif subdivide_mode == "pillar_subdivide":
		subdivides = random.randint(1, 4)
		step_size = width / 2 / subdivides
		lines = []
		for x in np.arange(step_size, width, step_size):
			lines += [
				[[x, 0], [x, valmap(x, 0, width / 2, 0, height) if x <= width / 2 else valmap(x, width / 2, width, height, 0)]],
			]
		top_points = [p[1] for p in lines]
		edge_lines = []
		for i in range(1, len(top_points)):
			edge_lines += [
				[top_points[i - 1], top_points[i]],
				[[top_points[i - 1][0], 0], [top_points[i][0], 0]],
			]
			if i < len(top_points) / 2:
				edge_lines += [
					[[top_points[i - 1][0], 0], top_points[i]],
				]
			else:
				edge_lines += [
					[top_points[i - 1], [top_points[i][0], 0]],
				]
		lines += [
			[[0, 0], top_points[0]],
			[[0, 0], [top_points[0][0], 0]],
			[[width, 0], top_points[-1]],
			[[width, 0], [top_points[-1][0], 0]],
		]
		lines += edge_lines
		for line in lines:
			ss.add_truss_element(location=line)

	ss.add_support_fixed(node_id=ss.find_node_id(vertex=[0, 0]))
	ss.add_support_fixed(node_id=ss.find_node_id(vertex=[width, 0]))
	return ss

ss = generate_truss()
ss.point_load(Fy=-500, node_id=ss.find_node_id(vertex=[MIN_WIDTH/2, MAX_HEIGHT]))
ss.solve(max_iter=500, geometrical_non_linear=True)
print(ss.get_node_results_system(node_id=ss.find_node_id(vertex=[MIN_WIDTH/2, 0])))

# ss.show_structure()
# ss.show_reaction_force()
# ss.show_axial_force()
# ss.show_shear_force()
# ss.show_bending_moment()
# ss.show_displacement()
