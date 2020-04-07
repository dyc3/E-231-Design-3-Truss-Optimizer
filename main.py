#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import math
from anastruct import SystemElements
import random
import copy

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
	# ss.add_truss_element(location=[[0, 0], [width, 0]])
	# ss.add_truss_element(location=[[0, 0], [width/2, height]])
	# ss.add_truss_element(location=[[width, 0], [width/2, height]])
	subdivide_mode = random.choice(["triangle_subdivide", "radial_subdivide", "pillar_subdivide"])
	subdivide_mode = "triangle_subdivide"
	if subdivide_mode == "triangle_subdivide":
		subdivides = random.randint(1, 3)
		triangles = [
			[
				[[0, 0], [width, 0]],
				[[width, 0], [width/2, height]],
				[[width/2, height], [0, 0]],
			],
		]
		# lines = copy.deepcopy(base_lines)
		# for line in lines:
		# 	ss.add_truss_element(location=line)
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

	# ss.add_support_fixed(node_id=ss.find_node_id(vertex=[0, 0]))
	# ss.add_support_fixed(node_id=ss.find_node_id(vertex=[width, 0]))
	return ss

generate_truss().show_structure()
