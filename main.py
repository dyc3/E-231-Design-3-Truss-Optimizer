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
# must hold >300 lbs but no more than 500

# extra credit: +2
# holds 320 lbs, less than 60 in

# extra credit: +2
# 6:1 strength to weight ratio

MIN_WIDTH = 15
MAX_HEIGHT = 4
MAX_POSSIBLE_LOAD = 500 # lbs
MIN_POSSIBLE_LOAD = 300 # lbs

MODULUS_OF_ELASTICITY = 15900000 # psi
BRASS_YIELD_STRESS = 59000 # psi
BRASS_CROSS_SECTION_AREA = 0.006216 # in^2
BRASS_DENSITY = 0.308 # lbs/in^3
MOMENT_OF_INERTIA = 1.2968e-05

JOHNSON_EULER_TRANSITION_lENGTH = 3.3 # in
END_CONDITION_FACTOR = 0.8 # in

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

def generate_truss(subdivide_mode=None, subdivides=None):
	"""
	Randomly generate a valid truss
	"""
	ss = SystemElements(EA=MODULUS_OF_ELASTICITY * BRASS_CROSS_SECTION_AREA, EI=MODULUS_OF_ELASTICITY * MOMENT_OF_INERTIA)
	width = MIN_WIDTH
	height = MAX_HEIGHT
	if not subdivide_mode:
		subdivide_mode = random.choice(["triangle_subdivide", "radial_subdivide", "pillar_subdivide"])
	if subdivide_mode == "triangle_subdivide":
		if not subdivides:
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
		raw_lines = np.reshape(triangles, (-1, 2, 2))
		# sort coordinates in each line
		raw_lines = [sorted(line, key=lambda p: p[0]) for line in raw_lines]
		# sort lines by first point's x value
		raw_lines = sorted(raw_lines, key=lambda l: l[0][0])
		# remove duplicate lines
		lines = []
		for line in raw_lines:
			is_duplicate = False
			for l in lines:
				if np.array_equal(line, l):
					is_duplicate = True
			if not is_duplicate:
				lines.append(line)
		for line in lines:
			ss.add_truss_element(location=line)
	elif subdivide_mode == "radial_subdivide":
		if not subdivides:
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
		if not subdivides:
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

ss = generate_truss("radial_subdivide", 2)
# ss.point_load(Fy=-500, node_id=ss.find_node_id(vertex=[MIN_WIDTH/2, MAX_HEIGHT]))
# ss.solve(max_iter=500, geometrical_non_linear=True)

# ss.show_structure()
# ss.show_reaction_force()
# ss.show_axial_force()
# ss.show_shear_force()
# ss.show_bending_moment()
# ss.show_displacement()

def is_truss_valid(truss):
	return len(truss.element_map.values()) == 2 * len(truss.node_map.values()) - 3

def check_for_failing_members(truss):
	def calculate_max_force(member):
		force = member['N']

		if force > 0:
			if member['length'] < JOHNSON_EULER_TRANSITION_lENGTH:
				# perform johnson calculation
				max_load = (
					BRASS_CROSS_SECTION_AREA *
					(BRASS_YIELD_STRESS - (MODULUS_OF_ELASTICITY ** -1 *
					(((BRASS_YIELD_STRESS / (2 * math.pi)) ** 2) *
					(END_CONDITION_FACTOR * member['length'] /
					math.sqrt(MOMENT_OF_INERTIA / BRASS_CROSS_SECTION_AREA)) ** 2
					))))

			else:
				# perfrom euler calculation
				max_load = (
					(math.pi ** 2 * MODULUS_OF_ELASTICITY * MOMENT_OF_INERTIA) /
					(END_CONDITION_FACTOR * member['length']) ** 2
				)

			return force > max_load
		else:
			return False

	return map(calculate_max_force, truss.get_element_results())

def score_truss(truss):
	member_lengths = [element.l for element in truss.element_map.values()]
	total_member_length = sum(member_lengths)
	material_weight = BRASS_CROSS_SECTION_AREA * total_member_length * BRASS_DENSITY

	load_node_id = truss.find_node_id(vertex=[MIN_WIDTH/2, MAX_HEIGHT])
	load_range_min, load_range_max = MIN_POSSIBLE_LOAD, MAX_POSSIBLE_LOAD
	max_load = 0
	while load_range_max - load_range_min > 0.001:
		mid = (load_range_min + load_range_max) / 2
		truss.point_load(Fy=-mid, node_id=load_node_id)
		truss.solve(max_iter=500)
		if not any(check_for_failing_members(truss)):
			load_range_min = mid
			max_load = mid
		else:
			load_range_max = mid
	print(f"all members: {total_member_length} in, {material_weight} lbs, holds max load {max_load}")
	return max_load / material_weight

print(f"truss score: {score_truss(ss)}")
