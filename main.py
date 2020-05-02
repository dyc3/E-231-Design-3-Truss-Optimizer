#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import math
from anastruct import SystemElements
import random
import copy
import numpy as np
from tqdm import tqdm
import os, pickle
from itertools import combinations
from scipy.spatial.distance import euclidean
from multiprocessing import Pool
import argparse
from scipy.io import loadmat, savemat
import os
import hyperopt

# trusses must span 15 inches, and there must be a connection at the top center of the truss
# member length must not exceed 72 inches, as 2 lengths of 36 inches
# must hold >300 lbs but no more than 500

# extra credit: +2
# holds 320 lbs, less than 60 in

# extra credit: +2
# 6:1 strength to weight ratio

parser = argparse.ArgumentParser()
parser.add_argument('--generations', type=int, default=20)
parser.add_argument('--population-size', type=int, default=40)
parser.add_argument('--grid-size-x', type=int, default=3)
parser.add_argument('--grid-size-y', type=int, default=3)
parser.add_argument('--hyper-connected', default=False, action='store_true')
parser.add_argument('--cross-rate', type=float, default=0.5)
parser.add_argument('--mutation-rate', type=float, default=0.008)
parser.add_argument('--show-trusses', default=False, action='store_true')
parser.add_argument('--disable-parallel', dest='parallel', action='store_false')
parser.add_argument('--score-truss', type=str)
parser.add_argument('--optimize-truss', type=str)
parser.add_argument('--optimize-truss-iterations', type=int, default=1000)
args = parser.parse_args()

print(args)

MIN_WIDTH = 15
MAX_HEIGHT = 4
MAX_POSSIBLE_LOAD = 500 # lbs
MIN_POSSIBLE_LOAD = 0 # lbs

MODULUS_OF_ELASTICITY = 15900000 # psi
BRASS_YIELD_STRESS = 59000 # psi
BRASS_CROSS_SECTION_AREA = 0.006216 # in^2
BRASS_DENSITY = 0.308 # lbs/in^3
MOMENT_OF_INERTIA = 1.2968e-05

JOHNSON_EULER_TRANSITION_LENGTH = 3.3 # in
END_CONDITION_FACTOR = 0.8 # in

def dist(a, b):
	return math.sqrt((b[0] - a[0])**2 + (b[1] - a[1])**2)

def midpoint(a, b):
    return [(a[0]+b[0])/2, (a[1]+b[1])/2]

def valmap(value, istart, istop, ostart, ostop):
	return ostart + (ostop - ostart) * ((value - istart) / (istop - istart))

def lbsToN(lbs):
	return lbs * 4.4482216282509

def collinear(x1, y1, x2, y2, x3, y3):
	return (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) == 0

def slope(x1, y1, x2, y2):
	if x2 - x1 == 0:
		return float('inf')
	return (y2 - y1) / (x2 - x1)

class Truss:
	def __init__(self):
		# format:
		# list of tuples of x and y coordinates
		self.nodes = []
		self.members = []

	@property
	def total_width(self):
		left_most = min(self.nodes, key=lambda p: p[0])
		right_most = max(self.nodes, key=lambda p: p[0])
		return right_most[0] - left_most[0]

	@property
	def total_height(self):
		top_most = max(self.nodes, key=lambda p: p[1])
		bottom_most = min(self.nodes, key=lambda p: p[1])
		return top_most[1] - bottom_most[1]

	@property
	def top_node(self):
		return max(truss.nodes, key=lambda p: p[1])

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
		total_width = self.total_width
		total_height = self.total_height
		lines = list([list(self.nodes[idx] for idx in member) for member in self.members])
		total_length = sum([dist(*line) for line in lines])
		return total_width >= MIN_WIDTH and total_height <= MAX_HEIGHT and total_length <= 72 and len(self.members) >=  2 * len(self.nodes) - 3

	def calculate_member_forces(self):
		top_most_idx = self.nodes.index(max(self.nodes, key=lambda p: p[1])) # load will be placed on this node
		pass

	def to_anastruct(self) -> SystemElements:
		lines = list([list(self.nodes[idx] for idx in member) for member in self.members])
		truss = SystemElements(EA=MODULUS_OF_ELASTICITY * BRASS_CROSS_SECTION_AREA, EI=MODULUS_OF_ELASTICITY * MOMENT_OF_INERTIA)
		for line in lines:
			truss.add_truss_element(line)
		truss.add_support_hinged(node_id=truss.find_node_id(vertex=list(min(self.nodes, key=lambda p: p[0]))))
		truss.add_support_hinged(node_id=truss.find_node_id(vertex=list(max(self.nodes, key=lambda p: p[0]))))
		return truss

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

# print("is_valid", test_truss.is_valid())
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
		lines = exclude_duplicate_members(raw_lines)
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

	ss.add_support_hinged(node_id=ss.find_node_id(vertex=[0, 0]))
	ss.add_support_hinged(node_id=ss.find_node_id(vertex=[width, 0]))
	return ss

def generate_truss_grid(height, width, grid_size_x, grid_size_y, hyper_connected=False, exclude_known_useless=True):
	all_grid_points = np.array(np.meshgrid(np.arange(0, width + 0.01, width / grid_size_x), np.arange(0, height + 0.01, height / grid_size_y))).T.reshape(-1, 2)
	max_dist = euclidean([0, 0], [width / grid_size_x, height / grid_size_y])
	if exclude_known_useless:
		all_grid_points = np.array(list(filter(lambda point: point[1] < (height / grid_size_y * (grid_size_y - 1)) + 0.1 or point[0] > 1, all_grid_points)))
	all_possible_members = []
	if hyper_connected:
		for point1 in all_grid_points:
			for point2 in all_grid_points:
				if euclidean(point1, point2) < 1e-3:
					continue
				all_possible_members.append([point1, point2])
		all_possible_members = exclude_duplicate_members(all_possible_members)
	else:
		def is_adjacent(a, b):
			if a[0] == b[0]:
				return abs(a[1] - b[1]) <= height / grid_size_y + 0.01
			if a[1] == b[1]:
				return abs(a[0] - b[0]) <= width / grid_size_x + 0.01
			return euclidean(a, b) <= max_dist
		comb = np.array(list(filter(lambda x: is_adjacent(all_grid_points[x[1]], all_grid_points[x[0]]), combinations(range(len(all_grid_points)), 2))))
		all_possible_members = list(map(lambda idx: [all_grid_points[idx[0]], all_grid_points[idx[1]]], comb))

	# verify there are no zero length members, and no duplicate members
	for member in all_possible_members:
		assert euclidean(*member) > 0
		count = 0
		for memberb in all_possible_members:
			if are_members_equal(memberb, member):
				count += 1
			assert count <= 1

	all_possible_members = np.array(all_possible_members)
	return all_possible_members # + valmap(np.random.rand(*all_possible_members.shape), 0, 1, -0.2, 0.2) # for debugging

def are_members_equal(a, b):
	return np.allclose(sorted(a, key=lambda p: p[0]), sorted(b, key=lambda p: p[0]))

assert are_members_equal([[0, 1], [1, 0]], [[1, 0], [0, 1]])

def are_members_connected(a, b):
	for pointA in a:
		for pointB in b:
			if np.array_equal(pointA, pointB):
				return pointA
	return False

assert np.array_equal(are_members_connected([[0, 1], [1, 0]], [[1, 0], [2, 1]]), [1, 0])
assert not are_members_connected([[0, 1], [1, 0]], [[0, 0], [2, 1]])

def do_members_intersect(a, b):
	"""
	Returns True if the members intersect
	a1: [x, y] a point on the first line
	a2: [x, y] another point on the first line
	b1: [x, y] a point on the second line
	b2: [x, y] another point on the second line
	"""
	s = np.vstack([a[0],a[1],b[0],b[1]])        # s for stacked
	h = np.hstack((s, np.ones((4, 1)))) # h for homogeneous
	l1 = np.cross(h[0], h[1])           # get first line
	l2 = np.cross(h[2], h[3])           # get second line
	x, y, z = np.cross(l1, l2)          # point of intersection
	# if z == 0, lines are parallel
	return z != 0

def optimize_2_connection_nodes(grid, organism):
	"""
	Takes nodes that have 2 connections, and directly connects the end nodes
	Reduces the number of nodes used, and makes it easier to find duplicate members.
	"""
	assert len(grid) == len(organism)
	members = grid[organism]
	# because im lazy, let anastruct do the thinking
	truss = SystemElements()
	for member in members:
		truss.add_truss_element(member)
	assert len(members) == len(truss.element_map.values())
	two_connection_nodes = [truss.node_map[node_id] for node_id, connections in truss.node_element_map.items() if len(connections) == 2]

	for node in two_connection_nodes:
		elements = list(node.elements.values())
		connected_nodes = [
			elements[0].node_id1,
			elements[0].node_id2,
			elements[1].node_id1,
			elements[1].node_id2,
		]
		shortcut_member_nodes = [truss.node_map[node_id] for node_id in connected_nodes if node_id != node.id]
		shortcut_member = [node.vertex.coordinates for node in shortcut_member_nodes]
		if euclidean(*shortcut_member) < 1e-3:
			continue
		old_member_idx1 = -1
		old_member_idx2 = -1
		new_member_idx = -1
		for i, member in enumerate(grid):
			if old_member_idx1 >= 0 and old_member_idx2 >= 0 and new_member_idx >= 0:
				break
			if are_members_equal(shortcut_member, member):
				new_member_idx = i
				continue
			if not organism[i]:
				continue
			if are_members_equal([elements[0].vertex_1.coordinates, elements[0].vertex_2.coordinates], member):
				old_member_idx1 = i
				continue
			if are_members_equal([elements[1].vertex_1.coordinates, elements[1].vertex_2.coordinates], member):
				old_member_idx2 = i
				continue
		if old_member_idx1 == old_member_idx2:
			continue
		if new_member_idx == -1:
			raise Exception("member not found")
		organism[old_member_idx1] = False
		organism[old_member_idx2] = False
		organism[new_member_idx] = True
	return organism

def exclude_duplicate_members(members):
	"""
	Remove duplicate members
	"""
	lines = []
	for line in members:
		is_duplicate = False
		for l in lines:
			if are_members_equal(line, l):
				is_duplicate = True
				break
		if not is_duplicate:
			lines.append(line)
	return lines

def optimize_colinear_members(members):
	# members_by_slope = {}
	# for i, member in enumerate(members):
	# 	s = round(slope(*member.flatten()), 4)
	# 	if s not in members_by_slope:
	# 		members_by_slope[s] = []
	# 	members_by_slope[s].append((i, member))

	# for s in members_by_slope:
	# 	member_group = members_by_slope[s]
	# 	for a_i, a in member_group:
	# 		for b_i, b in member_group:
	# 			if are_members_equal(a, b):
	# 				continue
	# 			shared_point = are_members_connected(a, b)
	# 			if not np.any(shared_point):
	# 				continue
	# 			new_member = np.array(list(filter(lambda x: not np.array_equal(x, shared_point), list(a) + list(b))))
	# 			if new_member.shape == (2, 2):
	# 				members[a_i] = new_member
	# 				members[b_i] = None

	# return np.array(list(filter(lambda x: not np.any(np.isnan(x)), members)))
	# return members[~np.isnan(members)].reshape(-1, 2, 2)

	# for a_i, a in enumerate(members):
	# 	for b_i, b in enumerate(members):
	# 		if are_members_equal(a, b):
	# 			continue
	# 		shared_point = are_members_connected(a, b)
	# 		# make sure no other members are connected to this point
	# 		num_connections = 0
	# 		for c in members.reshape(-1, 2):
	# 			if np.array_equal(c, shared_point):
	# 				num_connections += 1
	# 			if num_connections > 2:
	# 				break
	# 		if num_connections != 2:
	# 			continue
	# 		new_member = np.array(list(filter(lambda x: not np.array_equal(x, shared_point), list(a) + list(b))))
	# 		if new_member.shape == (2, 2):
	# 			members[a_i] = new_member
	# 			members[b_i] = None

	# # return np.array(list(filter(lambda x: not np.any(np.isnan(x)), members)))
	# return members[~np.isnan(members)].reshape(-1, 2, 2)

	new_members = []
	for a_i, a in enumerate(members):
		for b_i, b in enumerate(members):
			if are_members_equal(a, b):
				continue
			shared_point = are_members_connected(a, b)
			if np.array_equal(shared_point, [0, 0]) or np.array_equal(shared_point, [MIN_WIDTH, 0]) or np.array_equal(shared_point, [MIN_WIDTH / 2, MAX_HEIGHT]):
				continue
			# make sure no other members are connected to this point
			num_connections = 0
			for c in members.reshape(-1, 2):
				if np.array_equal(c, shared_point):
					num_connections += 1
				if num_connections > 2:
					break
			if num_connections != 2:
				continue
			new_member = np.array(list(filter(lambda x: not np.array_equal(x, shared_point), list(a) + list(b))))
			if new_member.shape == (2, 2):
				members[a_i] = None
				members[b_i] = None
				new_members.append(new_member)

	# return np.array(list(filter(lambda x: not np.any(np.isnan(x)), members)))
	if len(new_members) > 0:
		return np.concatenate((members[~np.isnan(members)].reshape(-1, 2, 2), np.array(new_members)), axis=0)
	else:
		return members

# _test_optimize = optimize_colinear_members(np.array([ [[0, 0], [1, 0]], [[2, 0], [1, 0]] ]))
_test_optimize = optimize_colinear_members(np.array([ [[0, 0], [1, 0]], [[1, 0], [2, 0]] ], dtype="float64"))
assert len(_test_optimize) == 1

def generate_truss_by_grid(grid, enabled):
	"""
	enabled is a list of booleans indicating which members in the grid are enabled. Length must match the total possible members in the grid
	"""
	enabled = np.array(enabled)
	width = MIN_WIDTH / 2
	height = MAX_HEIGHT
	all_possible_members = grid
	# print(f"number of possible members: {len(all_possible_members)}")
	# assert len(all_possible_members) == len(enabled)
	members = all_possible_members[enabled]
	# print(f"members selected: {len(members)}")
	# members = optimize_colinear_members(members)
	# mirror the members to the right side
	members_mirror = np.copy(members)
	members_x, members_y = members_mirror[:,:,[True, False]], members_mirror[:,:,[False, True]]
	members_x *= -1
	members_x += width * 2
	members_mirror = np.append(members_x, members_y, axis=2)
	members = np.append(members, members_mirror, axis=0)
	members = exclude_duplicate_members(members)
	truss = SystemElements(EA=MODULUS_OF_ELASTICITY * BRASS_CROSS_SECTION_AREA, EI=MODULUS_OF_ELASTICITY * MOMENT_OF_INERTIA)
	for member in members:
		truss.add_truss_element(member)
	try:
		truss.add_support_hinged(node_id=truss.find_node_id(vertex=[0, 0]))
		truss.add_support_hinged(node_id=truss.find_node_id(vertex=[width * 2, 0]))
		return truss
	except:
		return None

def load_truss_from_file(mat_file_path):
	"""
	Loads a truss from the Stevens' Truss Analyzer program
	"""
	data = loadmat(mat_file_path)
	nodes = np.array(data['nodes']) - [1, 3]
	elements = np.array(data['elements']) - 1
	truss = Truss()
	truss.nodes = nodes
	truss.members = elements
	return truss

def load_truss_from_file_to_anastruct(mat_file_path):
	"""
	Loads a truss from the Stevens' Truss Analyzer program
	"""
	data = loadmat(mat_file_path)
	nodes = np.array(data['nodes']) - [1, 3]
	elements = np.array(data['elements']) - 1
	truss = SystemElements(EA=MODULUS_OF_ELASTICITY * BRASS_CROSS_SECTION_AREA, EI=MODULUS_OF_ELASTICITY * MOMENT_OF_INERTIA)
	for elementnodes in elements:
		member = nodes[elementnodes]
		truss.add_truss_element(member)
	truss.add_support_hinged(node_id=truss.find_node_id(vertex=[0, 0]))
	truss.add_support_hinged(node_id=truss.find_node_id(vertex=[MIN_WIDTH, 0]))
	return truss

def save_truss_for_truss_analyzer(truss, math_file_path):
	left_node_id = truss.find_node_id(vertex=[0, 0])
	right_node_id = truss.find_node_id(vertex=[MIN_WIDTH, 0])
	nodes = [None] * len(truss.node_map)
	elements = [None] * len(truss.element_map)
	for node_id in truss.node_map:
		nodes[node_id - 1] = truss.node_map[node_id].vertex.coordinates
	if left_node_id > 2:
		nodes[0], nodes[left_node_id - 1] = nodes[left_node_id - 1], nodes[0]
	if right_node_id > 2:
		nodes[1], nodes[right_node_id - 1] = nodes[right_node_id - 1], nodes[1]
	nodes = np.array(nodes)
	for element_id in truss.element_map:
		point1 = truss.element_map[element_id].node_1.vertex.coordinates
		point2 = truss.element_map[element_id].node_2.vertex.coordinates
		point1idx = -1
		point2idx = -1
		for i, node in enumerate(nodes):
			if np.array_equal(node, point1):
				point1idx = i
				break
		for i, node in enumerate(nodes):
			if np.array_equal(node, point2):
				point2idx = i
				break
		elements[element_id - 1] = [point1idx, point2idx]
	savemat(math_file_path, {
		"nodes": np.array(nodes) + [1, 3],
		"elements": np.array(elements) + 1,
		"loads": []
	}, appendmat=False)

# ss = generate_truss("radial_subdivide", 2)
# ss.point_load(Fy=-500, node_id=ss.find_node_id(vertex=[MIN_WIDTH/2, MAX_HEIGHT]))
# ss.solve(max_iter=500, geometrical_non_linear=True)

# ss.show_structure()
# ss.show_reaction_force()
# ss.show_axial_force()
# ss.show_shear_force()
# ss.show_bending_moment()
# ss.show_displacement()

def is_truss_valid(truss):
	# for node_id in truss.node_element_map:
	# 	if len(truss.node_element_map[node_id]) < 2:
	# 		return False
	return len(truss.element_map.values()) == 2 * len(truss.node_map.values()) - 3

def calculate_max_force(member):
	force = -member['N']

	if force > 0:
		if member['length'] < JOHNSON_EULER_TRANSITION_LENGTH:
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

		return max_load / force
	else:
		return False

def check_max_load(truss):
	loads = list(map(calculate_max_force, truss.get_element_results()))
	return min(list(filter(lambda x: x > 0, loads)) or [0])

def score_truss(truss, silent=False, load_node=None):
	member_lengths = [element.l for element in truss.element_map.values()]
	total_member_length = sum(member_lengths)
	material_weight = BRASS_CROSS_SECTION_AREA * total_member_length * BRASS_DENSITY
	num_hanging_members = sum([1 for connections in truss.node_element_map.values() if len(connections) == 1])

	load_node_id = truss.find_node_id(vertex=(list(load_node) if np.any(load_node) else [MIN_WIDTH / 2, MAX_HEIGHT]))
	if not load_node_id:
		load_node_id = truss.nearest_node("both", [MIN_WIDTH / 2, MAX_HEIGHT])
	load_range_min, load_range_max = MIN_POSSIBLE_LOAD, MAX_POSSIBLE_LOAD
	truss.point_load(Fy=-1, node_id=load_node_id)
	truss.solve(max_iter=500)
	max_load = check_max_load(truss)
	if not silent:
		print(f"all members: {total_member_length} in, {material_weight:.2f} lbs, holds max load {max_load}, {num_hanging_members} hanging members")
	return max_load / material_weight / total_member_length * ((total_member_length < 72) * 2 + 1)

# for mode in ["triangle_subdivide", "radial_subdivide", "pillar_subdivide"]:
# 	for subdivides in range(1, 5):
# 		truss = generate_truss(mode, subdivides)
# 		is_valid = is_truss_valid(truss)
# 		save_truss_for_truss_analyzer(truss, f"./trusses/truss_{mode}_{subdivides}.mat")
# 		score = score_truss(truss)
# 		print(f"truss {mode}/{subdivides} valid: {is_valid} score: {score:.1f}")

np.random.seed(42)

# truss = generate_truss_by_grid(grid, ([True, True, False] * 5000)[:len(grid)])
# truss.show_structure()

def generate_valid_truss(grid):
	truss = members = None
	attempt = 0
	while not truss or not is_truss_valid(truss):
		attempt += 1
		# print(f"{attempt}   ", end="\r")
		members = np.random.rand(len(grid)) < (0.06 if args.hyper_connected else 0.4)
		if args.hyper_connected:
			members = optimize_2_connection_nodes(grid, members)
		truss = generate_truss_by_grid(grid, members)
		# if truss:
		# 	truss.show_structure()
		if truss and not truss.find_node_id(vertex=[MIN_WIDTH / 2, MAX_HEIGHT]):
			truss = None
	return members

def mutate(pop, mutation_rate=args.mutation_rate):
	"""
	Vectorized random mutations.
	:param pop: (array)
	:param mutation_rate: (flt)
	:return: (array)
	"""
	idx = np.where(np.random.rand(pop.shape[0], pop.shape[1]) < mutation_rate)
	val = np.random.randint(0, 2, idx[0].shape[0])
	pop[idx] = val
	return pop

def crossover(pop, cross_rate=args.cross_rate):
	"""
	Vectorized crossover
	:param pop: (array)
	:param cross_rate: (flt)
	:return: (array)
	"""
	# [bool] Rows that will crossover.
	selection_rows = np.random.rand(pop.shape[0]) < cross_rate

	selection = pop[selection_rows]
	shuffle_seed = np.arange(selection.shape[0])
	np.random.shuffle(shuffle_seed)

	# 2d array with [rows of the (selected) population, bool]
	cross_idx = np.array(np.round(np.random.rand(selection.shape[0], pop.shape[1])), dtype=np.bool)
	idx = np.where(cross_idx)

	selection[idx] = selection[shuffle_seed][idx]
	pop[selection_rows] = selection

	return pop

def rank_selection(pop, fitness):
	"""
	Rank selection. And make a selection based on their ranking score. Note that this isn't the fitness.
	:param pop: (array) Population.
	:param fitness: (array) Fitness values.
	:return: (array) Population selection with replacement, selected for mating.
	"""
	order = np.argsort(fitness)
	# Population ordered by fitness.
	pop = np.array(pop)[order]

	# Rank probability is proportional to you position, not you fitness. So an ordered fitness array, would have these
	# probabilities [1, 1/2, 1/3 ... 1/n] / sum
	rank_p = 1 / np.arange(1, pop.shape[0] + 1)
	# Make a selection based on their ranking.
	idx = np.random.choice(np.arange(pop.shape[0]), size=pop.shape[0], replace=True, p=rank_p / np.sum(rank_p))
	return pop[idx]

def eliminate_zero_force_members(organism):
	organism = copy.deepcopy(organism)
	truss = generate_truss_by_grid(grid, organism)
	if not truss:
		return organism
	if not truss.find_node_id(vertex=[MIN_WIDTH/2, MAX_HEIGHT]):
		return organism
	truss.point_load(Fy=-100, node_id=truss.find_node_id(vertex=[MIN_WIDTH/2, MAX_HEIGHT]))
	truss.solve()
	member_idxs = np.where(organism)[0]
	loads = np.array(list(map(lambda x: x["N"], truss.get_element_results())))
	assert len(member_idxs) * 2 == len(loads)
	loads = loads[:len(member_idxs)]
	zero_force_idxs = member_idxs[np.where((loads == 0) | (loads == False))[0]]
	organism[zero_force_idxs] = False
	return organism

name = "grid_4_6"

if not os.path.exists("trusses"):
	os.mkdir("trusses")
if not os.path.exists("img"):
	os.mkdir("img")
if not os.path.exists("img/" + name):
	os.mkdir("img/" + name)

def save_organism_figure(organism, fitness, generation, suffix=""):
	truss = generate_truss_by_grid(grid, organism)
	fig = truss.show_structure(show=False, verbosity=1)
	plt.title(f"fitness = {round(fitness, 3)}")
	fig.savefig(os.path.join("./img", name, f"ga{generation}{suffix}.png"))
	save_truss_for_truss_analyzer(truss, f"./trusses/gen_{generation}{suffix}.mat")

def genetic_optimization(population):
	for generation in range(args.generations):
		print(f"GENERATION {generation}")
		fitness = []
		for organism in tqdm(population):
			truss = generate_truss_by_grid(grid, organism)
			try:
				fitness.append(score_truss(truss, True))
			except np.linalg.LinAlgError:
				fitness.append(0)
			except:
				fitness.append(0)
		fitness = np.array(fitness)
		max_idx = np.argmax(fitness)

		try:
			print(f"fitness = {round(fitness[max_idx], 3)}")
			save_organism_figure(population[max_idx], fitness[max_idx], generation)
			# try:
			# 	save_organism_figure(eliminate_zero_force_members(population[max_idx]), fitness[max_idx], generation, "_nozero")
			# except np.linalg.LinAlgError:
			# 	path = os.path.join("./img", name, f"ga{generation}_nozero.png")
			# 	if os.path.exists(path):
			# 		os.remove(path)
			with open(os.path.join("./img", name, "save.pkl"), "wb") as f:
				pickle.dump(population, f)
		except AttributeError:
			pass

		pop = copy.deepcopy(rank_selection(population, fitness))
		valid_pop = []
		while len(valid_pop) < len(population):
			new_population = mutate(crossover(pop))
			for organism in new_population:
				truss = generate_truss_by_grid(grid, organism)
				if truss and truss.find_node_id(vertex=[MIN_WIDTH / 2, MAX_HEIGHT]) and is_truss_valid(truss):
					valid_pop.append(organism)
					if len(valid_pop) >= len(population):
						break
			print(f"\r mutating {len(valid_pop)}/{len(population)}   ", end="")
		print("done")
		population = np.array(valid_pop)

	return population

def get_index_of_node(nodes, targetnode):
	for i, node in enumerate(truss.nodes):
		if np.array_equal(node, targetnode):
			return i
	return -1

def optimize_member_lengths(truss: Truss) -> Truss:
	"""
	Moves the nodes around to find the optimal lengths. Keeps the support nodes in place, and does not move the load node's x axis.
	Uses hyperparameter optimizeation to maximize the output of the truss score
	"""
	truss.nodes = np.array(truss.nodes)
	load_node_idx = get_index_of_node(truss.nodes, max(truss.nodes, key=lambda p: p[1]))
	support_node_idxs = [
		get_index_of_node(truss.nodes, min(truss.nodes, key=lambda p: p[0])),
		get_index_of_node(truss.nodes, max(truss.nodes, key=lambda p: p[0])),
	]
	x_bounds = (0, MIN_WIDTH)
	y_bounds = (MAX_HEIGHT * 0.1, truss.nodes[load_node_idx][1] - 0.1)
	trials = hyperopt.Trials()
	problem_space = {}
	for i, node in enumerate(truss.nodes):
		if i in support_node_idxs:
			continue
		if i == load_node_idx:
			# problem_space[f"{i}:1"] = hyperopt.hp.uniform(f"{i}:1", *(y_bounds[1] * 0.85, y_bounds[1]))
			pass
		else:
			problem_space[f"{i}:0"] = hyperopt.hp.uniform(f"{i}:0", node[0] - 0.5, node[0] + 0.5)
			if node[1] != 0:
				problem_space[f"{i}:1"] = hyperopt.hp.uniform(f"{i}:1", *y_bounds)
	def truss_loss(parameters):
		tmp_truss = Truss()
		tmp_truss.nodes = np.copy(truss.nodes)
		tmp_truss.members = truss.members
		tmp_truss.nodes[load_node_idx] = truss.nodes[load_node_idx]
		for support_idx in support_node_idxs:
			tmp_truss.nodes[support_idx] = truss.nodes[support_idx]
		for key, value in parameters.items():
			node_idx, dimension_idx = map(int, key.split(":"))
			tmp_truss.nodes[node_idx][dimension_idx] = value
		return 1 / score_truss(tmp_truss.to_anastruct(), silent=True, load_node=[truss.nodes[load_node_idx][0], tmp_truss.nodes[load_node_idx][1]])
	best = hyperopt.fmin(fn=truss_loss, space=problem_space, algo=hyperopt.tpe.suggest, trials=trials, max_evals=args.optimize_truss_iterations)
	print(best)
	for key, value in best.items():
		node_idx, dimension_idx = map(int, key.split(":"))
		truss.nodes[node_idx][dimension_idx] = value
	return truss

if args.optimize_truss:
	truss = load_truss_from_file(args.optimize_truss)
	truss = optimize_member_lengths(truss)
	print(f"{score_truss(truss.to_anastruct(), load_node=truss.top_node)}")
	save_truss_for_truss_analyzer(truss.to_anastruct(), args.optimize_truss.replace(".mat", "") + "_optimized.mat")
	if args.show_trusses:
		truss.to_anastruct().show_structure()
	os._exit(0)

grid = generate_truss_grid(MAX_HEIGHT, MIN_WIDTH / 2, args.grid_size_x, args.grid_size_y, hyper_connected=args.hyper_connected)

if args.score_truss:
	if args.score_truss == "gui":
		import PySimpleGUI as sg
		from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, FigureCanvasAgg

		sg.theme('DarkAmber')   # Add a touch of color
		# All the stuff inside your window.
		layout = [  [sg.Canvas(size=(10, 5), key='canvas')],
					[sg.Checkbox(f'member {i}') for i in range(len(grid) // 2)],
					[sg.Checkbox(f'member {i}') for i in range(len(grid) // 2, len(grid))],
					[sg.Button('Ok'), sg.Button('Cancel')] ]

		def draw_figure(canvas, figure, loc=(0, 0)):
			figure_canvas_agg = FigureCanvasTkAgg(figure, window['canvas'].TKCanvas)
			figure_canvas_agg.draw()
			figure_canvas_agg.get_tk_widget().grid(row=0, column=0)
			return figure_canvas_agg

		# Create the Window
		window = sg.Window('Window Title', layout)
		window.Finalize()

		# Event Loop to process "events" and get the "values" of the inputs
		while True:
			print("loop")
			event, values = window.read()
			print(event, values)
			if event in (None, 'Cancel'):   # if user closes window or clicks cancel
				break
			plt.clf()

			truss = generate_truss_by_grid(grid, list([values[i] for i in range(len(grid))]))
			if truss:
				print(len(truss.element_map.values()))
				print(2 * len(truss.node_map.values()) - 3)
				fig = truss.show_structure(show=False)
				fig_canvas_agg = draw_figure(window['canvas'].TKCanvas, fig)
			else:
				print("couldnt build truss")
		window.close()
	else:
		truss = load_truss_from_file(args.score_truss)
	ana = truss.to_anastruct()
	print(f'valid: {is_truss_valid(ana)}')
	print(f'score: {score_truss(ana, load_node=truss.top_node)}')
	if args.show_trusses:
		ana.show_structure()
	os._exit(0)

print("generating initial population...")
truss_population = None
if args.parallel:
	with Pool() as p:
		truss_population = list(tqdm(p.imap(generate_valid_truss, [grid] * args.population_size), total=args.population_size))
else:
	truss_population = list(tqdm(map(generate_valid_truss, [grid] * args.population_size), total=args.population_size))

for members in genetic_optimization(truss_population):
	truss = generate_truss_by_grid(grid, members)
	try:
		print(f"truss score: {score_truss(truss)}")
		if args.show_trusses:
			truss.show_results()
	except np.linalg.LinAlgError:
		print(f"can't score invalid truss")
		if args.show_trusses:
			truss.show_structure()
	except Exception as e:
		print(f"can't score: {e}")
		if args.show_trusses:
			truss.show_structure()
