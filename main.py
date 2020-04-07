#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import math

# trusses must span 15 inches, and there must be a connection at the top center of the truss
# member length must not exceed 72 inches, as 2 lengths of 36 inches

# extra credit: +2
# holds 320 lbs, less than 60 in

# extra credit: +2
# 6:1 strength to weight ratio

def dist(a, b):
	return math.sqrt((b[0] - a[0])**2 + (b[1] - a[1])**2)

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
		return total_width >= 15 and total_height <= 4 and total_length <= 72 and len(self.members) >=  2 * len(self.nodes) - 3

	def calculate_member_forces(self):
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
test_truss.draw()
