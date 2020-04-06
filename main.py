#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt

# trusses must span 15 inches, and there must be a connection at the top center of the truss

# extra credit: +2
# holds 320 lbs, less than 60 in

# extra credit: +2
# 6:1 strength to weight ratio

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
		return len(self.members) >=  2 * len(self.nodes) - 3

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
