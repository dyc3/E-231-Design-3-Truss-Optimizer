#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt

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

# trusses must span 15 inches, and there must be a connection at the top center of the truss
def is_valid(truss):
	pass

test_truss.draw()
