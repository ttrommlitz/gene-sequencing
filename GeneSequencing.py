#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		# map (i, j) to prev (i, j) along the optimal path
		prev = {}
		# n is the length of the shorter sequence
		smaller = seq2 if len(seq2) < len(seq1) else seq1
		larger = seq1 if len(seq1) > len(seq2) else seq2

		n = len(smaller)
		m = len(larger)

		if m > align_length:
			m = align_length
		
		if n > align_length:
			n = align_length

		if not banded:
			# initialize an m+1 by n+1 matrix
			dp = [[0 for i in range(m+1)] for j in range(n+1)]
			
			# initialize the first row and column

			prev[(0, 0)] = None

			for i in range(1, n+1):
				dp[i][0] = INDEL * i
				prev[(i, 0)] = (i-1, 0)

			for j in range(1, m+1):
				dp[0][j] = INDEL * j
				prev[(0, j)] = (0, j-1)

			# fill in the rest of the matrix
			for i in range(1, n + 1):
				for j in range(1, m + 1):
					left = dp[i][j-1] + INDEL
					top = dp[i-1][j] + INDEL
					diag = dp[i-1][j-1] + (MATCH if larger[j-1] == smaller[i-1] else SUB)

					# break ties in the following order: left, top, diag
					if left <= top and left <= diag:
						dp[i][j] = left
						prev[(i, j)] = (i, j-1)
					elif top <= left and top <= diag:
						dp[i][j] = top
						prev[(i, j)] = (i-1, j)
					else:
						dp[i][j] = diag
						prev[(i, j)] = (i-1, j-1)

		else:
			# banded algorithm with a band of 7
			d = 3
			dp = {}

			dp[(0,0)] = 0
			prev[(0, 0)] = None
			
			# initialize the first row and column
			for i in range(1, d + 1):
				dp[(i, 0)] = INDEL * i
				prev[(i, 0)] = (i-1, 0)
			
			for j in range(1, d + 1):
				dp[(0, j)] = INDEL * j
				prev[(0, j)] = (0, j-1)

			for i in range(1, n + 1):
				for j in range(i - d, i + d + 1):
					if j <= 0:
						continue

					if j > m:
						continue

					top = dp[(i-1, j)] + INDEL if (i-1, j) in dp else None
					left = dp[(i, j-1)] + INDEL if (i, j-1) in dp else None
					diag = dp[(i-1, j-1)] + (MATCH if larger[j-1] == smaller[i-1] else SUB)

					if j == i - d:
						if top <= diag:
							dp[(i, j)] = top
							prev[(i, j)] = (i-1, j)
						else:
							dp[(i, j)] = diag
							prev[(i, j)] = (i-1, j-1)

					elif j == i + d:
						if left <= diag:
							dp[(i, j)] = left
							prev[(i, j)] = (i, j-1)
						else:
							dp[(i, j)] = diag
							prev[(i, j)] = (i-1, j-1)
					
					else:
						# break ties in the following order: left, top, diag
						if left <= top and left <= diag:
							dp[(i, j)] = left
							prev[(i, j)] = (i, j-1)
						elif top <= left and top <= diag:
							dp[(i, j)] = top
							prev[(i, j)] = (i-1, j)
						else:
							dp[(i, j)] = diag
							prev[(i, j)] = (i-1, j-1)

		self.banded = banded
		self.MaxCharactersToAlign = align_length

		# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		score = None
		alignment1 = ''
		alignment2 = ''
		if not banded or (n, m) in dp:
			score = dp[n][m] if not banded else dp[(n, m)]

			i = n
			j = m

			while prev[(i, j)] is not None:
				prev_i, prev_j = prev[(i, j)]
				if prev_i == i - 1 and prev_j == j - 1:
					alignment1 += larger[j-1]
					alignment2 += smaller[i-1]
				elif prev_i == i - 1:
					alignment1 += '-'
					alignment2 += smaller[i-1]
				else:
					alignment1 += larger[j-1]
					alignment2 += '-'

				i = prev_i
				j = prev_j
			
			alignment1 = alignment1[::-1][:100]
			alignment2 = alignment2[::-1][:100]

		else:
			score = float('inf')
			alignment1 = 'No Alignment Possible'
			alignment2 = 'No Alignment Possible'

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
