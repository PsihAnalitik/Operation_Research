import numpy as np
from Symplex_method import _SymplexTable

'''
z = np.array([2, -3, 1])
A = np.array([[1, 2, 1], 
              [-1, 1, 3], 
              [0, -1, 2]])
b = np.array([5, 1, 4])

problem = {'A': A, 'c': z, 'b': b, 'type': 'min', 'signs': ['=', '=', '=']}
symplex_ex = _SymplexTable(problem)
res = symplex_ex.iterate()
'''

z = np.array([8, -2, 4, -7, 4], dtype=float)
A = np.array([[4, -7, 6, 5, -2],
              [9, 2, -9, 3, 9],
              [-1, 5, 2, -4, 5]], dtype=float)
b = np.array([9, 6, 6], dtype=float)

signs = np.array(['=', '=', '='])

problem = {'A': A, 'c': z, 'b': b, 'type': 'min', 'signs': ['=', '=', '=']}
symplex_ex = _SymplexTable(problem)
res = symplex_ex.iterate()
print(res['basis'])
print(res['table'])