# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 23:11:04 2024

@author: mmanc
"""

# =============================================================================
# Schwarzchild metric and Christoffel symbol
# =============================================================================

import sympy
from einsteinpy.symbolic import MetricTensor, ChristoffelSymbols, RiemannCurvatureTensor

sympy.init_printing() 
#%%

syms = sympy.symbols('t r θ φ')
G, M, c, a = sympy.symbols("G M c a")
#where a is the schwarzchild radius i.e. a = 2GM/c^2
list2d = [[0 for i in range(4)] for i in range(4)]
list2d[0][0] = 1 - (a / syms[1])
list2d[1][1] = -1 / (1 - (a / syms[1]))
list2d[2][2] = -1 * (syms[1] ** 2) 
list2d[3][3] = -1 * (syms[1] ** 2) * (sympy.sin(syms[2]) ** 2) 
sch = MetricTensor(list2d, syms)
sch.tensor()
print(sch)

#%%

ch = ChristoffelSymbols.from_metric(sch)
ch.tensor()

# note 0 = r, 1 = theta and 2 = phi

christoffel_dict = {}

index_name = {'t':0, 'r':1, 'θ':2, 'φ':3}
print('The non-zero Christoffel symbols are:')
print('')
for i in index_name:
    for j in index_name:
        for k in index_name:
            symbol = ch.tensor()[index_name[i], index_name[j], index_name[k]]
            indices = [i,j,k]
            
            if symbol != 0 and symbol not in christoffel_dict.values():
                # Add symbol to dictionary
                christoffel_dict[str(indices)] = symbol
                print(indices, '=', symbol)
                print('')
     
            
print('')
print('The rest of the Christoffel symbols are zero')





