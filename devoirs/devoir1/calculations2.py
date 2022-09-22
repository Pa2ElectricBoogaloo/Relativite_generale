#%%
from sympy import *
from sympy import Array
from IPython.display import display

init_printing()

r, phi = symbols(r'r \varphi')
coords = Array([r,phi])
X =  Array([2 * r * cos(phi), 2 * r * sin(phi), r**2 - 1])/(r**2 + 1)
dX = derive_by_array(X, coords).simplify()

# Metrique 
gd = tensorcontraction(tensorproduct(dX, dX), (1, 3)).simplify()

# Metrique inverse
gu = Array([[1/gd[0, 0], 0], [0, 1/gd[1, 1]]]) 

# Symbole Christoffel 
dg = derive_by_array(gd, coords)

gammad = (permutedims(dg,(2,0,1)) + permutedims(dg,(2,1,0))-dg)/2

gamma = tensorcontraction(tensorproduct(gu, gammad), (1,2)).simplify()

# Tenseur de Riemann
dgamma = derive_by_array(gamma, coords).simplify()

gg = tensorcontraction(tensorproduct(gamma, gamma), (0,4)).simplify()

Riemann = (
    permutedims(dgamma,(1,2,0,3)) 
    - permutedims(dgamma,(1,2,3,0)) 
    + permutedims(gg,(2,0,3,1)) 
    - permutedims(gg,(2,0,1,3))
).simplify()


# Tenseur de Ricci
Ricci = tensorcontraction(Riemann, (0, 2)).simplify()

# Scalaire de Ricci 
R = tensorcontraction(tensorproduct(gu, Ricci), (0, 2), (1, 3)).simplify()

# Display 
display(gd, gu, gamma, Riemann, Ricci, R)

# %%
