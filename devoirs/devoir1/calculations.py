#%%
import sympy as sp 

theta, phi, psi = sp.symbols(r'\theta \phi \psi')
that = sp.Matrix([1, 0, 0])
xhat = sp.Matrix([0, 1, 0])

R = sp.Matrix([[1, 0, 0], [0, sp.cos(theta), -sp.sin(theta)], [0, sp.sin(theta), sp.cos(theta)]])

B1 = sp.Matrix([[sp.cosh(phi), sp.sinh(phi), 0], [sp.sinh(phi), sp.cosh(phi), 0], [0, 0, 1]])

B2 = sp.Matrix([[sp.cosh(psi), 0, sp.sinh(psi)], [0, 1, 0], [sp.sinh(psi), 0, sp.cosh(psi)]])

l =  sp.sqrt(sp.sinh(phi)**2 * sp.cosh(psi)**2 + sp.sinh(psi)**2)
c = sp.sinh(phi) * sp.cosh(psi)/l
s = sp.sinh(psi)/l

Rthat = sp.Matrix([[1, 0, 0], [0, c, s], [0, -s, c]])


s1 = sp.sqrt(sp.cosh(phi)**2 * sp.cosh(psi)**2 - 1)
c1 = sp.cosh(phi) * sp.cosh(psi)

B3 = sp.simplify(sp.Matrix([[c1, -s1, 0], [-s1, c1, 0], [0, 0, 1]]))

c2 = sp.simplify(Rthat**(-1) * B3 * Rthat * B1 * B2 * xhat)
sp.simplify(c2[2]/c2[1])


# %%
import sympy as sp 

theta, phi, psi = sp.symbols(r'\theta \phi \psi')

R = sp.Matrix([[1, 0, 0], [0, sp.cos(theta), -sp.sin(theta)], [0, sp.sin(theta), sp.cos(theta)]])

B1 = sp.Matrix([[sp.cosh(phi), sp.sinh(phi), 0], [sp.sinh(phi), sp.cosh(phi), 0], [0, 0, 1]])

B2 = sp.Matrix([[sp.cosh(psi), 0, sp.sinh(psi)], [0, 1, 0], [sp.sinh(psi), 0, sp.cosh(psi)]])

M = R * B1 * B2 
M[2, 1]-M[1, 2]

print(sp.latex(R), '\n')
print(sp.latex(B1), '\n')
print(sp.latex(B2), '\n')
print(sp.latex(M), '\n')
print(sp.latex(M[2, 1]-M[1, 2]))
# %%
