import sympy as sp
from sympy.printing import ccode

# 1. Define Symbols
# n = current time step, p = previous time step
xn = sp.symbols('xn[0:4]')  # Current states: x1n, x2n, x3n, x4n
xp = sp.symbols('xp[0:4]')  # Previous states: x1p, x2p, x3p, x4p
T = sp.symbols('T')         # Sampling period
H, C, beta = sp.symbols('H C beta')

# 2. Continuous-time Equations (f(x)) from the paper
f = [
    xn[0] - xn[1] - C*xn[0]*xn[1] - xn[0]*xn[1]**2,
    H*xn[0] - 3*xn[1] + C*xn[0]*xn[1] + xn[0]*xn[1]**2 + beta*(xn[3] - xn[1]),
    xn[2] - xn[3] - C*xn[2]*xn[3] - xn[2]*xn[3]**2,
    H*xn[2] - 3*xn[3] + C*xn[2]*xn[3] + xn[2]*xn[3]**2 + 2*beta*(xn[1] - xn[3])
]

fp = [
    xp[0] - xp[1] - C*xp[0]*xp[1] - xp[0]*xp[1]**2,
    H*xp[0] - 3*xp[1] + C*xp[0]*xp[1] + xp[0]*xp[1]**2 + beta*(xp[3] - xp[1]),
    xp[2] - xp[3] - C*xp[2]*xp[3] - xp[2]*xp[3]**2,
    H*xp[2] - 3*xp[3] + C*xp[2]*xp[3] + xp[2]*xp[3]**2 + 2*beta*(xp[1] - xp[3])
]

# 3. Formulate the Tustin Residuals: R = x[n] - x[p] - (T/2)*(f(x[n]) + f(x[p]))
R = []
for i in range(4):
    R.append(xn[i] - xp[i] - (T/2) * (f[i] + fp[i]))

# 4. Calculate the Jacobian Matrix (J = dR/dxn)
J = sp.Matrix(R).jacobian(xn)

# 5. Generate C Code
print("// --- Tustin Residuals ---")
for i in range(4):
    print(f"R[{i}] = {ccode(R[i])};")

print("\n// --- Jacobian Matrix Elements ---")
for i in range(4):
    for j in range(4):
        if J[i,j] != 0:
            print(f"J[{i}][{j}] = {ccode(J[i,j])};")