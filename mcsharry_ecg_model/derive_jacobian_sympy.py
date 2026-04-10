import sympy as sp
from sympy.printing import ccode

T = sp.symbols('T', real=True)
tn = sp.symbols('tn', real=True)
tp = sp.symbols('tp', real=True)

xn = sp.symbols('xn0:3', real=True)
xp = sp.symbols('xp0:3', real=True)

ti = sp.symbols('ti0:5', real=True)
ai = sp.symbols('ai0:5', real=True)
bi = sp.symbols('bi0:5', real=True)

omega_n, omega_p = sp.symbols('omega_n omega_p', real=True)
fhi, A = sp.symbols('fhi A', real=True)

pi = sp.pi

def wrap_pm_pi(expr):
    return sp.atan2(sp.sin(expr), sp.cos(expr))

def ecgsyn_rhs(x, y, z, t, omega):
    alpha = 1 - sp.sqrt(x**2 + y**2)
    theta = sp.atan2(y, x)
    z0 = A * sp.sin(2 * pi * fhi * t)

    zsum = 0
    for k in range(5):
        dth = wrap_pm_pi(theta - ti[k])
        zsum += ai[k] * dth * sp.exp(-dth**2 / (2 * bi[k]**2))

    fx = alpha * x - omega * y
    fy = alpha * y + omega * x
    fz = -zsum - (z - z0)
    return [fx, fy, fz]

fp = ecgsyn_rhs(xp[0], xp[1], xp[2], tp, omega_p)
fn = ecgsyn_rhs(xn[0], xn[1], xn[2], tn, omega_n)

R = [xn[i] - xp[i] - (T/2)*(fn[i] + fp[i]) for i in range(3)]
J = sp.Matrix(R).jacobian(xn)

print("// Residuals")
for i in range(3):
    print(f"R[{i}] = {ccode(sp.simplify(R[i]))};")

print("\n// Jacobian")
for i in range(3):
    for j in range(3):
        if J[i, j] != 0:
            print(f"J[{i}][{j}] = {ccode(sp.factor(J[i, j]))};")