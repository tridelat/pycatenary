from pycatenary import cable
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt
import numpy as np

# Cable length [m]
L = 10
R1x = 0.2 * 5 # Horizontal force at fairlead
R1y = 7 # Vertical force at fairlead
p = 1 # Lineic weight

# Indices
i_x = 0
i_y = 1
i_T = 2
i_theta = 3


# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_bvp.html
def fun(s, state):
    '''Returns the derivative of state vector which describes the ordinary differential equation of catenary'''
    param = {
        "f": 1,
        "p": 1
        }

    return np.vstack((np.cos(state[i_theta]), \
                      np.sin(state[i_theta]), \
                      param["p"] * np.sin(state[i_theta]),\
           (param["f"] + param["p"] * np.cos(state[i_theta])) / state[i_T]))

def bc(ya, yb):
    '''Returns a vector of functions have to zero to satisfy boundary conditions '''
    return np.array(
      [ya[i_x] - 0,\ # x coordinate of first extremity of cable
       ya[i_y] - 0,\ # y coordinate of first extremity of cable
       yb[i_T] * np.cos(yb[i_theta]) - R1x,\ # Horizontal equilibrium of force at other extremity
       yb[i_T] * np.sin(yb[i_theta]) - R1y]) # Vertical equilibrium of force at other extremity

def analytical_solution(s, p, R1x, R1y, L):
    """
    :param p: linear weight [N/m]
    :param R1x: horizontal force at end of catenary [N]
    :param R1y: vertical force at end of catenary [N]
    :param L: length of catenary [m]
    :return: horizontal and vertical coordinates describing catenary curve with weight
    only according to https://api-mecaspa.pagesperso-orange.fr/aero_energies/eolien_volant/idees_sur_le_calcul_du_cable.htm

    """
    K = p / R1x # Ratio of lineic weight to constant horizontal force in catenary
    ttheta0 = (R1y - p * L) / R1x  # Tangent of cable inclination at starting extremity
    lambda0 = -1 / K * np.arcsinh(ttheta0)
    x = 1 / K * np.arcsinh(s * p / R1x - np.sinh(K * lambda0)) + lambda0
    y = 1 / K * (np.cosh(np.arcsinh(s * p / R1x - np.sinh(K * lambda0))) - np.cosh(K * lambda0))
    return (x, y)

def init_guess(s, p, R1x, R1y, L):
    (x, y) = analytical_solution(s, p, R1x, R1y, L)

    ttheta0 = (R1y - p * L) / R1x
    guess0 = np.vstack((x, y, R1y - p * (L - s), np.arctan(ttheta0) * (L - s) / L + np.arctan(R1y / R1x) * s / L))
    return guess0

# curvilinear abscissa
curvi = np.linspace(0, L, 100)
(x, y) = analytical_solution(curvi, p, R1x, R1y, L)
guess0 = init_guess(curvi, p, R1x, R1y, L)
print(guess0)

print(x, y)
plt.plot(x, y, label='Geom')
plt.show()
sol_init = np.zeros((4, curvi.size))
sol_init[i_x, :] = x
sol_init[i_y, :] = y
sol_init[i_T, :] = guess0[i_T, :]
sol_init[i_theta, :] = guess0[i_theta, :]

param = {
  "f": 10,
  "p": 10
}
p = np.array([0, 0])
p[0] = param["f"]
p[1] = param["p"]
res_b = solve_bvp(fun, bc, curvi, sol_init)
x_plot = np.linspace(0, 1, L)
y_plot_b = res_b.sol(x_plot)[0]
print(res_b)
plt.plot(res_b.y[i_x, :], res_b.y[i_y, :], label='Geom')

plt.legend()
plt.xlabel("x")
plt.ylabel("y")

plt.show()

# define properties of cable
length = L    # length of line
w = p    # submerged weight
EA = 560e3    # axial stiffness (for elastic catenary)
floor = False    # if True, contact is possible at the level of the anchor
anchor = [0., 0., 0.]
fairlead = [x[-1], 0., y[-1]]

# create cable instance
l1 = cable.MooringLine(L=length,
                       w = w,
                       EA = None,
                       anchor = anchor,
                       fairlead = fairlead,
                       floor = floor)

# compute calculation
l1.computeSolution()

# change fairlead position
l1.setFairleadCoords([x[-1], 0., y[-1]])

# recompute solution
l1.computeSolution()

# get tension along line (between 0. and total line length)
T = curvi.copy()
for index, item in enumerate(curvi):
    print(item)
    print(index)
    print(l1.getTension(item))
    T[index] = np.linalg.norm(l1.getTension(item))
print(T)

# get xyz coordinates along line
print(L)
print(curvi[-1])
xyz = 0 * guess0[0:3, :]
for index, item in enumerate(curvi[0:-2]):
     xyz[:, index] = l1.s2xyz(item)
plt.plot(xyz[0, :], xyz[2, :], label='Geom')
plt.show()

