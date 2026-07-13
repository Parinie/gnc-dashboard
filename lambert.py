import numpy as np

MU_SUN = 1.32712440018e11  # Gravitational parameter of the Sun (km^3/s^2

def lambert(r1, r2, tof, mu=MU_SUN):
    """
    Solve Lambert's problem: given two position vectors and time of flight,
    find the velocities at departure and arrival.

    r1 : departure position vector (km)
    r2 : arrival position vector (km)
    tof : time of flight (sconds)

    Returns v1 (km/s) and v2 (km/s), or None if no solution found.
    """

    r1_norm = np.linalg.norm(r1)
    r2_norm = np.linalg.norm(r2)

    # Cosine of transfer angle
    cos_dnu = np.dot(r1, r2) / (r1_norm * r2_norm)
    cos_dnu = np.clip(cos_dnu, -1, 1)

    #Assume prograde transfer 
    cross = np.cross(r1, r2)
    if cross[2] >= 0:
        dnu = np.arccos(cos_dnu)
    else:
        dnu = 2 * np.pi - np.arccos(cos_dnu)

    # Stumpuff functions for universal variable formulation
    def stumpff_c(psi):
        if psi > 1e-6:
            return (1 - np.cos(np.sqrt(psi))) / psi
        elif psi < -1e-6:
            return (np.cosh(np.sqrt(-psi)) - 1) / (-psi)
        return 0.5
    
    def stumpff_s(psi):
        if psi > 1e-6:
            sp = np.sqrt(psi)
            return (sp - np.sin(sp)) / (psi * sp)
        elif psi < -1e-6:
            sp = np.sqrt(-psi)
            return (np.sinh(sp) - sp) / ((-psi) * sp)
        return 1/6
    
    # A paramter
    A = np.sin(dnu) * np.sqrt(r1_norm * r2_norm / (1 - cos_dnu))

    if A == 0:
        return None
    
    # Iterate to find psi
    psi = 0.0
    psi_low = -4 * np.pi**2
    psi_high = 4 * np.pi**2

    for i in range(1000):
        c2 = stumpff_c(psi)
        c3 = stumpff_s(psi)

        if c2 == 0:
            return None
        
        y = r1_norm + r2_norm + A * (psi * c3 -1) / np.sqrt(c2)

        if A > 0 and y < 0:
            psi_low = psi
            psi += 0.1
            continue

        if y < 0:
            return None

        chi = np.sqrt(y / c2)
        tof_calc = (chi**3 * c3 + A * np.sqrt(y)) / np.sqrt(mu)

        if abs(tof_calc - tof) < 1e-6:
            break
        elif tof_calc < tof:
            psi_low = psi
        else:
            psi_high = psi
        
        psi = (psi_low + psi_high) / 2

    # Compute Lagrange coeeficients
    f = 1 - y / r1_norm
    g = A * np.sqrt(y / mu)
    g_dot = 1 - y / r2_norm

    v1 = (r2 - f * r1) / g
    v2 = (g_dot * r2 - r1) / g

    return v1, v2

if __name__ == "__main__":
    from datetime import datetime
    from bodies import get_position

    # Test: Earth to Mars, 200 day transfer
    date1 = datetime(2026, 11, 27)
    date2 = datetime(2027, 6, 15)  # ~200 days later

    r1, v1_body = get_position("Earth", date1)
    r2, v2_body = get_position("Mars", date2)

    tof = (date2 - date1).total_seconds()

    result = lambert(r1, r2, tof)

    if result:
        v1_transfer, v2_transfer = result

        # Delta-V at departure (difference between transfer orbit and Earth's velocity)
        dv1 = np.linalg.norm(v1_transfer - v1_body)
        dv2 = np.linalg.norm(v2_transfer - v2_body)
        total_dv = dv1 + dv2

        print(f"Earth departure delta-V:  {dv1:.3f} km/s")
        print(f"Mars arrival delta-V:     {dv2:.3f} km/s")
        print(f"Total delta-V:            {total_dv:.3f} km/s")
    else:
        print("No solution found")

"""
Understanding the variables:

psi
It's a universal variable that captures the geometry of the transfer orbit. 
It's related to the change in eccentric anomaly between departure and arrival. 
Positive psi means an elliptical transfer, negative means hyperbolic, zero 
means parabolic. 
We iterate on psi until the computed time of flight matches the one we want.

Stumpff functions (c2 and c3)
These are special mathematical functions that handle all three orbit types 
(elliptical, parabolic, hyperbolic) in one unified equation.
Without them you'd need three separate sets of equations depending on orbit type.

A (the geometry parameter)
A captures the geometry of the transfer — specifically the relationship between 
the two position vectors and the transfer angle. It's derived from the cross 
product of r1 and r2. 

Chi
universal variable itself — it's a generalized anomaly that works for any 
orbit type. 
Once we've found psi, chi is computed from it and used to find the Lagrange 
coefficients f, g, and g_dot.

Lagrange coefficients (f, g, g_dot)
They let you express the final position and velocity purely in terms of the 
initial position and velocity without needing to know the orbit type explicitly.
r2 = f * r1 + g * v1
v2 = f_dot * r1 + g_dot * v1

Overall flow of the algorithm:
1. Compute transfer angle between r1 and r2
2. Guess psi (start at 0)
3. Compute Stumpff functions from psi
4. Compute what time of flight that psi implies
5. Compare to desired time of flight
6. Adjust psi up or down (bisection search)
7. Repeat until converged
8. Use converged psi to compute Lagrange coefficients
9. Solve for v1 and v2
"""

