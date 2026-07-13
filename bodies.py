import numpy as np
from datetime import datetime

# Gravitational parameter of the Sum (km^3/s^2)
MU_SUN = 1.32712440018e11

# Orbital elements for earth and mars (J2000 epoch)
# [semi_major_axis (km), eccentricity, inclination (rad),
# RAAN (rad), arg_perihelion (rad), mean_anomaly_at_epoch (rad)]

BODIES = {
    "Earth": {
        "a": 1.000 * 1.496e8,    #km
        "e": 0.0167,
        "i": 0.0,
        "T": 365.25,              # orbital period in days
        "L0": 100.464            # mean longitude at J2000 (degrees)
    },
    "Mars": {
        "a": 1.524 * 1.496e8,    # km
        "e": 0.0934,
        "i": 0.0322,             # radians
        "T": 686.97,              # orbital period in days
        "L0": 355.453
    }
}

J2000 = datetime(2000, 1, 1, 12, 0, 0)

def days_since_j2000(date):
    return (date - J2000).total_seconds() / 86400

def get_position(body_name, date):
    """
    Returns position vector (km) and velocity vector (km/s)
    for a solar system body at a given date.
    Uses simplified ciruclar orbit approximation.
    """

    body = BODIES[body_name]
    a = body["a"]
    e = body["e"]
    T = body["T"]
    L0 = np.radians(body["L0"])

    days = days_since_j2000(date)

    # Mean motion (rad/day)
    n = 2 * np.pi / (T * 86400) # rad/s

    # Mean anomaly with correct starting longitude
    M = (L0 + n * days * 86400) % (2 * np.pi)

    # Solve Kepler's equation iteratively for eccentric anomaly 
    E = M
    for i in range(100):
        E_new = M + e * np.sin(E)
        if abs(E_new - E) < 1e-12:
            break
        E = E_new
    
    # True anomaly 
    nu = 2 * np.arctan2(
        np.sqrt(1 + e) * np.sin(E / 2),
        np.sqrt(1 - e) * np.cos(E / 2)
    )

    # Distance from sun
    r = a * (1 - e * np.cos(E))

    # Position in orbital plane 
    x = r * np.cos(nu)
    y = r * np.sin(nu)

    # Velocity components in orbital plane
    h = np.sqrt(MU_SUN * a * (1 - e**2))  # specific angular momentum
    vx = -(MU_SUN / h) * np.sin(nu)
    vy = (MU_SUN / h) * (e + np.cos(nu))

    return np.array([x, y, 0.0]), np.array([vx, vy, 0.0])

# Quick sanity check
if __name__ == "__main__":
    date = datetime(2026, 1, 1)
    r_earth, v_earth = get_position("Earth", date)
    r_mars, v_mars = get_position("Mars", date)

    print(f"Earth position: {np.linalg.norm(r_earth)/1.496e8:.3f} AU from Sun")
    print(f"Mars position:  {np.linalg.norm(r_mars)/1.496e8:.3f} AU from Sun")