import numpy as np
from astropy import units as u
from poliastro.iod import izzo

MU_SUN = 1.32712440018e11

def lambert(r1, r2, tof, mu=MU_SUN):
    """
    Lambert solver using poliastro's Izzo algorithm.
    Validated implementation used in real mission planning.
    """
    r1_vec = r1 * u.km
    r2_vec = r2 * u.km
    tof_q  = tof * u.s
    mu_q   = mu * u.km**3 / u.s**2

    try:
        solutions = izzo.lambert(mu_q, r1_vec, r2_vec, tof_q, prograde=True)
        sol_list = list(solutions)
        v1 = sol_list[0].to(u.km/u.s).value
        v2 = sol_list[1].to(u.km/u.s).value
        return v1, v2
    except Exception as e:
        print(f"Lambert solver error: {e}")
        return None


if __name__ == "__main__":
    from datetime import datetime
    from bodies import get_position

    date1 = datetime(2026, 11, 27)
    date2 = datetime(2027, 6, 15)  # ~200 days later

    r1, v1_body = get_position("Earth", date1)
    r2, v2_body = get_position("Mars", date2)

    tof = (date2 - date1).total_seconds()
    result = lambert(r1, r2, tof)

    if result:
        v1_transfer, v2_transfer = result

        dv1 = np.linalg.norm(v1_transfer - v1_body)
        dv2 = np.linalg.norm(v2_transfer - v2_body)
        total_dv = dv1 + dv2

        cos_angle = np.dot(v1_transfer, v1_body) / (np.linalg.norm(v1_transfer) * np.linalg.norm(v1_body))
        angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))

        print(f"Earth departure delta-V:  {dv1:.3f} km/s")
        print(f"Mars arrival delta-V:     {dv2:.3f} km/s")
        print(f"Total delta-V:            {total_dv:.3f} km/s")
        print(f"Angle between v1 and Earth v: {angle:.2f} degrees")
        print(f"Transfer v1 magnitude: {np.linalg.norm(v1_transfer):.3f} km/s")
        print(f"Earth v magnitude:     {np.linalg.norm(v1_body):.3f} km/s")
    else:
        print("No solution found")