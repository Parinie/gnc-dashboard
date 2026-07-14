import numpy as np
from datetime import datetime, timedelta
from lambert import lambert
from bodies import get_position

def compute_porkchop(departure_start, departure_end, arrival_start, arrival_end,
                     departure_steps=30, arrival_steps=30):
    """
    Compute delta-V grid for all departure/arrival date combinations.
    Returns departure dates, arrival dates, and delta-V matrix.
    """
    departure_dates = [departure_start + timedelta(days=i * (departure_end - departure_start).days / departure_steps)
                       for i in range(arrival_steps + 1)]
    
    arrival_dates = [arrival_start + timedelta(days=i * (arrival_end - arrival_start).days / arrival_steps)
                     for i in range(arrival_steps + 1)]
    
    dv_matrix = np.full((len(arrival_dates), len(departure_dates)), np.nan)

    print(f"Computing {len(departure_dates) * len(arrival_dates)} trajectories...")
    
    for i, d_date in enumerate(departure_dates):
        for j, a_date in enumerate(arrival_dates):
            # Arrival must be after departure
            if a_date <= d_date:
                continue
            
            tof = (a_date - d_date).total_seconds()

            # Skip very short or very long transfers
            tof_days = tof / 86400
            if tof_days < 60 or tof_days > 500:
                continue
            
            r1, v1_body = get_position('Earth', d_date)
            r2, v2_body = get_position('Mars', a_date)

            result = lambert(r1, r2, tof)

            if result:
                v1_transfer, v2_transfer = result
                dv1 = np.linalg.norm(v1_transfer - v1_body)
                dv2 = np.linalg.norm(v2_transfer - v2_body)
                total_dv = dv1 + dv2

                # Cap unrealistic values
                if total_dv < 20:
                    dv_matrix[j, i] = total_dv
    
    return departure_dates, arrival_dates, dv_matrix

if __name__ == "__main__":
    import plotly.graph_objects as go

    # Scan around the good 2026-2027 window
    dep_start = datetime(2026, 9, 1)
    dep_end   = datetime(2027, 2, 1)
    arr_start = datetime(2027, 1, 1)
    arr_end   = datetime(2027, 10, 1)

    dep_dates, arr_dates, dv_matrix = compute_porkchop(
        dep_start, dep_end, arr_start, arr_end,
        departure_steps=40, arrival_steps=40
    )

    dep_strings = [d.strftime("%Y-%m-%d") for d in dep_dates]
    arr_strings = [d.strftime("%Y-%m-%d") for d in arr_dates]

fig = go.Figure(data=go.Heatmap(
        z=dv_matrix,
        x=dep_strings,
        y=arr_strings,
        colorscale="Viridis_r",
        colorbar=dict(title="Total ΔV (km/s)"),
        zmin=3,
        zmax=15
    ))

fig.update_layout(
    title="Earth → Mars porkchop plot (2026–2027)",
    xaxis_title="Departure date",
    yaxis_title="Arrival date",
    width=900,
    height=700
)

fig.show()
print("\nDone! Look for the dark region — that's your optimal launch window.")

