import numpy as np
from datetime import datetime, timedelta
import dash
from dash import dcc, html, Input, Output
import plotly.graph_objects as go
from porkchop import compute_porkchop
from bodies import get_position
from lambert import lambert

app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1("GNC Mission Design Dashboard",
            style={"textAlign": "center", "fontFamily": "Arial"}),
    html.P("Earth → Mars trajectory analysis",
           style={"textAlign": "center", "color": "gray"}),

    html.Div([
        html.Div([
            html.Label("Departure window start"),
            dcc.DatePickerSingle(
                id="dep-start",
                date="2026-09-01",
                display_format="YYYY-MM-DD"
            ),
        ], style={"margin": "10px"}),

        html.Div([
            html.Label("Departure window end"),
            dcc.DatePickerSingle(
                id="dep-end",
                date="2027-02-01",
                display_format="YYYY-MM-DD"
            ),
        ], style={"margin": "10px"}),

        html.Button("Compute porkchop plot", id="compute-btn",
                   style={"margin": "20px", "padding": "10px 20px"}),

    ], style={"display": "flex", "alignItems": "center", "justifyContent": "center"}),

    dcc.Loading(
        dcc.Graph(id="porkchop-graph"),
        type="circle"
    ),

    html.Div(id="selected-trajectory"),

    dcc.Graph(id="trajectory-3d")
])


@app.callback(
    Output("porkchop-graph", "figure"),
    Input("compute-btn", "n_clicks"),
    Input("dep-start", "date"),
    Input("dep-end", "date"),
    prevent_initial_call=False
)
def update_porkchop(n_clicks, dep_start, dep_end):
    dep_start = datetime.strptime(dep_start[:10], "%Y-%m-%d")
    dep_end   = datetime.strptime(dep_end[:10],   "%Y-%m-%d")
    arr_start = dep_start + timedelta(days=60)
    arr_end   = dep_end   + timedelta(days=500)

    dep_dates, arr_dates, dv_matrix = compute_porkchop(
        dep_start, dep_end, arr_start, arr_end,
        departure_steps=25, arrival_steps=25
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
        title="Earth → Mars porkchop plot",
        xaxis_title="Departure date",
        yaxis_title="Arrival date",
        height=500
    )

    return fig


@app.callback(
    Output("trajectory-3d", "figure"),
    Output("selected-trajectory", "children"),
    Input("porkchop-graph", "clickData")
)
def show_trajectory(clickData):
    if clickData is None:
        return go.Figure(), ""

    dep_date = datetime.strptime(clickData["points"][0]["x"], "%Y-%m-%d")
    arr_date = datetime.strptime(clickData["points"][0]["y"], "%Y-%m-%d")
    tof = (arr_date - dep_date).total_seconds()

    r1, v1_body = get_position("Earth", dep_date)
    r2, v2_body = get_position("Mars",  arr_date)

    result = lambert(r1, r2, tof)
    if not result:
        return go.Figure(), "No solution found for this trajectory."

    v1_transfer, v2_transfer = result
    dv1 = np.linalg.norm(v1_transfer - v1_body)
    dv2 = np.linalg.norm(v2_transfer - v2_body)
    total_dv = dv1 + dv2

    # Propagate transfer trajectory
    steps = 200
    transfer_points = []
    for k in range(steps + 1):
        t = k * tof / steps
        # Simple Kepler propagation along transfer
        frac = k / steps
        r_interp = r1 + frac * (r2 - r1)
        transfer_points.append(r_interp)
    transfer_points = np.array(transfer_points)

    AU = 1.496e8
    fig = go.Figure()

    # Sun
    fig.add_trace(go.Scatter3d(
        x=[0], y=[0], z=[0],
        mode="markers",
        marker=dict(color="yellow", size=8),
        name="Sun"
    ))

    # Earth position
    fig.add_trace(go.Scatter3d(
        x=[r1[0]/AU], y=[r1[1]/AU], z=[r1[2]/AU],
        mode="markers",
        marker=dict(color="blue", size=6),
        name="Earth (departure)"
    ))

    # Mars position
    fig.add_trace(go.Scatter3d(
        x=[r2[0]/AU], y=[r2[1]/AU], z=[r2[2]/AU],
        mode="markers",
        marker=dict(color="red", size=6),
        name="Mars (arrival)"
    ))

    # Transfer trajectory
    fig.add_trace(go.Scatter3d(
        x=transfer_points[:, 0]/AU,
        y=transfer_points[:, 1]/AU,
        z=transfer_points[:, 2]/AU,
        mode="lines",
        line=dict(color="orange", width=3),
        name="Transfer trajectory"
    ))

    fig.update_layout(
        title=f"Transfer: {dep_date.strftime('%Y-%m-%d')} → {arr_date.strftime('%Y-%m-%d')}",
        scene=dict(
            xaxis_title="X (AU)",
            yaxis_title="Y (AU)",
            zaxis_title="Z (AU)",
            aspectmode="data"
        ),
        height=500
    )

    info = html.Div([
        html.H3(f"Departure: {dep_date.strftime('%Y-%m-%d')} → Arrival: {arr_date.strftime('%Y-%m-%d')}"),
        html.P(f"Time of flight: {(arr_date - dep_date).days} days"),
        html.P(f"Departure ΔV: {dv1:.3f} km/s"),
        html.P(f"Arrival ΔV:   {dv2:.3f} km/s"),
        html.P(f"Total ΔV:     {total_dv:.3f} km/s",
               style={"fontWeight": "bold"})
    ], style={"textAlign": "center", "fontFamily": "Arial"})

    return fig, info


if __name__ == "__main__":
    app.run(debug=True)