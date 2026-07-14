# GNC Mission Design Dashboard

Interactive web tool for Earth-Mars trajectory analysis. Select a departure window, compute optimal launch dates, and visualize transfer trajectories 
with delta-V budgets.

## What it does

- Generates porkchop plots showing total delta-V across all departure and arrival date combinations
- Identifies optimal launch windows with minimum fuel cost
- Displays time of flight, departure delta-V, and arrival delta-V for any selected trajectory
- Renders the transfer trajectory in 3D with Earth, Mars, and Sun positions

## Core math

- Lambert's problem (universal variable method with Stumpff functions) — given two position vectors and a time of flight, computes the transfer orbit velocities
- Two-body orbital mechanics for planetary ephemeris (Kepler's equation, true anomaly propagation)
- Delta-V computed as the difference between transfer orbit velocity and planetary velocity at departure and arrival

## Why it matters

Delta-V is the fundamental currency of mission design — it determines propellant mass and therefore mission cost and feasibility. 
This tool replicates the core workflow mission designers use to identify viable launch windows before committing to a trajectory.

## Stack

- Python, NumPy
- poliastro — validated Lambert solver
- Plotly / Dash — interactive porkchop heatmap and 3D trajectory visualization

## Background

Built as part of a summer portfolio project targeting GNC and mission design roles in the space and eVTOL industry. 
Aerospace Engineering senior at CU Boulder.
