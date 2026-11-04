import streamlit as st
import pandas as pd
import plotly.express as px

from src.simulation import KesslerSimulation
from src.physics.collisions import calculate_impact_metrics

st.set_page_config(page_title="Kessler Syndrome Simulator", layout="wide")

@st.cache_data
def run_simulation(years: int, base_launches: int, growth_rate_pct: float, global_lifespan: int, laser_removals: int):
    sim = KesslerSimulation()
    
    sim.base_annual_launches = base_launches
    sim.launch_growth_rate = growth_rate_pct / 100.0
    sim.sat_lifespan_years = float(global_lifespan)
    sim.laser_removals_per_year = laser_removals
    
    history = []
    
    for shell in sim.shells:
        history.append({
            "Year": 0,
            "Altitude (km)": shell.altitude_km,
            "Debris Count": shell.debris_count,
            "Active Satellites": shell.sat_count,
            "Volume (km3)": shell.volume_km3
        })
        
    for year in range(1, years + 1):
        sim.step(dt_years=1.0)
        for shell in sim.shells:
            history.append({
                "Year": year,
                "Altitude (km)": shell.altitude_km,
                "Debris Count": shell.debris_count,
                "Active Satellites": shell.sat_count,
                "Volume (km3)": shell.volume_km3
            })
            
    return pd.DataFrame(history)

# --- SIDEBAR CONFIGURATION ---
st.sidebar.header("Global Space Economy")
base_launches = st.sidebar.slider("Base Annual Launches", min_value=0, max_value=10000, value=2500, step=500)
launch_growth = st.sidebar.slider("Annual Launch Growth (%)", min_value=0.0, max_value=30.0, value=8.0, step=1.0)
global_lifespan = st.sidebar.slider("Constellation Lifespan (Years)", min_value=1, max_value=15, value=5)

st.sidebar.header("Active Debris Removal")
laser_removals = st.sidebar.slider("Laser Engagements / Year", min_value=0, max_value=50000, value=0, step=1000)

with st.sidebar.expander("ADR Financial Assumptions"):
    capex_cost = st.number_input("Facility Build Cost ($)", value=1_000_000_000, step=50_000_000)
    fixed_opex = st.number_input("Fixed Annual Running Cost ($)", value=50_000_000, step=5_000_000)
    cost_per_shot = st.number_input("Cost per Laser Blast ($)", value=2500, step=500)

st.sidebar.header("Simulation Parameters")
sim_years = st.sidebar.slider("Years to Simulate", min_value=10, max_value=200, value=100, step=10)

st.sidebar.header("Single Mission Profile")
sat_altitude = st.sidebar.selectbox("Target Orbit (km)", range(400, 1050, 50), index=3) 
sat_area_m2 = st.sidebar.slider("Cross-Section Area (m²)", min_value=1.0, max_value=100.0, value=10.0, step=1.0)
sat_lifespan = st.sidebar.slider("Planned Mission Duration (Years)", min_value=1, max_value=25, value=5)

# --- RUN SIMULATION ---
with st.spinner("Calculating orbital cascade..."):
    df_history = run_simulation(sim_years, base_launches, launch_growth, global_lifespan, laser_removals)

# --- MAIN DASHBOARD ---
st.title("🌍 Orbital Debris Cascade Simulator")
st.markdown("Modeling the Kessler Syndrome with continuous mega-constellation deployment and active laser remediation.")

df_start = df_history[df_history["Year"] == 0]
df_end = df_history[df_history["Year"] == sim_years]

total_debris_start = df_start["Debris Count"].sum()
total_debris_end = df_end["Debris Count"].sum()
net_debris = total_debris_end - total_debris_start

total_sats_start = df_start["Active Satellites"].sum()
total_sats_end = df_end["Active Satellites"].sum()

col1, col2, col3, col4, col5 = st.columns(5)
col1.metric("Initial Debris (>1cm)", f"{total_debris_start:,.0f}")
col2.metric(f"Debris (Year {sim_years})", f"{total_debris_end:,.0f}", delta=f"{net_debris:,.0f} frags", delta_color="inverse")

# The new metric explicitly tracking the CelesTrak starting population
col3.metric("Initial Active Sats", f"{total_sats_start:,.0f}")

col4.metric(f"Active Sats (Year {sim_years})", f"{total_sats_end:,.0f}", delta=f"{total_sats_end - total_sats_start:,.0f} in orbit")

# --- FINANCIAL CALCULATIONS ---
if laser_removals > 0:
    annual_laser_cost = fixed_opex + (laser_removals * cost_per_shot)
    total_program_cost = capex_cost + (annual_laser_cost * sim_years)
    
    col5.metric(
        label="ADR Program Cost (Avg/yr)", 
        value=f"${annual_laser_cost / 1_000_000:.1f}M / yr", 
        delta=f"${total_program_cost / 1_000_000_000:.2f}B Total (inc Build)", 
        delta_color="off"
    )
else:
    col5.metric("ADR Program Cost", "$0 / yr", "System Offline")
    
st.divider()

# --- TABS ---
tab_heatmap, tab_sats, tab_lifespan = st.tabs(["Debris Cascade Heatmap", "Satellite Population", "Mission Viability"])

with tab_heatmap:
    st.subheader("Orbital Debris Contamination Over Time")
    pivot_df = df_history.pivot(index="Altitude (km)", columns="Year", values="Debris Count")
    
    fig_heat = px.imshow(
        pivot_df,
        labels=dict(x="Simulation Year", y="Altitude (km)", color="Debris Count"),
        x=pivot_df.columns,
        y=pivot_df.index,
        color_continuous_scale="Inferno",
        aspect="auto",
        origin="lower"
    )
    fig_heat.update_layout(template="plotly_dark", height=500)
    st.plotly_chart(fig_heat, use_container_width=True)
    
with tab_sats:
    st.subheader("Active Satellite Population by Altitude")
    st.markdown("Observe the equilibrium between continuous launches and planned de-orbits, and watch for population crashes caused by debris cascades.")
    pivot_sats = df_history.pivot(index="Altitude (km)", columns="Year", values="Active Satellites")
    
    fig_sats = px.imshow(
        pivot_sats,
        labels=dict(x="Simulation Year", y="Altitude (km)", color="Active Satellites"),
        x=pivot_sats.columns,
        y=pivot_sats.index,
        color_continuous_scale="Viridis",
        aspect="auto",
        origin="lower"
    )
    fig_sats.update_layout(template="plotly_dark", height=500)
    st.plotly_chart(fig_sats, use_container_width=True)

with tab_lifespan:
    st.subheader(f"Mission Viability at {sat_altitude} km")
    df_alt = df_history[df_history["Altitude (km)"] == sat_altitude].copy()
    
    mtbi_list = []
    prob_list = []
    
    for _, row in df_alt.iterrows():
        spatial_density = row["Debris Count"] / row["Volume (km3)"]
        metrics = calculate_impact_metrics(
            spatial_density_km3=spatial_density,
            cross_section_m2=sat_area_m2,
            time_horizon_years=sat_lifespan
        )
        
        mtbi = metrics["mtbi_years"] if metrics["mtbi_years"] < 10000 else 10000
        mtbi_list.append(mtbi)
        prob_list.append(metrics["probability_of_impact"] * 100.0)
        
    df_alt["MTBI (Years)"] = mtbi_list
    df_alt["Catastrophic Breach Risk (%)"] = prob_list
    
    fig_prob = px.bar(
        df_alt, x="Year", y="Catastrophic Breach Risk (%)",
        title=f"Probability of being hit during a {sat_lifespan}-year mission",
        color="Catastrophic Breach Risk (%)",
        color_continuous_scale="Reds"
    )
    fig_prob.update_layout(template="plotly_dark", height=400, yaxis_range=[0, 100])
    st.plotly_chart(fig_prob, use_container_width=True)