import streamlit as st
import pandas as pd
import plotly.express as px
import numpy as np # Make sure numpy is imported at the top of your app.py

from src.simulation import KesslerSimulation
from src.physics.collisions import calculate_impact_metrics

st.set_page_config(page_title="Kessler Syndrome Simulator", layout="wide")

@st.cache_data
def run_simulation(years: int, base_launches: int, growth_rate_pct: float, global_lifespan: int, pmd_failure_pct: float, laser_removals: int, enable_dynamic_highway_routing: bool):
    sim = KesslerSimulation()
    
    sim.base_annual_launches = base_launches
    sim.launch_growth_rate = growth_rate_pct / 100.0
    sim.sat_lifespan_years = float(global_lifespan)
    sim.pmd_failure_rate = pmd_failure_pct / 100.0   # Added this line
    sim.laser_removals_per_year = laser_removals
    sim.enable_dynamic_highway_routing = enable_dynamic_highway_routing
    
    history = []
    
# Capture initial state
    for shell in sim.shells:
            history.append({
                "Year": 0,
                "Altitude (km)": shell.altitude_km,
                "Circular Debris": shell.debris_count,
                "Heavy Debris": shell.heavy_debris_count,
                "Transiting Debris": sim.get_transiting_debris_in_shell(shell.altitude_km),
                "Derelict Satellites": shell.derelict_count,  # NEW
                "Active Satellites": shell.sat_count,
                "Diverted Launches": sim.launches_diverted_this_step,
                "Volume (km3)": shell.volume_km3
            })
        
    # Run integration
    for year in range(1, years + 1):
        sim.step(dt_years=1.0)
        for shell in sim.shells:
            history.append({
                "Year": year, # (or 0 for the first loop)
                "Altitude (km)": shell.altitude_km,
                "Circular Debris": shell.debris_count,
                "Heavy Debris": shell.heavy_debris_count,
                "Transiting Debris": sim.get_transiting_debris_in_shell(shell.altitude_km),
                "Derelict Satellites": shell.derelict_count,  # NEW
                "Active Satellites": shell.sat_count,
                "Diverted Launches": sim.launches_diverted_this_step,
                "Volume (km3)": shell.volume_km3
            })
            
    df = pd.DataFrame(history)
    # Effective debris is what satellites actually collide with
    df["Effective Debris"] = (
        df["Circular Debris"] + df["Heavy Debris"] + df["Transiting Debris"]
    )
    return df

# --- SIDEBAR CONFIGURATION ---
st.sidebar.header("Global Space Economy")
base_launches = st.sidebar.slider("Base Annual Launches", min_value=0, max_value=10000, value=2500, step=500)
launch_growth = st.sidebar.slider("Annual Launch Growth (%)", min_value=0.0, max_value=30.0, value=8.0, step=1.0)
global_lifespan = st.sidebar.slider("Constellation Lifespan (Years)", min_value=1, max_value=15, value=5)
pmd_failure_pct = st.sidebar.slider("PMD Failure Rate (%)", min_value=0.0, max_value=100.0, value=10.0, step=1.0, help="Percentage of satellites that fail to de-orbit at end of life, becoming permanent derelict debris.")

st.sidebar.header("Active Debris Removal")
num_facilities = st.sidebar.slider("Global Laser Facilities", min_value=1, max_value=10, value=1, step=1, help="Number of ground stations built worldwide.")
laser_removals = st.sidebar.slider("Total Laser Engagements / Year", min_value=0, max_value=1_000_000, value=0, step=50_000)
st.sidebar.caption("Fires on highest-risk shell. Creates 2-year transiting decay orbits.")

with st.sidebar.expander("ADR Financial Assumptions (Per Facility)"):
    laser_base_cost = st.number_input("Ground Laser CAPEX ($)", value=1_000_000_000, step=50_000_000, key="capex", help="Upfront construction cost per facility.")
    personnel_cost = st.number_input("Annual Personnel Cost ($/yr)", value=30_000_000, step=1_000_000, key="personnel", help="Wages per facility.")
    maintenance_cost = st.number_input("Annual Maintenance ($/yr)", value=20_000_000, step=1_000_000, key="maintenance", help="Hardware upkeep per facility.")
    cost_per_shot = st.number_input("Cost per Laser Blast ($)", value=2500, step=500, key="cost_shot", help="Energy and variable costs per ablation pulse.")

with st.sidebar.expander("ADR Financial Assumptions"):
    laser_base_cost = st.number_input("Ground Laser CAPEX ($)", value=1_000_000_000, step=50_000_000, help="Upfront construction cost of the facility.")
    personnel_cost = st.number_input("Annual Personnel Cost ($/yr)", value=30_000_000, step=1_000_000, help="Wages for scientists, engineers, and facility operators.")
    maintenance_cost = st.number_input("Annual Maintenance ($/yr)", value=20_000_000, step=1_000_000, help="Hardware upkeep, optics replacement, and facility maintenance.")
    cost_per_shot = st.number_input("Cost per Laser Blast ($)", value=2500, step=500, help="Energy and variable costs per ablation pulse.")

st.sidebar.header("Simulation Parameters")
sim_years = st.sidebar.slider("Years to Simulate", min_value=10, max_value=200, value=100, step=10)

st.sidebar.header("Mission Viability Profile")
enable_dynamic_highway_routing = st.sidebar.toggle(
    "Enable Dynamic Highway Routing",
    value=True,
    help="Divert launches from the 450 km Highway to the Expansion Shelf when lethal-impact survival falls below 95%.",
)
sat_altitude = st.sidebar.selectbox(
    "Viability Chart Altitude (km)", range(300, 1050, 50), index=5
)
sat_area_m2 = st.sidebar.slider("Cross-Section Area (m²)", min_value=1.0, max_value=100.0, value=10.0, step=1.0)
sat_lifespan = st.sidebar.slider("Planned Mission Duration (Years)", min_value=1, max_value=25, value=5)

# --- RUN SIMULATION ---
with st.spinner("Integrating ODEs & crossing hazards..."):
    df_history = run_simulation(
        sim_years,
        base_launches,
        launch_growth,
        global_lifespan,
        pmd_failure_pct,
        laser_removals,
        enable_dynamic_highway_routing,
    )

# --- MAIN DASHBOARD ---
st.title("🌍 Orbital Debris Cascade Simulator")
st.markdown("Modeling the Kessler Syndrome with continuous mega-constellation deployment and realistic transiting laser mitigation.")

df_start = df_history[df_history["Year"] == 0]
df_end = df_history[df_history["Year"] == sim_years]

total_debris_start = df_start["Circular Debris"].sum()
total_debris_end = df_end["Circular Debris"].sum()
net_debris = total_debris_end - total_debris_start
total_heavy_debris_end = df_end["Heavy Debris"].sum()
total_transiting_end = df_end["Transiting Debris"].sum()
total_derelicts_end = df_end["Derelict Satellites"].sum()  # NEW

total_sats_start = df_start["Active Satellites"].sum()
total_sats_end = df_end["Active Satellites"].sum()

# Adjust columns to fit the debris population metrics
col1, col2, col3, col4, col5, col6 = st.columns(6)
col1.metric("Initial Debris (>1cm)", f"{total_debris_start:,.0f}")
col2.metric(f"Circular Debris (Yr {sim_years})", f"{total_debris_end:,.0f}", delta=f"{net_debris:,.0f} frags", delta_color="inverse")
col3.metric(
    f"Heavy Debris (Yr {sim_years})",
    f"{total_heavy_debris_end:,.0f}",
    help="Tracked 10 kg structural chunks, with their own atmospheric drag profile.",
)

# The new metric tracking dead satellite buses
col4.metric(
    "Derelict Satellites", 
    f"{total_derelicts_end:,.0f}", 
    help="Massive dead satellites that failed PMD. Extremely high cross-section targets."
)

col5.metric(f"Active Sats (Yr {sim_years})", f"{total_sats_end:,.0f}")

# --- FINANCIAL CALCULATIONS ---
if laser_removals > 0:
    # Multiply fixed costs by number of global facilities
    annual_fixed_opex = (personnel_cost + maintenance_cost) * num_facilities
    annual_laser_cost = annual_fixed_opex + (laser_removals * cost_per_shot)
    total_program_cost = (laser_base_cost * num_facilities) + (annual_laser_cost * sim_years)
    
    col6.metric(
        label="ADR Program Cost", 
        value=f"${annual_laser_cost / 1_000_000:.1f}M / yr", 
        delta=f"${total_program_cost / 1_000_000_000:.2f}B Total", 
        delta_color="off"
    )
else:
    total_program_cost = 0.0
    col6.metric("ADR Program Cost", "$0 / yr", "System Offline")

diverted_launches = df_history.loc[
    df_history["Year"] == sim_years, "Diverted Launches"
].iloc[0]
st.metric(
    f"Launches Diverted Due to Congestion (Yr {sim_years})",
    f"{diverted_launches:,.0f}",
)

st.divider()

# --- TABS ---
tab_heatmap, tab_debris, tab_sats, tab_lifespan, tab_finance = st.tabs([
    "Effective Debris Heatmap", 
    "Debris Population",
    "Satellite Population", 
    "Mission Viability",
    "Program Funding"
])

with tab_heatmap:
    st.subheader("Effective Collision Hazard Over Time")
    st.caption("Combines 1 cm shrapnel, heavy structural chunks, and temporary transiting debris.")
    pivot_df = df_history.pivot(index="Altitude (km)", columns="Year", values="Effective Debris")
    fig_heat = px.imshow(pivot_df, labels=dict(x="Simulation Year", y="Altitude (km)", color="Hazard Count"), x=pivot_df.columns, y=pivot_df.index, color_continuous_scale="Inferno", aspect="auto", origin="lower")
    fig_heat.update_layout(template="plotly_dark", height=500)
    st.plotly_chart(fig_heat, use_container_width=True)

with tab_debris:
    st.subheader("Global Debris Inventory")
    st.caption("Shows 1 cm shrapnel, 10 kg heavy chunks, and temporary transiting debris. (Logarithmic Scale)")
    df_agg = df_history.groupby("Year")[[
        "Circular Debris", "Heavy Debris", "Transiting Debris"
    ]].sum().reset_index()
    
    # ADDED log_y=True here
    fig_inv = px.area(
        df_agg,
        x="Year",
        y=["Circular Debris", "Heavy Debris", "Transiting Debris"],
        color_discrete_sequence=["#EF553B", "#AB63FA", "#FFA15A"],
    )
    
    # Update layout to explicitly set the Y-axis type to log
    fig_inv.update_layout(
        template="plotly_dark", 
        height=450, 
        yaxis_title="Debris Object Count (Log Scale)", 
        legend_title="Population",
        yaxis_type="log"  # Enforces log scale on the area chart
    )
    st.plotly_chart(fig_inv, use_container_width=True)
    
with tab_sats:
    st.subheader("Active Satellite Population by Altitude")
    st.caption("Color scale is logarithmic to show detail across extreme density differences.")
    
    # 1. Pivot the raw data
    pivot_sats = df_history.pivot(index="Altitude (km)", columns="Year", values="Active Satellites")
    
    # 2. Create a log10 version of the data for the colors (add +1 to avoid log(0) errors)
    log_pivot_sats = np.log10(pivot_sats + 1)
    
    # 3. Plot using the log data for color
    fig_sats = px.imshow(
        log_pivot_sats, 
        labels=dict(x="Simulation Year", y="Altitude (km)", color="Log10(Satellites)"), 
        x=pivot_sats.columns, 
        y=pivot_sats.index, 
        color_continuous_scale="Viridis", 
        aspect="auto", 
        origin="lower"
    )
    
    # 4. Inject the raw, un-logged numbers into the hover tooltip
    fig_sats.update_traces(
        customdata=pivot_sats.values,
        hovertemplate="Year: %{x}<br>Altitude: %{y} km<br>Active Satellites: %{customdata:,.0f}<extra></extra>"
    )
    
    fig_sats.update_layout(template="plotly_dark", height=500)
    st.plotly_chart(fig_sats, use_container_width=True)

with tab_lifespan:
    st.subheader(f"Mission Viability at {sat_altitude} km")
    df_alt = df_history[df_history["Altitude (km)"] == sat_altitude].copy()
    mtbi_list = []
    prob_list = []
    for _, row in df_alt.iterrows():
        volume_km3 = row["Volume (km3)"]
        spatial_density = (
            row["Effective Debris"] + row["Derelict Satellites"]
        ) / volume_km3
        lethal_spatial_density = (
            row["Heavy Debris"] + row["Derelict Satellites"]
        ) / volume_km3
        metrics = calculate_impact_metrics(
            spatial_density,
            sat_area_m2,
            sat_lifespan,
            10.0,
            lethal_spatial_density_km3=lethal_spatial_density,
        )
        mtbi_list.append(metrics["mtbi_years"] if metrics["mtbi_years"] < 10000 else 10000)
        prob_list.append(
            metrics["probability_of_catastrophic_breach"] * 100.0
        )
    df_alt["MTBI (Years)"] = mtbi_list
    df_alt["Catastrophic Breach Risk (%)"] = prob_list
    
    st.caption(
        "Lethal-risk estimate uses heavy debris and derelict satellites only; "
        "shrapnel is excluded, although sub-threshold impacts on exposed satellite "
        "areas can cause mission kills."
    )
    fig_prob = px.bar(df_alt, x="Year", y="Catastrophic Breach Risk (%)", title=f"Probability of catastrophic breach during a {sat_lifespan}-year mission", color="Catastrophic Breach Risk (%)", color_continuous_scale="Reds")
    fig_prob.update_layout(template="plotly_dark", height=400, yaxis_range=[0, 100])
    st.plotly_chart(fig_prob, use_container_width=True)

with tab_finance:
    st.subheader("Space Sustainability Tax Calculator")
    st.markdown("Distribute the total cost of the ADR Laser program between new satellite launches and an annual residency fee for active satellites in orbit.")
    
    if laser_removals == 0:
        st.warning("Activate the Laser Broom in the sidebar to calculate funding requirements.")
    else:
        st.info(f"**Total Program Cost ({sim_years} Years):** ${total_program_cost:,.0f}")
        
        # Display the granular breakdown
        b_col1, b_col2, b_col3, b_col4 = st.columns(4)
        b_col1.metric("Total Facilities CAPEX", f"${(laser_base_cost * num_facilities) / 1_000_000_000:.2f}B")
        b_col2.metric("Total Personnel OPEX", f"${(personnel_cost * num_facilities * sim_years) / 1_000_000:.0f}M")
        b_col3.metric("Total Maintenance OPEX", f"${(maintenance_cost * num_facilities * sim_years) / 1_000_000:.0f}M")
        b_col4.metric("Total Laser Ablation Cost", f"${(cost_per_shot * laser_removals * sim_years) / 1_000_000_000:.2f}B")
        st.divider()
        
        cumulative_launches = sum(base_launches * ((1.0 + (launch_growth/100.0)) ** y) for y in range(sim_years))
        cumulative_sat_years = df_history[df_history["Year"] > 0]["Active Satellites"].sum()
        
        funding_split = st.slider(
            "Funding Split (% from Launch Tax vs % from Residency Fee)", 
            min_value=0, max_value=100, value=50, 
            help="0 = 100% Residency Fee, 100 = 100% Launch Tax"
        )
        
        launch_burden = total_program_cost * (funding_split / 100.0)
        residency_burden = total_program_cost * ((100.0 - funding_split) / 100.0)
        
        tax_per_launch = launch_burden / cumulative_launches if cumulative_launches > 0 else 0
        fee_per_sat_year = residency_burden / cumulative_sat_years if cumulative_sat_years > 0 else 0
        
        fcol1, fcol2 = st.columns(2)
        fcol1.metric("Required Tax Per New Launch", f"${tax_per_launch:,.2f}")
        fcol2.metric("Annual Residency Fee (Per Active Sat)", f"${fee_per_sat_year:,.2f} / yr")
        
        st.markdown(f"> **Simulation Context:** Over {sim_years} years, there will be **{cumulative_launches:,.0f}** new satellite deployments and **{cumulative_sat_years:,.0f}** active satellite-years tracked in orbit.")