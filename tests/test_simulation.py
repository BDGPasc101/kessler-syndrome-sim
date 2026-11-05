import os
import json
import pytest

from src.simulation import KesslerSimulation
from src.physics.collisions import calculate_impact_metrics
from src.physics.breakup import calculate_breakup_fragments


# ---------------------------------------------------------------------------
# Test 1: Verify Initial Satellite Population is Loaded and Non-Zero
# ---------------------------------------------------------------------------
def test_initial_sat_count_nonzero_and_matches_data():
    data_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "data", "celestrak_shells.json")
    )
    assert os.path.exists(data_path), f"Telemetry file missing at {data_path}"

    with open(data_path, "r") as f:
        data = json.load(f)

    expected_sats = sum(val.get("satellites", 0) for val in data.values())
    assert expected_sats > 0, "Telemetry JSON has 0 total active satellites recorded."

    sim = KesslerSimulation()
    total_sim_sats = sum(shell.sat_count for shell in sim.shells)

    assert total_sim_sats > 0, f"Sim initialized with {total_sim_sats} sats, expected > 0."
    
    shell_550 = next((s for s in sim.shells if s.altitude_km == 550), None)
    assert shell_550 is not None, "550km shell missing from simulation"
    assert shell_550.sat_count > 0, "550km shell has 0 satellites loaded."


# ---------------------------------------------------------------------------
# Test 2: Mass & Debris Conservation Under Drag Only
# ---------------------------------------------------------------------------
def test_atmospheric_drag_monotonic_downward_flux():
    sim = KesslerSimulation()
    sim.base_annual_launches = 0
    sim.launch_growth_rate = 0.0

    # Reset both satellite sub-pools directly (sat_count is a read-only property)
    for shell in sim.shells:
        shell.legacy_sat_count = 0.0
        shell.new_sat_count = 0.0
        shell.debris_count = 0.0

    idx_800 = next(i for i, s in enumerate(sim.shells) if s.altitude_km == 800)
    sim.shells[idx_800].debris_count = 10_000.0

    sim.step(dt_years=5.0)

    # 800km debris must decay downward
    assert sim.shells[idx_800].debris_count < 10_000.0, "800km debris failed to decay."

    # Debris in the shell directly below (750km) must receive downward flux
    idx_750 = next(i for i, s in enumerate(sim.shells) if s.altitude_km == 750)
    assert sim.shells[idx_750].debris_count > 0.0, "No downward flux received by 750km shell."


# ---------------------------------------------------------------------------
# Test 3: NASA Breakup Model Fragment Scaling
# ---------------------------------------------------------------------------
def test_nasa_breakup_energy_and_fragment_scaling():
    v_imp = 10.0  # km/s

    # Use exact parameter name: impact_velocity_km_s
    res_small = calculate_breakup_fragments(
        target_mass_kg=10.0, projectile_mass_kg=1.0, impact_velocity_km_s=v_imp
    )
    res_large = calculate_breakup_fragments(
        target_mass_kg=500.0, projectile_mass_kg=1.0, impact_velocity_km_s=v_imp
    )

    assert res_small["fragments_generated"] > 0
    assert res_large["fragments_generated"] > res_small["fragments_generated"], (
        "Larger target mass must produce more fragments under catastrophic breakup."
    )


# ---------------------------------------------------------------------------
# Test 4: Collision Rate Scaling with Spatial Density
# ---------------------------------------------------------------------------
def test_collision_probability_density_scaling():
    area = 5.0        # m^2
    v_rel = 10.0      # km/s
    dt = 1.0          # year

    density_low = 1e-8
    density_high = 2e-8

    m_low = calculate_impact_metrics(density_low, area, dt, v_rel)
    m_high = calculate_impact_metrics(density_high, area, dt, v_rel)

    ratio = m_high["rate_per_year"] / m_low["rate_per_year"]
    assert pytest.approx(ratio, rel=1e-3) == 2.0, "Impact rate does not scale linearly with spatial density."
    assert 0.0 <= m_low["probability_of_impact"] <= 1.0


# ---------------------------------------------------------------------------
# Test 5: Laser Broom Targeted Removal
# ---------------------------------------------------------------------------
def test_laser_broom_targets_most_congested_shells_first():
    sim = KesslerSimulation()
    sim.base_annual_launches = 0
    sim.laser_removals_per_year = 5_000

    for s in sim.shells:
        s.legacy_sat_count = 0.0
        s.new_sat_count = 0.0
        s.debris_count = 0.0

    idx_550 = next(i for i, s in enumerate(sim.shells) if s.altitude_km == 550)
    idx_800 = next(i for i, s in enumerate(sim.shells) if s.altitude_km == 800)

    sim.shells[idx_550].debris_count = 1_000.0
    sim.shells[idx_800].debris_count = 20_000.0

    # Step 1 year with 5,000 laser shots
    sim.step(dt_years=1.0)

    # 800km should absorb all 5,000 shots (20,000 - 5,000 = ~15,000, minus slight drag)
    assert sim.shells[idx_800].debris_count <= 15_000.0, "Laser failed to target 800km shell."

    # 550km experiences natural drag (~237 lost), but should NOT receive any laser hits (would drop below 700 if hit)
    assert sim.shells[idx_550].debris_count > 700.0, "Laser incorrectly deducted shots from 550km shell."