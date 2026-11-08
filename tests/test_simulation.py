import os
import json
import pytest

from src.simulation import KesslerSimulation, TransitingDebrisBatch
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

    for shell in sim.shells:
        shell.legacy_sat_count = 0.0
        shell.new_sat_count = 0.0
        shell.debris_count = 0.0

    idx_800 = next(i for i, s in enumerate(sim.shells) if s.altitude_km == 800)
    sim.shells[idx_800].debris_count = 10_000.0

    sim.step(dt_years=5.0)

    assert sim.shells[idx_800].debris_count < 10_000.0, "800km debris failed to decay."
    idx_750 = next(i for i, s in enumerate(sim.shells) if s.altitude_km == 750)
    assert sim.shells[idx_750].debris_count > 0.0, "No downward flux received by 750km shell."
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-7)


# ---------------------------------------------------------------------------
# Test 3: NASA Breakup Model Fragment Scaling
# ---------------------------------------------------------------------------
def test_nasa_breakup_energy_and_fragment_scaling():
    v_imp = 10.0  # km/s

    res_small = calculate_breakup_fragments(
        target_mass_kg=10.0, projectile_mass_kg=1.0, impact_velocity_km_s=v_imp
    )
    res_large = calculate_breakup_fragments(
        target_mass_kg=500.0, projectile_mass_kg=1.0, impact_velocity_km_s=v_imp
    )

    assert res_small["fragments_generated"] > 0
    assert res_large["fragments_generated"] > res_small["fragments_generated"]


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
    assert pytest.approx(ratio, rel=1e-3) == 2.0
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

    sim.step(dt_years=1.0)

    assert sim.shells[idx_800].debris_count <= 15_000.0, "Laser failed to target 800km shell."
    assert sim.shells[idx_550].debris_count > 700.0, "Laser incorrectly deducted shots from 550km shell."


# ---------------------------------------------------------------------------
# Test 6: Laser Transit Debris Generates Temporary Lower-Shell Flux
# ---------------------------------------------------------------------------
def test_laser_debris_transit_flux_and_burnup():
    sim = KesslerSimulation()
    sim.base_annual_launches = 0
    sim.launch_growth_rate = 0.0
    sim.laser_removals_per_year = 2_000.0
    sim.laser_transit_duration_years = 2.0

    for s in sim.shells:
        s.legacy_sat_count = 0.0
        s.new_sat_count = 0.0
        s.debris_count = 0.0

    idx_800 = next(i for i, s in enumerate(sim.shells) if s.altitude_km == 800)
    sim.shells[idx_800].debris_count = 10_000.0

    assert sim.get_transiting_debris_in_shell(550.0) == 0.0

    # Year 1: Laser fires 2,000 shots at 800km
    sim.step(dt_years=1.0)

    assert sim.shells[idx_800].debris_count < 10_000.0
    assert len(sim.transiting_batches) > 0
    assert sim.total_active_transiting_debris > 0.0
    assert sim.get_transiting_debris_in_shell(550.0) > 0.0

    # Step Year 2 and Year 3 with laser disabled to isolate burnup of Year 1 cohort
    sim.laser_removals_per_year = 0.0
    sim.step(dt_years=1.0)
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-7)
    sim.step(dt_years=1.0)
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-7)

    # After full transit duration, the batch has completed atmospheric reentry
    assert sim.get_transiting_debris_in_shell(550.0) == 0.0
    assert len(sim.transiting_batches) == 0

from src.physics.breakup import calculate_breakup_fragments

def test_breakup_mass_conservation_limit():
    # Setup a high-energy impact that mathematically generates ~323 raw fragments
    target_mass = 250.0
    projectile_mass = 1.0
    fragment_mass = 1.0
    
    result = calculate_breakup_fragments(
        target_mass_kg=target_mass,
        projectile_mass_kg=projectile_mass,
        impact_velocity_km_s=10.0,
        fragment_mass_kg=fragment_mass
    )
    
    total_mass = target_mass + projectile_mass
    max_allowed_fragments = total_mass / fragment_mass
    
    # Verify the output was successfully clamped to the physical limit (251)
    assert result["fragments_generated"] <= max_allowed_fragments, \
        f"Mass leak: {result['fragments_generated']} fragments exceed {max_allowed_fragments} max allowed."
    
    # Explicitly check the exact clamp value
    assert result["fragments_generated"] == 251
    assert (
        result["fragments_generated"] * fragment_mass
        + result["untracked_mass_kg"]
        == pytest.approx(result["ejecta_mass_kg"])
    )

def test_breakup_cratering_mass_conservation():
    # Setup a low-energy cratering impact (1g paint fleck vs 500kg target)
    result = calculate_breakup_fragments(
        target_mass_kg=500.0,
        projectile_mass_kg=0.001,
        impact_velocity_km_s=10.0,
        fragment_mass_kg=1.0
    )
    
    # Ejecta mass is 0.01kg, which cannot generate a full 1.0kg fragment
    assert result["fragments_generated"] == 0
    assert result["untracked_mass_kg"] == pytest.approx(result["ejecta_mass_kg"])


def test_oversized_step_conserves_mass_through_pmd_failure():
    sim = KesslerSimulation()
    sim.base_annual_launches = 0.0
    sim.pmd_failure_rate = 0.1
    for shell in sim.shells:
        shell.legacy_sat_count = 0.0
        shell.new_sat_count = 0.0
        shell.derelict_count = 0.0
        shell.debris_count = 0.0

    shell = next(s for s in sim.shells if s.altitude_km == 800)
    shell.legacy_sat_count = 10.0
    sim.step(dt_years=10.0)

    assert shell.sat_count == pytest.approx(0.0)
    assert shell.derelict_count == pytest.approx(1.0)
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-8)


def test_launch_pmd_derelict_and_atmospheric_reentry_lifecycle():
    sim = KesslerSimulation(data_path="missing-telemetry.json")
    sim.shells = [next(shell for shell in sim.shells if shell.altitude_km == 450)]
    sim.base_annual_launches = 10.0
    sim.launch_growth_rate = 0.0
    sim.sat_lifespan_years = 1.0
    sim.pmd_failure_rate = 1.0

    sim.step()
    sim.base_annual_launches = 0.0
    shell = sim.shells[0]
    assert shell.new_sat_count == pytest.approx(10.0)
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-8)

    sim.step()
    assert shell.new_sat_count == pytest.approx(0.0)
    assert shell.derelict_count == pytest.approx(10.0)
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-8)

    sim.step()
    assert shell.derelict_count == pytest.approx(0.0)
    assert sim.reentered_mass_kg == pytest.approx(2_500.0)
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-8)


def test_catastrophic_collision_consumes_target_and_projectile_mass_once():
    sim = KesslerSimulation()
    sim.base_annual_launches = 0.0
    for shell in sim.shells:
        shell.legacy_sat_count = 0.0
        shell.new_sat_count = 0.0
        shell.derelict_count = 0.0
        shell.debris_count = 0.0

    shell = next(s for s in sim.shells if s.altitude_km == 550)
    shell.legacy_sat_count = 1.0
    shell.debris_count = 20_000_000.0

    sim.step()

    assert shell.sat_count == pytest.approx(0.0)
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-8)


def test_invalid_step_parameters_raise_without_division_errors():
    sim = KesslerSimulation()
    with pytest.raises(ValueError, match="dt_years"):
        sim.step(dt_years=-1.0)

    sim.sat_lifespan_years = 0.0
    with pytest.raises(ValueError, match="sat_lifespan_years"):
        sim.step()

    sim = KesslerSimulation(data_path="missing-telemetry.json")
    sim.transiting_batches.append(
        TransitingDebrisBatch(
            count=1.0,
            origin_altitude_km=800.0,
            remaining_years=1.0,
            transit_duration_years=0.0,
        )
    )
    with pytest.raises(ValueError, match="Transit batches"):
        sim.step()


def test_breakup_rejects_unrepresentable_derived_energy():
    with pytest.raises(ValueError, match="kinetic energy"):
        calculate_breakup_fragments(
            target_mass_kg=250.0,
            projectile_mass_kg=1.0,
            impact_velocity_km_s=1e308,
        )


def test_subthreshold_transit_mass_is_not_discarded():
    sim = KesslerSimulation(data_path="missing-telemetry.json")
    sim.base_annual_launches = 0.0
    sim.transiting_batches = [
        # A small fractional cohort used to be purged before its scheduled reentry.
        TransitingDebrisBatch(
            count=0.05,
            origin_altitude_km=800.0,
            remaining_years=1.0,
            mass_kg=0.05,
        )
    ]
    sim.initial_mass_kg += 0.05

    sim.step()

    assert sim.total_active_transiting_debris == pytest.approx(0.0)
    assert sim.reentered_mass_kg == pytest.approx(0.05)
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-8)