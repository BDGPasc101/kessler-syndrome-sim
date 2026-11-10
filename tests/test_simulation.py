import os
import json
import math
import pytest

from src.simulation import KesslerSimulation, TransitingDebrisBatch
from src.physics.collisions import (
    calculate_impact_metrics,
    calculate_relative_velocity_km_s,
)
from src.physics.breakup import (
    HEAVY_CHUNK_CROSS_SECTION_M2,
    HEAVY_CHUNK_MASS_KG,
    SHRAPNEL_CROSS_SECTION_M2,
    SHRAPNEL_MASS_KG,
    WHIPPLE_SHIELD_LETHAL_ENERGY_J,
    calculate_breakup_fragments,
    calculate_kinetic_energy_j,
)


def test_transiting_debris_occupancy_weights_time_near_apogee():
    batch = TransitingDebrisBatch(
        count=1.0,
        origin_altitude_km=800.0,
        remaining_years=1.0,
    )

    perigee_shell_fraction = batch.occupancy_fraction(400.0)
    apogee_shell_fraction = batch.occupancy_fraction(800.0)

    assert 0.0 < perigee_shell_fraction < 1.0
    assert 0.0 < apogee_shell_fraction < 1.0
    assert apogee_shell_fraction > perigee_shell_fraction


def test_transiting_debris_occupancy_clips_shells_to_orbit_bounds():
    batch = TransitingDebrisBatch(
        count=1.0,
        origin_altitude_km=430.0,
        remaining_years=1.0,
    )

    assert batch.occupancy_fraction(450.0) > 0.0
    assert batch.occupancy_fraction(500.0) == 0.0
    assert batch.occupancy_fraction(350.0) == 0.0
    assert (
        TransitingDebrisBatch(
            count=1.0,
            origin_altitude_km=200.0,
            remaining_years=1.0,
        ).occupancy_fraction(400.0)
        == 0.0
    )


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


def _disable_physics_during_launch_routing_test(sim, monkeypatch):
    monkeypatch.setattr(sim, "_adaptive_substep_years", lambda remaining: remaining)
    monkeypatch.setattr(sim, "_integrate_collision_and_drag_substep", lambda dt: None)


def test_initial_legacy_grid_derelicts_are_seeded_and_mass_accounted():
    sim = KesslerSimulation(data_path="missing-telemetry.json")

    derelict_shells = {
        shell.altitude_km: shell.derelict_count
        for shell in sim.shells
        if shell.altitude_km in (700.0, 750.0, 800.0, 850.0)
    }

    assert derelict_shells == {
        700.0: 200.0,
        750.0: 200.0,
        800.0: 200.0,
        850.0: 200.0,
    }
    assert {shell.altitude_km for shell in sim.shells} >= {300.0, 350.0}
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-8)


def test_launches_follow_base_highway_regime_shares(monkeypatch):
    sim = KesslerSimulation(data_path="missing-telemetry.json")
    sim.base_annual_launches = 1000.0
    sim.launch_growth_rate = 0.0
    _disable_physics_during_launch_routing_test(sim, monkeypatch)

    sim.step()

    launches_by_altitude = {
        shell.altitude_km: shell.new_sat_count for shell in sim.shells
    }
    assert sum(launches_by_altitude[alt] for alt in (300, 350, 400)) == pytest.approx(50)
    assert launches_by_altitude[450] == pytest.approx(600)
    assert sum(launches_by_altitude[alt] for alt in (500, 550, 600)) == pytest.approx(250)
    assert sum(launches_by_altitude[alt] for alt in (650, 700, 750, 800, 850)) == pytest.approx(100)
    assert sum(launches_by_altitude.values()) == pytest.approx(1000)
    assert sim.launches_diverted_this_step == pytest.approx(0.0)
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-7)


def test_dynamic_highway_routing_diverts_excess_lethal_risk(monkeypatch):
    sim = KesslerSimulation(data_path="missing-telemetry.json")
    sim.base_annual_launches = 1000.0
    sim.launch_growth_rate = 0.0
    highway_shell = next(shell for shell in sim.shells if shell.altitude_km == 450)
    highway_shell.heavy_debris_count = 100_000.0
    _disable_physics_during_launch_routing_test(sim, monkeypatch)

    from src.physics.collisions import calculate_impact_metrics

    lethal_density = highway_shell.heavy_debris_count / highway_shell.volume_km3
    expected_risk = calculate_impact_metrics(
        lethal_density,
        10.0,
        sim.sat_lifespan_years + 5.0,
        calculate_relative_velocity_km_s(450.0),
        lethal_spatial_density_km3=lethal_density,
    )["probability_of_catastrophic_breach"]
    expected_diverted_share = min(0.60, max(0.0, expected_risk - 0.05))
    assert expected_diverted_share > 0.0

    sim.step()

    assert sim.launches_diverted_this_step == pytest.approx(
        1000.0 * expected_diverted_share
    )
    assert highway_shell.new_sat_count == pytest.approx(
        1000.0 * (0.60 - expected_diverted_share)
    )
    assert sum(
        shell.new_sat_count for shell in sim.shells
        if shell.altitude_km in (500, 550, 600)
    ) == pytest.approx(1000.0 * (0.25 + expected_diverted_share))
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-7)


def test_disabled_dynamic_highway_routing_keeps_fixed_allocations(monkeypatch):
    sim = KesslerSimulation(data_path="missing-telemetry.json")
    sim.base_annual_launches = 1000.0
    sim.launch_growth_rate = 0.0
    sim.enable_dynamic_highway_routing = False
    highway_shell = next(shell for shell in sim.shells if shell.altitude_km == 450)
    highway_shell.heavy_debris_count = 100_000.0
    _disable_physics_during_launch_routing_test(sim, monkeypatch)

    sim.step()

    assert sim.launches_diverted_this_step == pytest.approx(0.0)
    assert highway_shell.new_sat_count == pytest.approx(600.0)
    assert sum(
        shell.new_sat_count for shell in sim.shells
        if shell.altitude_km in (500, 550, 600)
    ) == pytest.approx(250.0)


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


def test_adaptive_drag_substeps_transport_debris_through_lower_shells():
    sim = KesslerSimulation(data_path="missing-telemetry.json")
    sim.base_annual_launches = 0.0
    for shell in sim.shells:
        shell.legacy_sat_count = 0.0
        shell.new_sat_count = 0.0
        shell.debris_count = 0.0

    shell_450 = next(s for s in sim.shells if s.altitude_km == 450)
    shell_400 = next(s for s in sim.shells if s.altitude_km == 400)
    shell_450.debris_count = 1_000.0

    first_substep = sim._adaptive_substep_years(1.0)
    assert first_substep <= 10.0 / 99.6133

    sim.step(dt_years=1.0)

    assert shell_450.debris_count < 1_000.0
    assert shell_400.debris_count > 0.0
    assert sim.reentered_mass_kg > 0.0
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-7)


def test_dense_collision_states_limit_per_substep_exposure():
    sim = KesslerSimulation(data_path="missing-telemetry.json")
    sim.base_annual_launches = 0.0
    for shell in sim.shells:
        shell.legacy_sat_count = 0.0
        shell.new_sat_count = 0.0
        shell.derelict_count = 0.0
        shell.debris_count = 0.0

    shell = next(s for s in sim.shells if s.altitude_km == 800)
    shell.legacy_sat_count = 100.0
    shell.debris_count = 100_000_000.0
    rate = calculate_impact_metrics(
        shell.debris_count / shell.volume_km3,
        5.0,
        1.0,
        calculate_relative_velocity_km_s(shell.altitude_km),
    )["rate_per_year"]

    assert sim._adaptive_substep_years(1.0) <= 0.1 / rate

    sim.step(dt_years=1.0)

    assert all(
        math.isfinite(value) and value >= 0.0
        for value in (shell.sat_count, shell.debris_count, shell.debris_mass_kg)
    )
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-7)


def test_relative_velocity_uses_circular_speed_and_collision_geometry():
    velocity_400 = calculate_relative_velocity_km_s(400.0)
    velocity_1000 = calculate_relative_velocity_km_s(1000.0)

    assert velocity_400 == pytest.approx(
        math.sqrt(2.0 * 3.986e5 / (6371.0 + 400.0))
    )
    assert velocity_1000 == pytest.approx(
        math.sqrt(2.0 * 3.986e5 / (6371.0 + 1000.0))
    )
    assert velocity_400 > velocity_1000


def test_collision_breakup_uses_shell_relative_velocity(monkeypatch):
    import src.physics.breakup as breakup

    sim = KesslerSimulation(data_path="missing-telemetry.json")
    shell = next(s for s in sim.shells if s.altitude_km == 550)
    shell.legacy_sat_count = 100.0
    shell.heavy_debris_count = 1_000_000.0
    shell.sync_mass_after_external_count_change()

    original_calculate_breakup = breakup.calculate_breakup_fragments
    impact_velocities = []

    def record_impact_velocity(**kwargs):
        impact_velocities.append(kwargs["impact_velocity_km_s"])
        return original_calculate_breakup(**kwargs)

    monkeypatch.setattr(
        breakup, "calculate_breakup_fragments", record_impact_velocity
    )
    sim._integrate_collision_and_drag_substep(dt_years=1.0)

    assert impact_velocities
    assert impact_velocities == pytest.approx(
        [calculate_relative_velocity_km_s(shell.altitude_km)]
        * len(impact_velocities)
    )


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


def test_heavy_projectile_can_force_lethal_breakup_below_energy_threshold():
    result = calculate_breakup_fragments(
        target_mass_kg=500.0,
        projectile_mass_kg=HEAVY_CHUNK_MASS_KG,
        impact_velocity_km_s=0.01,
        force_catastrophic=True,
    )

    assert result["specific_energy_j_kg"] < 40_000.0
    assert result["is_catastrophic"]
    assert result["ejecta_mass_kg"] == pytest.approx(
        500.0 + HEAVY_CHUNK_MASS_KG
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
    assert pytest.approx(ratio, rel=1e-3) == 2.0
    assert 0.0 <= m_low["probability_of_impact"] <= 1.0


def test_impact_metrics_separate_lethal_breach_risk_from_any_impact():
    metrics = calculate_impact_metrics(
        spatial_density_km3=1e-5,
        cross_section_m2=10.0,
        time_horizon_years=5.0,
        relative_velocity_km_s=10.0,
        lethal_spatial_density_km3=2e-7,
    )
    no_lethal_projectiles = calculate_impact_metrics(
        spatial_density_km3=1e-5,
        cross_section_m2=10.0,
        time_horizon_years=5.0,
        relative_velocity_km_s=10.0,
        lethal_spatial_density_km3=0.0,
    )

    assert metrics["probability_of_impact"] > 0.0
    assert metrics["probability_of_catastrophic_breach"] > 0.0
    assert (
        metrics["probability_of_catastrophic_breach"]
        < metrics["probability_of_impact"]
    )
    assert no_lethal_projectiles["probability_of_catastrophic_breach"] == 0.0


def _make_shrapnel_collision_simulation():
    sim = KesslerSimulation(data_path="missing-telemetry.json")
    shell = next(s for s in sim.shells if s.altitude_km == 550)
    shell.legacy_sat_count = 100.0
    shell.satellite_mass_kg = 100.0 * 250.0
    shell.debris_count = 100_000_000.0
    shell.debris_mass_kg = shell.debris_count * SHRAPNEL_MASS_KG
    shell.update_mass_count_snapshot()
    sim.initial_mass_kg = sim.total_mass_in_orbit_kg
    return sim, shell


def test_subthreshold_shrapnel_is_absorbed_and_conserved():
    sim, shell = _make_shrapnel_collision_simulation()
    sim.whipple_shield_lethal_energy_j = WHIPPLE_SHIELD_LETHAL_ENERGY_J
    sim.shield_coverage_fraction = 1.0
    satellite_count_before = shell.sat_count
    satellite_mass_before = shell.satellite_mass_kg

    sim._integrate_collision_and_drag_substep(dt_years=1.0)

    assert shell.sat_count == pytest.approx(satellite_count_before)
    assert shell.satellite_mass_kg > satellite_mass_before
    assert shell.collision_lifespan_damage_years > 0.0
    assert (
        sim.total_mass_in_orbit_kg + sim.reentered_mass_kg
        == pytest.approx(sim.initial_mass_kg, rel=1e-10)
    )


def test_shield_coverage_defaults_to_thirty_percent_and_validates_range():
    assert (
        KesslerSimulation(data_path="missing-telemetry.json").shield_coverage_fraction
        == 0.3
    )

    for invalid_coverage, error_message in (
        (-0.01, "shield_coverage_fraction"),
        (1.01, "shield_coverage_fraction"),
        (float("nan"), "must be finite"),
    ):
        sim = KesslerSimulation(data_path="missing-telemetry.json")
        sim.shield_coverage_fraction = invalid_coverage
        with pytest.raises(ValueError, match=error_message):
            sim.step()


def test_subthreshold_shrapnel_mission_kills_exposed_area_without_breakup(
    monkeypatch,
):
    import src.physics.collisions as collisions
    import src.physics.drag as drag

    sim, shell = _make_shrapnel_collision_simulation()
    sim.whipple_shield_lethal_energy_j = 1e9
    monkeypatch.setattr(
        collisions,
        "calculate_impact_metrics",
        lambda *args, **kwargs: {"rate_per_year": 0.1},
    )
    monkeypatch.setattr(
        drag,
        "calculate_orbital_decay_rate",
        lambda *args, **kwargs: {"decay_rate_km_year": 0.0},
    )

    initial_active_count = shell.sat_count
    initial_satellite_mass = shell.satellite_mass_kg
    initial_derelict_mass = shell.derelict_mass_kg
    initial_shrapnel_count = shell.debris_count
    expected_events = initial_active_count * -math.expm1(-0.1)
    expected_mission_kills = expected_events * (
        1.0 - sim.shield_coverage_fraction
    )
    expected_shielded_events = expected_events * sim.shield_coverage_fraction

    sim._integrate_collision_and_drag_substep(dt_years=1.0)

    assert shell.sat_count == pytest.approx(
        initial_active_count - expected_mission_kills
    )
    assert shell.derelict_count == pytest.approx(expected_mission_kills)
    assert shell.collision_lifespan_damage_years == pytest.approx(
        expected_shielded_events
        * sim.satellite_lifespan_damage_years_per_absorbed_impact
    )
    assert shell.debris_count == pytest.approx(
        initial_shrapnel_count - expected_events
    )
    assert shell.heavy_debris_count == 0.0
    expected_target_mass_transferred = expected_mission_kills * (
        initial_satellite_mass / initial_active_count
    )
    assert shell.satellite_mass_kg == pytest.approx(
        initial_satellite_mass
        - expected_target_mass_transferred
        + expected_shielded_events * SHRAPNEL_MASS_KG
    )
    assert shell.derelict_mass_kg == pytest.approx(
        initial_derelict_mass
        + expected_target_mass_transferred
        + expected_mission_kills * SHRAPNEL_MASS_KG
    )
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-8)


@pytest.mark.parametrize(
    ("coverage", "expected_mission_kill_fraction"),
    [(0.0, 1.0), (1.0, 0.0)],
)
def test_shield_coverage_endpoints(
    coverage, expected_mission_kill_fraction, monkeypatch
):
    import src.physics.collisions as collisions
    import src.physics.drag as drag

    sim, shell = _make_shrapnel_collision_simulation()
    sim.whipple_shield_lethal_energy_j = 1e9
    sim.shield_coverage_fraction = coverage
    monkeypatch.setattr(
        collisions,
        "calculate_impact_metrics",
        lambda *args, **kwargs: {"rate_per_year": 0.1},
    )
    monkeypatch.setattr(
        drag,
        "calculate_orbital_decay_rate",
        lambda *args, **kwargs: {"decay_rate_km_year": 0.0},
    )
    expected_events = shell.sat_count * -math.expm1(-0.1)

    sim._integrate_collision_and_drag_substep(dt_years=1.0)

    assert shell.derelict_count == pytest.approx(
        expected_events * expected_mission_kill_fraction
    )
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-8)


def test_shrapnel_at_shield_threshold_is_not_absorbed():
    sim, shell = _make_shrapnel_collision_simulation()
    relative_velocity_km_s = calculate_relative_velocity_km_s(shell.altitude_km)
    satellite_mass_before = shell.satellite_mass_kg
    sim.whipple_shield_lethal_energy_j = calculate_kinetic_energy_j(
        SHRAPNEL_MASS_KG, relative_velocity_km_s
    )

    sim._integrate_collision_and_drag_substep(dt_years=1.0)

    assert shell.collision_lifespan_damage_years == 0.0
    assert shell.satellite_mass_kg < satellite_mass_before
    assert (
        sim.total_mass_in_orbit_kg + sim.reentered_mass_kg
        == pytest.approx(sim.initial_mass_kg, rel=1e-10)
    )


def test_above_threshold_shrapnel_is_not_absorbed():
    sim, shell = _make_shrapnel_collision_simulation()
    relative_velocity_km_s = calculate_relative_velocity_km_s(shell.altitude_km)
    satellite_mass_before = shell.satellite_mass_kg
    kinetic_energy_j = calculate_kinetic_energy_j(
        SHRAPNEL_MASS_KG, relative_velocity_km_s
    )
    sim.whipple_shield_lethal_energy_j = kinetic_energy_j * 0.9

    sim._integrate_collision_and_drag_substep(dt_years=1.0)

    assert shell.collision_lifespan_damage_years == 0.0
    assert shell.satellite_mass_kg < satellite_mass_before
    assert (
        sim.total_mass_in_orbit_kg + sim.reentered_mass_kg
        == pytest.approx(sim.initial_mass_kg, rel=1e-10)
    )


def test_heavy_debris_remains_catastrophically_lethal_and_conserved():
    sim = KesslerSimulation(data_path="missing-telemetry.json")
    shell = next(s for s in sim.shells if s.altitude_km == 550)
    shell.legacy_sat_count = 100.0
    shell.satellite_mass_kg = 100.0 * 250.0
    shell.heavy_debris_count = 10_000.0
    shell.heavy_debris_mass_kg = (
        shell.heavy_debris_count * HEAVY_CHUNK_MASS_KG
    )
    shell.update_mass_count_snapshot()
    sim.initial_mass_kg = sim.total_mass_in_orbit_kg
    satellite_count_before = shell.sat_count

    sim._integrate_collision_and_drag_substep(dt_years=1.0)

    assert shell.sat_count < satellite_count_before
    assert (
        sim.total_mass_in_orbit_kg + sim.reentered_mass_kg
        == pytest.approx(sim.initial_mass_kg, rel=1e-10)
    )


def test_mixed_shrapnel_and_heavy_impacts_conserve_mass():
    sim, shell = _make_shrapnel_collision_simulation()
    shell.heavy_debris_count = 1_000.0
    shell.heavy_debris_mass_kg = (
        shell.heavy_debris_count * HEAVY_CHUNK_MASS_KG
    )
    shell.update_mass_count_snapshot()
    sim.initial_mass_kg = sim.total_mass_in_orbit_kg

    sim._integrate_collision_and_drag_substep(dt_years=1.0)

    assert shell.collision_lifespan_damage_years > 0.0
    assert shell.sat_count < 100.0
    assert (
        sim.total_mass_in_orbit_kg + sim.reentered_mass_kg
        == pytest.approx(sim.initial_mass_kg, rel=1e-10)
    )


def test_absorbed_impact_lifespan_wear_causes_aggregate_early_retirement():
    def make_simulation(wear_years):
        sim = KesslerSimulation(data_path="missing-telemetry.json")
        sim.base_annual_launches = 0.0
        sim.sat_lifespan_years = 10_000.0
        sim.pmd_failure_rate = 0.0
        shell = next(s for s in sim.shells if s.altitude_km == 450)
        shell.legacy_sat_count = 100.0
        shell.satellite_mass_kg = 25_000.0
        shell.collision_lifespan_damage_years = wear_years
        shell.update_mass_count_snapshot()
        sim.initial_mass_kg = sim.total_mass_in_orbit_kg
        return sim, shell

    control, control_shell = make_simulation(0.0)
    damaged, damaged_shell = make_simulation(100.0)

    control.step(dt_years=0.01)
    damaged.step(dt_years=0.01)

    assert damaged_shell.sat_count < control_shell.sat_count
    assert damaged_shell.collision_lifespan_damage_years == 0.0
    assert damaged.mass_balance_error_kg == pytest.approx(0.0, abs=1e-8)


# ---------------------------------------------------------------------------
# Test 5: Laser Broom Targeted Removal
# ---------------------------------------------------------------------------
def test_laser_broom_targets_most_congested_shells_first():
    sim = KesslerSimulation()
    sim.base_annual_launches = 0
    sim.laser_removals_per_year = 5_000
    control = KesslerSimulation()
    control.base_annual_launches = 0

    for s in sim.shells:
        s.legacy_sat_count = 0.0
        s.new_sat_count = 0.0
        s.debris_count = 0.0
    for s in control.shells:
        s.legacy_sat_count = 0.0
        s.new_sat_count = 0.0
        s.debris_count = 0.0

    idx_550 = next(i for i, s in enumerate(sim.shells) if s.altitude_km == 550)
    idx_800 = next(i for i, s in enumerate(sim.shells) if s.altitude_km == 800)

    sim.shells[idx_550].debris_count = 1_000.0
    sim.shells[idx_800].debris_count = 20_000.0
    control.shells[idx_550].debris_count = 1_000.0
    control.shells[idx_800].debris_count = 20_000.0

    sim.step(dt_years=1.0)
    control.step(dt_years=1.0)

    assert sim.shells[idx_800].debris_count <= 15_000.0, "Laser failed to target 800km shell."
    assert sim.shells[idx_550].debris_count == pytest.approx(
        control.shells[idx_550].debris_count
    ), "Laser incorrectly deducted shots from 550km shell."


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

def test_breakup_mass_conservation_limit():
    result = calculate_breakup_fragments(
        target_mass_kg=250.0,
        projectile_mass_kg=1.0,
        impact_velocity_km_s=10.0,
    )

    assert result["shrapnel_count"] > 1_000
    assert result["heavy_chunk_count"] == 20
    assert (
        result["shrapnel_count"] * SHRAPNEL_MASS_KG
        + result["heavy_chunk_count"] * HEAVY_CHUNK_MASS_KG
        + result["untracked_mass_kg"]
        == pytest.approx(result["ejecta_mass_kg"])
    )

def test_breakup_cratering_mass_conservation():
    # Setup a low-energy cratering impact (1g paint fleck vs 500kg target)
    result = calculate_breakup_fragments(
        target_mass_kg=500.0,
        projectile_mass_kg=0.001,
        impact_velocity_km_s=10.0
    )
    
    assert result["shrapnel_count"] > 0
    assert result["heavy_chunk_count"] == 0
    assert (
        result["shrapnel_count"] * SHRAPNEL_MASS_KG
        + result["untracked_mass_kg"]
        == pytest.approx(result["ejecta_mass_kg"])
    )


def test_heavy_debris_has_distinct_drag_and_transports_mass_between_shells():
    from src.physics.drag import calculate_orbital_decay_rate

    sim = KesslerSimulation(data_path="missing-telemetry.json")
    sim.base_annual_launches = 0.0
    for shell in sim.shells:
        shell.legacy_sat_count = 0.0
        shell.new_sat_count = 0.0
        shell.debris_count = 0.0
        shell.heavy_debris_count = 0.0

    shell_450 = next(s for s in sim.shells if s.altitude_km == 450)
    shell_400 = next(s for s in sim.shells if s.altitude_km == 400)
    shell_450.heavy_debris_count = 10.0
    light_drag = calculate_orbital_decay_rate(
        450.0, SHRAPNEL_MASS_KG, SHRAPNEL_CROSS_SECTION_M2
    )["decay_rate_km_year"]
    heavy_drag = calculate_orbital_decay_rate(
        450.0, HEAVY_CHUNK_MASS_KG, HEAVY_CHUNK_CROSS_SECTION_M2
    )["decay_rate_km_year"]
    assert heavy_drag < light_drag

    sim.step()

    assert shell_450.heavy_debris_count < 10.0
    assert shell_400.heavy_debris_count > 0.0
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-7)


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
    sim.initial_mass_kg = sum(shell.total_mass_kg for shell in sim.shells)
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
    assert 0.0 < shell.derelict_count < 10.0
    assert sim.reentered_mass_kg > 0.0
    assert sim.mass_balance_error_kg == pytest.approx(0.0, abs=1e-8)


def test_catastrophic_collision_consumes_target_and_projectile_mass_once():
    sim = KesslerSimulation()
    sim.base_annual_launches = 0.0
    for shell in sim.shells:
        shell.legacy_sat_count = 0.0
        shell.new_sat_count = 0.0
        shell.derelict_count = 0.0
        shell.debris_count = 0.0
        shell.heavy_debris_count = 0.0

    shell = next(s for s in sim.shells if s.altitude_km == 550)
    shell.legacy_sat_count = 1.0
    shell.heavy_debris_count = 1_000_000.0

    sim.step()

    assert 0.0 < shell.sat_count < 0.8
    assert shell.debris_count > 0.0
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