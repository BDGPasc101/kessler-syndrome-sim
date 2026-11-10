from dataclasses import dataclass
from typing import List, Optional
import math
import sys

from src.physics.breakup import (
    HEAVY_CHUNK_CROSS_SECTION_M2,
    HEAVY_CHUNK_MASS_KG,
    SHRAPNEL_CROSS_SECTION_M2,
    SHRAPNEL_MASS_KG,
    WHIPPLE_SHIELD_LETHAL_ENERGY_J,
    calculate_kinetic_energy_j,
)

MU = 398600.4418
EARTH_RADIUS_KM = 6378.137
SATELLITE_MASS_KG = 250.0
FRAGMENT_MASS_KG = SHRAPNEL_MASS_KG
MAX_DRAG_DISPLACEMENT_KM = 10.0
MAX_COLLISION_EXPOSURE = 0.1
ROUTING_REFERENCE_CROSS_SECTION_M2 = 10.0
HIGHWAY_ALTITUDE_KM = 450.0
HIGHWAY_BASE_ALLOCATION = 0.60
HIGHWAY_SURVIVAL_THRESHOLD = 0.95
LAUNCH_REGIMES = (
    (300.0, 350.0, 400.0),
    (450.0,),
    (500.0, 550.0, 600.0),
    (650.0, 700.0, 750.0, 800.0, 850.0),
)
LAUNCH_REGIME_SHARES = (0.05, HIGHWAY_BASE_ALLOCATION, 0.25, 0.10)
HISTORICAL_DERELICT_ALTITUDES_KM = {700.0, 750.0, 800.0, 850.0}
HISTORICAL_DERELICTS_PER_SHELL = 200.0

@dataclass
class TransitingDebrisBatch:
    """Represents a cohort of debris in elliptical transfer orbits undergoing de-orbit decay."""
    count: float
    origin_altitude_km: float
    remaining_years: float
    transit_duration_years: float = 2.0
    initial_count: float = 0.0
    mass_kg: Optional[float] = None

    def __post_init__(self):
        if self.initial_count == 0.0:
            self.initial_count = self.count
        if self.mass_kg is None:
            self.mass_kg = self.count * FRAGMENT_MASS_KG

    @property
    def is_expired(self) -> bool:
        return self.remaining_years <= 0.0

    def occupancy_fraction(self, shell_altitude_km: float) -> float:
        """Fraction of orbit spent within a specific 50km shell."""
        perigee_altitude_km = 200.0
        if (
            shell_altitude_km < 400.0
            or not math.isfinite(shell_altitude_km)
            or not math.isfinite(self.origin_altitude_km)
        ):
            return 0.0

        perigee_radius_km = EARTH_RADIUS_KM + perigee_altitude_km
        apogee_radius_km = EARTH_RADIUS_KM + self.origin_altitude_km
        if apogee_radius_km <= perigee_radius_km:
            return 0.0

        shell_inner_radius_km = max(
            perigee_radius_km,
            EARTH_RADIUS_KM + shell_altitude_km - 25.0,
        )
        shell_outer_radius_km = min(
            apogee_radius_km,
            EARTH_RADIUS_KM + shell_altitude_km + 25.0,
        )
        if shell_inner_radius_km >= shell_outer_radius_km:
            return 0.0

        semi_major_axis_km = (perigee_radius_km + apogee_radius_km) / 2.0
        eccentricity = (
            (apogee_radius_km - perigee_radius_km)
            / (apogee_radius_km + perigee_radius_km)
        )
        if eccentricity <= 0.0:
            return 0.0

        def mean_anomaly(radius_km: float) -> float:
            cos_eccentric_anomaly = (
                1.0 - radius_km / semi_major_axis_km
            ) / eccentricity
            eccentric_anomaly = math.acos(
                max(-1.0, min(1.0, cos_eccentric_anomaly))
            )
            return eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly)

        return (
            mean_anomaly(shell_outer_radius_km)
            - mean_anomaly(shell_inner_radius_km)
        ) / math.pi


class OrbitalShell:
    def __init__(
        self,
        altitude_km: float,
        debris_count: float = 0.0,
        legacy_sats: float = 0.0,
        derelict_count: float = 0.0,
    ):
        self.altitude_km = float(altitude_km)
        self.debris_count = float(debris_count)
        self.derelict_count = float(derelict_count)
        self.heavy_debris_count = 0.0
        self.legacy_sat_count = float(legacy_sats)
        self.new_sat_count = 0.0
        self.satellite_mass_kg = self.sat_count * SATELLITE_MASS_KG
        self.derelict_mass_kg = self.derelict_count * SATELLITE_MASS_KG
        self.debris_mass_kg = self.debris_count * FRAGMENT_MASS_KG
        self.heavy_debris_mass_kg = (
            self.heavy_debris_count * HEAVY_CHUNK_MASS_KG
        )
        self.collision_lifespan_damage_years = 0.0
        self._mass_count_snapshot = self._get_count_snapshot()

        # Shell volume: V = 4/3 * pi * (r_outer^3 - r_inner^3)
        r_inner = EARTH_RADIUS_KM + self.altitude_km - 25.0
        r_outer = EARTH_RADIUS_KM + self.altitude_km + 25.0
        self.volume_km3 = (4.0 / 3.0) * math.pi * (r_outer**3 - r_inner**3)

    @property
    def sat_count(self) -> float:
        return self.legacy_sat_count + self.new_sat_count

    @property
    def total_mass_kg(self) -> float:
        return (
            self.satellite_mass_kg
            + self.derelict_mass_kg
            + self.debris_mass_kg
            + self.heavy_debris_mass_kg
        )

    def _get_count_snapshot(self) -> tuple:
        return (
            self.legacy_sat_count,
            self.new_sat_count,
            self.derelict_count,
            self.debris_count,
            self.heavy_debris_count,
        )

    def sync_mass_after_external_count_change(self) -> float:
        """Reset nominal mass pools when callers directly edit public counts."""
        snapshot = self._get_count_snapshot()
        if snapshot == self._mass_count_snapshot:
            return 0.0

        old_mass = self.total_mass_kg
        self.satellite_mass_kg = self.sat_count * SATELLITE_MASS_KG
        self.derelict_mass_kg = self.derelict_count * SATELLITE_MASS_KG
        self.debris_mass_kg = self.debris_count * FRAGMENT_MASS_KG
        self.heavy_debris_mass_kg = (
            self.heavy_debris_count * HEAVY_CHUNK_MASS_KG
        )
        self._mass_count_snapshot = snapshot
        return self.total_mass_kg - old_mass

    def update_mass_count_snapshot(self):
        self._mass_count_snapshot = self._get_count_snapshot()


class KesslerSimulation:
    def __init__(self, data_path: str = "data/celestrak_shells.json"):
        self.shells: List[OrbitalShell] = []
        self.transiting_batches: List[TransitingDebrisBatch] = []
        
        # Operational & mitigation settings
        self.base_annual_launches: float = 2000.0
        self.launch_growth_rate: float = 0.05
        self.sat_lifespan_years: float = 5.0
        self.pmd_failure_rate: float = 0.10    # 10% fail to de-orbit
        self.whipple_shield_lethal_energy_j = WHIPPLE_SHIELD_LETHAL_ENERGY_J
        self.shield_coverage_fraction: float = 0.3
        self.satellite_lifespan_damage_years_per_absorbed_impact = 0.1
        self.laser_removals_per_year: float = 0.0
        self.laser_transit_duration_years: float = 2.0
        self.enable_dynamic_highway_routing: bool = True
        self.launches_diverted_this_step: float = 0.0
        
        self.current_year: float = 0.0
        self._load_telemetry(data_path)
        self.initial_mass_kg = sum(shell.total_mass_kg for shell in self.shells)
        self.launched_mass_kg = 0.0
        self.reentered_mass_kg = 0.0

    def _load_telemetry(self, path: str):
        import os, json
        if os.path.exists(path):
            with open(path, "r") as f:
                data = json.load(f)
            for alt in range(300, 1050, 50):
                s_alt = str(alt)
                sats = data.get(s_alt, {}).get("satellites", 0.0)
                deb = data.get(s_alt, {}).get("debris", 0.0)
                derelicts = (
                    HISTORICAL_DERELICTS_PER_SHELL
                    if float(alt) in HISTORICAL_DERELICT_ALTITUDES_KM
                    else 0.0
                )
                self.shells.append(
                    OrbitalShell(
                        altitude_km=alt,
                        debris_count=deb,
                        legacy_sats=sats,
                        derelict_count=derelicts,
                    )
                )
        else:
            for alt in range(300, 1050, 50):
                derelicts = (
                    HISTORICAL_DERELICTS_PER_SHELL
                    if float(alt) in HISTORICAL_DERELICT_ALTITUDES_KM
                    else 0.0
                )
                self.shells.append(
                    OrbitalShell(altitude_km=alt, derelict_count=derelicts)
                )

    def get_transiting_debris_in_shell(self, altitude_km: float) -> float:
        """Calculates instantaneous transiting debris inside an altitude shell."""
        flux = 0.0
        for batch in self.transiting_batches:
            flux += batch.count * batch.occupancy_fraction(altitude_km)
        return flux

    @property
    def total_active_transiting_debris(self) -> float:
        """Total debris currently transiting in eccentric decay orbits."""
        return sum(b.count for b in self.transiting_batches)

    @property
    def total_debris_in_orbit(self) -> float:
        """Total circular debris plus transiting debris."""
        return (
            sum(s.debris_count + s.heavy_debris_count for s in self.shells)
            + self.total_active_transiting_debris
        )

    @property
    def total_mass_in_orbit_kg(self) -> float:
        """Mass currently represented by orbital populations and transit batches."""
        return sum(shell.total_mass_kg for shell in self.shells) + sum(
            batch.mass_kg for batch in self.transiting_batches
        )

    @property
    def mass_balance_error_kg(self) -> float:
        """Current orbital mass plus reentry sinks minus all initial/launched mass."""
        accounted_mass = self.initial_mass_kg + self.launched_mass_kg
        return self.total_mass_in_orbit_kg + self.reentered_mass_kg - accounted_mass

    def _validate_step_inputs(self, dt_years: float):
        finite_values = {
            "dt_years": dt_years,
            "sat_lifespan_years": self.sat_lifespan_years,
            "pmd_failure_rate": self.pmd_failure_rate,
            "whipple_shield_lethal_energy_j": self.whipple_shield_lethal_energy_j,
            "shield_coverage_fraction": self.shield_coverage_fraction,
            "satellite_lifespan_damage_years_per_absorbed_impact":
                self.satellite_lifespan_damage_years_per_absorbed_impact,
            "base_annual_launches": self.base_annual_launches,
            "launch_growth_rate": self.launch_growth_rate,
            "laser_removals_per_year": self.laser_removals_per_year,
            "laser_transit_duration_years": self.laser_transit_duration_years,
        }
        if not all(math.isfinite(value) for value in finite_values.values()):
            raise ValueError("Simulation step parameters must be finite.")
        if dt_years < 0:
            raise ValueError("dt_years must not be negative.")
        if dt_years == 0:
            return
        if self.sat_lifespan_years <= 0:
            raise ValueError("sat_lifespan_years must be greater than zero.")
        if not 0.0 <= self.pmd_failure_rate <= 1.0:
            raise ValueError("pmd_failure_rate must be between zero and one.")
        if self.whipple_shield_lethal_energy_j <= 0:
            raise ValueError("whipple_shield_lethal_energy_j must be greater than zero.")
        if not 0.0 <= self.shield_coverage_fraction <= 1.0:
            raise ValueError("shield_coverage_fraction must be between zero and one.")
        if self.satellite_lifespan_damage_years_per_absorbed_impact < 0:
            raise ValueError(
                "satellite_lifespan_damage_years_per_absorbed_impact "
                "must not be negative."
            )
        if self.base_annual_launches < 0:
            raise ValueError("base_annual_launches must not be negative.")
        if self.launch_growth_rate < -1.0:
            raise ValueError("launch_growth_rate must be at least -1.")
        if self.laser_removals_per_year < 0:
            raise ValueError("laser_removals_per_year must not be negative.")
        if self.laser_removals_per_year > 0 and self.laser_transit_duration_years <= 0:
            raise ValueError("laser_transit_duration_years must be greater than zero.")
        for batch in self.transiting_batches:
            if (
                not math.isfinite(batch.remaining_years)
                or not math.isfinite(batch.transit_duration_years)
                or not math.isfinite(batch.count)
                or not math.isfinite(batch.mass_kg)
                or batch.remaining_years < 0
                or batch.transit_duration_years <= 0
                or batch.count < 0
                or batch.mass_kg < 0
            ):
                raise ValueError("Transit batches must have finite, nonnegative state and positive durations.")
        if not math.isfinite(self.current_year) or self.current_year < 0:
            raise ValueError("current_year must be finite and nonnegative.")

    def _project_launches(self, dt_years: float) -> float:
        growth_factor = 1.0 + self.launch_growth_rate
        if self.base_annual_launches == 0:
            annual_launches = 0.0
        elif growth_factor == 0 and self.current_year > 0:
            annual_launches = 0.0
        elif growth_factor > 1:
            log_launches = (
                math.log(self.base_annual_launches)
                + self.current_year * math.log(growth_factor)
            )
            if log_launches > math.log(sys.float_info.max):
                raise ValueError("Annual launch count overflowed; adjust launch parameters.")
            annual_launches = math.exp(log_launches)
        else:
            annual_launches = self.base_annual_launches * (
                growth_factor ** self.current_year
            )

        launches = annual_launches * dt_years
        if not math.isfinite(launches) or not math.isfinite(self.current_year + dt_years):
            raise ValueError("Step would overflow launch count or simulation time.")
        return launches

    def _highway_diversion_share(self) -> float:
        if not self.enable_dynamic_highway_routing:
            return 0.0

        from src.physics.collisions import (
            calculate_impact_metrics,
            calculate_relative_velocity_km_s,
        )

        highway_shell = next(
            (shell for shell in self.shells
             if shell.altitude_km == HIGHWAY_ALTITUDE_KM),
            None,
        )
        if highway_shell is None:
            raise ValueError("450 km Highway shell is required for launch routing.")

        lethal_density = (
            highway_shell.heavy_debris_count + highway_shell.derelict_count
        ) / highway_shell.volume_km3
        risk = calculate_impact_metrics(
            lethal_density,
            ROUTING_REFERENCE_CROSS_SECTION_M2,
            self.sat_lifespan_years + 5.0,
            calculate_relative_velocity_km_s(HIGHWAY_ALTITUDE_KM),
            lethal_spatial_density_km3=lethal_density,
        )["probability_of_catastrophic_breach"]
        excess_risk = max(0.0, risk - (1.0 - HIGHWAY_SURVIVAL_THRESHOLD))
        return min(HIGHWAY_BASE_ALLOCATION, excess_risk)

    def _consume_projectiles(
        self,
        shell: OrbitalShell,
        amount: float,
        projectile_type: str = "all",
    ) -> tuple[float, float]:
        """Remove consumed shrapnel or heavy chunks from this shell and transit."""
        if projectile_type not in {"all", "shrapnel", "heavy"}:
            raise ValueError("projectile_type must be all, shrapnel, or heavy.")
        remaining = amount
        mass_consumed = 0.0
        if projectile_type in {"all", "shrapnel"}:
            light_consumed = min(
                remaining,
                shell.debris_count,
                shell.debris_mass_kg / FRAGMENT_MASS_KG,
            )
            shell.debris_count -= light_consumed
            shell.debris_mass_kg -= light_consumed * FRAGMENT_MASS_KG
            remaining -= light_consumed
            mass_consumed += light_consumed * FRAGMENT_MASS_KG

            for batch in self.transiting_batches:
                if remaining <= 0:
                    break
                occupancy = batch.occupancy_fraction(shell.altitude_km)
                if occupancy <= 0:
                    continue
                available = min(
                    batch.count * occupancy,
                    batch.mass_kg / FRAGMENT_MASS_KG,
                )
                consumed = min(remaining, available, batch.count)
                batch.count -= consumed
                batch.mass_kg -= consumed * FRAGMENT_MASS_KG
                remaining -= consumed
                mass_consumed += consumed * FRAGMENT_MASS_KG

        if projectile_type in {"all", "heavy"}:
            heavy_consumed = min(
                remaining,
                shell.heavy_debris_count,
                shell.heavy_debris_mass_kg / HEAVY_CHUNK_MASS_KG,
            )
            shell.heavy_debris_count -= heavy_consumed
            shell.heavy_debris_mass_kg -= heavy_consumed * HEAVY_CHUNK_MASS_KG
            remaining -= heavy_consumed
            mass_consumed += heavy_consumed * HEAVY_CHUNK_MASS_KG
        return amount - remaining, mass_consumed

    def _adaptive_substep_years(self, remaining_years: float) -> float:
        """Choose a step that resolves shell-to-shell drag and collision exposure."""
        from src.physics.collisions import (
            calculate_impact_metrics,
            calculate_relative_velocity_km_s,
        )
        from src.physics.drag import calculate_orbital_decay_rate

        max_substep_years = remaining_years
        max_collision_rate = 0.0

        for shell in self.shells:
            transiting_in_shell = self.get_transiting_debris_in_shell(
                shell.altitude_km
            )
            if shell.debris_count + transiting_in_shell > 0:
                fragment_rate = calculate_orbital_decay_rate(
                    shell.altitude_km,
                    mass_kg=FRAGMENT_MASS_KG,
                    cross_section_m2=SHRAPNEL_CROSS_SECTION_M2,
                )["decay_rate_km_year"]
                if fragment_rate > 0:
                    max_substep_years = min(
                        max_substep_years,
                        MAX_DRAG_DISPLACEMENT_KM / fragment_rate,
                    )
            if shell.heavy_debris_count > 0:
                heavy_rate = calculate_orbital_decay_rate(
                    shell.altitude_km,
                    mass_kg=HEAVY_CHUNK_MASS_KG,
                    cross_section_m2=HEAVY_CHUNK_CROSS_SECTION_M2,
                )["decay_rate_km_year"]
                if heavy_rate > 0:
                    max_substep_years = min(
                        max_substep_years,
                        MAX_DRAG_DISPLACEMENT_KM / heavy_rate,
                    )
            if shell.derelict_count > 0:
                derelict_rate = calculate_orbital_decay_rate(
                    shell.altitude_km, mass_kg=250.0, cross_section_m2=5.0
                )["decay_rate_km_year"]
                max_substep_years = min(
                    max_substep_years,
                    MAX_DRAG_DISPLACEMENT_KM / derelict_rate,
                )

            effective_debris = (
                shell.debris_count
                + shell.heavy_debris_count
                + transiting_in_shell
                + shell.derelict_count
            )
            relative_velocity_km_s = calculate_relative_velocity_km_s(
                shell.altitude_km
            )
            impact_metrics = calculate_impact_metrics(
                effective_debris / shell.volume_km3,
                5.0,
                remaining_years,
                relative_velocity_km_s,
            )
            max_collision_rate = max(
                max_collision_rate, impact_metrics["rate_per_year"]
            )

        if max_collision_rate > 0:
            max_substep_years = min(
                max_substep_years,
                MAX_COLLISION_EXPOSURE / max_collision_rate,
            )

        if max_substep_years <= 0 or not math.isfinite(max_substep_years):
            raise ValueError("Adaptive integration produced an invalid substep.")
        return max_substep_years

    def _integrate_collision_and_drag_substep(self, dt_years: float):
        from src.physics.collisions import (
            calculate_impact_metrics,
            calculate_relative_velocity_km_s,
        )
        from src.physics.drag import calculate_orbital_decay_rate
        from src.physics.breakup import calculate_breakup_fragments

        downward_debris_flux = 0.0
        downward_debris_mass_flux = 0.0
        downward_heavy_debris_flux = 0.0
        downward_heavy_debris_mass_flux = 0.0
        downward_derelict_flux = 0.0
        downward_derelict_mass_flux = 0.0

        for i in reversed(range(len(self.shells))):
            shell = self.shells[i]
            transiting_in_shell = self.get_transiting_debris_in_shell(
                shell.altitude_km
            )
            relative_velocity_km_s = calculate_relative_velocity_km_s(
                shell.altitude_km
            )
            shrapnel_metrics = calculate_impact_metrics(
                (shell.debris_count + transiting_in_shell) / shell.volume_km3,
                5.0,
                dt_years,
                relative_velocity_km_s,
            )
            heavy_metrics = calculate_impact_metrics(
                shell.heavy_debris_count / shell.volume_km3,
                5.0,
                dt_years,
                relative_velocity_km_s,
            )

            active_targets = shell.sat_count
            derelict_targets = shell.derelict_count
            active_target_mass = shell.satellite_mass_kg
            derelict_target_mass = shell.derelict_mass_kg
            target_capacity = active_targets + derelict_targets
            raw_shrapnel_events = target_capacity * -math.expm1(
                -shrapnel_metrics["rate_per_year"] * dt_years
            )
            transit_projectiles = sum(
                min(
                    batch.count * batch.occupancy_fraction(shell.altitude_km),
                    batch.mass_kg / FRAGMENT_MASS_KG,
                )
                for batch in self.transiting_batches
            )
            shrapnel_capacity = min(
                shell.debris_count,
                shell.debris_mass_kg / FRAGMENT_MASS_KG,
            ) + transit_projectiles
            heavy_capacity = min(
                shell.heavy_debris_count,
                shell.heavy_debris_mass_kg / HEAVY_CHUNK_MASS_KG,
            )
            raw_heavy_events = target_capacity * -math.expm1(
                -heavy_metrics["rate_per_year"] * dt_years
            )
            shrapnel_events = min(raw_shrapnel_events, shrapnel_capacity)
            heavy_events = min(raw_heavy_events, heavy_capacity)
            requested_events = shrapnel_events + heavy_events
            if requested_events > target_capacity and requested_events > 0:
                target_scale = target_capacity / requested_events
                shrapnel_events *= target_scale
                heavy_events *= target_scale

            shrapnel_events, shrapnel_mass_consumed = self._consume_projectiles(
                shell, shrapnel_events, projectile_type="shrapnel"
            )
            heavy_events, heavy_mass_consumed = self._consume_projectiles(
                shell, heavy_events, projectile_type="heavy"
            )

            if target_capacity > 0:
                active_target_fraction = active_targets / target_capacity
                active_shrapnel_events = shrapnel_events * active_target_fraction
                derelict_shrapnel_events = (
                    shrapnel_events - active_shrapnel_events
                )
                active_heavy_events = heavy_events * active_target_fraction
                derelict_heavy_events = (
                    heavy_events - active_heavy_events
                )

                def apply_collision(
                    target_count: float,
                    target_mass: float,
                    events: float,
                    projectile_mass: float,
                    projectile_is_heavy: bool,
                    shielded_target: bool,
                ):
                    if events <= 0 or target_count <= 0:
                        return 0.0, 0.0, 0.0, 0.0, 0.0

                    average_target_mass = target_mass / target_count
                    kinetic_energy_j = calculate_kinetic_energy_j(
                        projectile_mass, relative_velocity_km_s
                    )
                    if (
                        shielded_target
                        and not projectile_is_heavy
                        and kinetic_energy_j < self.whipple_shield_lethal_energy_j
                    ):
                        shielded_events = (
                            events * self.shield_coverage_fraction
                        )
                        mission_kill_events = events - shielded_events
                        mission_kill_mass = (
                            mission_kill_events * average_target_mass
                        )
                        shell.collision_lifespan_damage_years += (
                            shielded_events
                            * self.satellite_lifespan_damage_years_per_absorbed_impact
                        )
                        return (
                            0.0,
                            0.0,
                            shielded_events * projectile_mass,
                            mission_kill_events,
                            mission_kill_mass,
                        )

                    breakup_data = calculate_breakup_fragments(
                        target_mass_kg=average_target_mass,
                        projectile_mass_kg=projectile_mass,
                        impact_velocity_km_s=relative_velocity_km_s,
                        fragment_mass_kg=FRAGMENT_MASS_KG,
                        force_catastrophic=projectile_is_heavy,
                    )
                    shell.debris_count += events * breakup_data["shrapnel_count"]
                    shell.debris_mass_kg += events * (
                        breakup_data["shrapnel_mass_kg"]
                        + breakup_data["untracked_mass_kg"]
                    )
                    shell.heavy_debris_count += (
                        events * breakup_data["heavy_chunk_count"]
                    )
                    shell.heavy_debris_mass_kg += (
                        events * breakup_data["heavy_chunk_mass_kg"]
                    )

                    if breakup_data["is_catastrophic"]:
                        destroyed_count = min(events, target_count)
                        destroyed_mass = destroyed_count * average_target_mass
                        return destroyed_count, destroyed_mass, 0.0, 0.0, 0.0

                    retained_mass_change = (
                        projectile_mass - breakup_data["ejecta_mass_kg"]
                    ) * events
                    return 0.0, 0.0, retained_mass_change, 0.0, 0.0

                shrapnel_mass_per_projectile = (
                    shrapnel_mass_consumed / shrapnel_events
                    if shrapnel_events > 0 else 0.0
                )
                heavy_mass_per_projectile = (
                    heavy_mass_consumed / heavy_events
                    if heavy_events > 0 else 0.0
                )
                active_destroyed = 0.0
                active_mass_removed = 0.0
                active_mission_kills = 0.0
                active_mission_kill_mass = 0.0
                derelict_destroyed = 0.0
                derelict_mass_removed = 0.0

                for events, projectile_mass, is_heavy, targets_are_active in (
                    (
                        active_shrapnel_events,
                        shrapnel_mass_per_projectile,
                        False,
                        True,
                    ),
                    (
                        derelict_shrapnel_events,
                        shrapnel_mass_per_projectile,
                        False,
                        False,
                    ),
                    (
                        active_heavy_events,
                        heavy_mass_per_projectile,
                        True,
                        True,
                    ),
                    (
                        derelict_heavy_events,
                        heavy_mass_per_projectile,
                        True,
                        False,
                    ),
                ):
                    target_count = active_targets if targets_are_active else derelict_targets
                    target_mass = (
                        active_target_mass
                        if targets_are_active
                        else derelict_target_mass
                    )
                    (
                        destroyed,
                        mass_removed,
                        mass_change,
                        mission_kills,
                        mission_kill_mass,
                    ) = apply_collision(
                        target_count,
                        target_mass,
                        events,
                        projectile_mass,
                        is_heavy,
                        targets_are_active,
                    )
                    if targets_are_active:
                        active_destroyed += destroyed
                        active_mass_removed += mass_removed
                        active_mission_kills += mission_kills
                        active_mission_kill_mass += mission_kill_mass
                        shell.satellite_mass_kg += mass_change
                        shell.derelict_count += mission_kills
                        shell.derelict_mass_kg += (
                            mission_kill_mass
                            + mission_kills * projectile_mass
                        )
                    else:
                        derelict_destroyed += destroyed
                        derelict_mass_removed += mass_removed
                        shell.derelict_mass_kg += mass_change

                active_retired = active_destroyed + active_mission_kills
                if active_retired > 0:
                    legacy_ratio = shell.legacy_sat_count / active_targets
                    legacy_destroyed = active_retired * legacy_ratio
                    shell.legacy_sat_count -= legacy_destroyed
                    shell.new_sat_count -= active_retired - legacy_destroyed
                shell.satellite_mass_kg -= (
                    active_mass_removed + active_mission_kill_mass
                )
                shell.derelict_count -= derelict_destroyed
                shell.derelict_mass_kg -= derelict_mass_removed

            fragment_rate = calculate_orbital_decay_rate(
                shell.altitude_km,
                mass_kg=FRAGMENT_MASS_KG,
                cross_section_m2=SHRAPNEL_CROSS_SECTION_M2,
            )["decay_rate_km_year"]
            decay_frac_frag = -math.expm1(
                -(fragment_rate * dt_years) / 50.0
            )
            debris_leaving = shell.debris_count * decay_frac_frag
            debris_mass_leaving = shell.debris_mass_kg * decay_frac_frag
            shell.debris_count -= debris_leaving
            shell.debris_mass_kg -= debris_mass_leaving

            heavy_rate = calculate_orbital_decay_rate(
                shell.altitude_km,
                mass_kg=HEAVY_CHUNK_MASS_KG,
                cross_section_m2=HEAVY_CHUNK_CROSS_SECTION_M2,
            )["decay_rate_km_year"]
            decay_frac_heavy = -math.expm1(
                -(heavy_rate * dt_years) / 50.0
            )
            heavy_debris_leaving = (
                shell.heavy_debris_count * decay_frac_heavy
            )
            heavy_debris_mass_leaving = (
                shell.heavy_debris_mass_kg * decay_frac_heavy
            )
            shell.heavy_debris_count -= heavy_debris_leaving
            shell.heavy_debris_mass_kg -= heavy_debris_mass_leaving

            derelict_rate = calculate_orbital_decay_rate(
                shell.altitude_km, mass_kg=250.0, cross_section_m2=5.0
            )["decay_rate_km_year"]
            decay_frac_derelict = -math.expm1(
                -(derelict_rate * dt_years) / 50.0
            )
            derelicts_leaving = shell.derelict_count * decay_frac_derelict
            derelict_mass_leaving = (
                shell.derelict_mass_kg * decay_frac_derelict
            )
            shell.derelict_count -= derelicts_leaving
            shell.derelict_mass_kg -= derelict_mass_leaving

            if i + 1 < len(self.shells):
                shell.debris_count += downward_debris_flux
                shell.debris_mass_kg += downward_debris_mass_flux
                shell.heavy_debris_count += downward_heavy_debris_flux
                shell.heavy_debris_mass_kg += downward_heavy_debris_mass_flux
                shell.derelict_count += downward_derelict_flux
                shell.derelict_mass_kg += downward_derelict_mass_flux

            downward_debris_flux = debris_leaving
            downward_debris_mass_flux = debris_mass_leaving
            downward_heavy_debris_flux = heavy_debris_leaving
            downward_heavy_debris_mass_flux = heavy_debris_mass_leaving
            downward_derelict_flux = derelicts_leaving
            downward_derelict_mass_flux = derelict_mass_leaving

            if i == 0:
                self.reentered_mass_kg += (
                    debris_mass_leaving
                    + heavy_debris_mass_leaving
                    + derelict_mass_leaving
                )

    def step(self, dt_years: float = 1.0):
        self._validate_step_inputs(dt_years)
        if dt_years == 0:
            return
        launch_count = self._project_launches(dt_years)
        diverted_share = self._highway_diversion_share()
        self.launches_diverted_this_step = launch_count * diverted_share

        for shell in self.shells:
            self.initial_mass_kg += shell.sync_mass_after_external_count_change()

        # 1. Update transiting batches (perigee drag decay)
        for batch in self.transiting_batches:
            elapsed_fraction = (
                1.0 if batch.remaining_years <= 0
                else min(1.0, dt_years / batch.remaining_years)
            )
            mass_before = batch.mass_kg
            batch.count *= 1.0 - elapsed_fraction
            batch.mass_kg *= 1.0 - elapsed_fraction
            batch.remaining_years = max(0.0, batch.remaining_years - dt_years)
            self.reentered_mass_kg += mass_before - batch.mass_kg

        # Purge batches that have completed atmospheric reentry
        self.transiting_batches = [
            batch for batch in self.transiting_batches if not batch.is_expired
        ]

        # 2. Risk-Based & Mass-Constrained Laser Broom Targeting
        if self.laser_removals_per_year > 0:
            # Assume 10% of tracked debris are massive objects (rocket bodies, dead buses)
            # with too much inertia for a laser to perturb. Only 90% are targetable fragments.
            TARGETABLE_FRACTION = 0.90

            # Evaluate threat level: prioritize shells where debris and active satellites mix
            def risk_score(shell):
                targetable_debris = shell.debris_count * TARGETABLE_FRACTION
                # If a shell has active sats, risk scales by potential collisions.
                # If barren, fallback to cleaning the highest debris count (weighted low).
                return targetable_debris * shell.sat_count if shell.sat_count > 0 else targetable_debris * 1e-6

            target_shell = max(self.shells, key=risk_score)

            targetable_count = target_shell.debris_count * TARGETABLE_FRACTION
            removals = min(targetable_count, self.laser_removals_per_year * dt_years)

            if removals > 0:
                removed_mass = (
                    target_shell.debris_mass_kg
                    * removals
                    / target_shell.debris_count
                )
                target_shell.debris_count -= removals
                target_shell.debris_mass_kg -= removed_mass

                # If targeted above 400km, object enters elliptical transit crossing lower shells
                if target_shell.altitude_km > 400.0:
                    self.transiting_batches.append(
                        TransitingDebrisBatch(
                            count=removals,
                            origin_altitude_km=target_shell.altitude_km,
                            remaining_years=self.laser_transit_duration_years,
                            transit_duration_years=self.laser_transit_duration_years,
                            mass_kg=removed_mass,
                        )
                    )
                else:
                    self.reentered_mass_kg += removed_mass

        # 3. Integrate collision and drag changes without altering outer-step events.
        elapsed_years = 0.0
        while elapsed_years < dt_years:
            remaining_years = dt_years - elapsed_years
            substep_years = self._adaptive_substep_years(remaining_years)
            if elapsed_years + substep_years == elapsed_years:
                raise ValueError("Adaptive integration step is below time precision.")
            self._integrate_collision_and_drag_substep(substep_years)
            if substep_years == remaining_years:
                elapsed_years = dt_years
            else:
                elapsed_years += substep_years

        # 4. Satellite lifespan turnover remains once per requested step.
        for shell in self.shells:
            # Convert accumulated satellite-years of wear to aggregate early retirements.
            lifespan_damage_fraction = (
                shell.collision_lifespan_damage_years
                / (self.sat_lifespan_years * shell.sat_count)
                if shell.sat_count > 0 else 0.0
            )
            loss_fraction = min(
                1.0,
                dt_years / self.sat_lifespan_years
                + lifespan_damage_fraction,
            )
            legacy_loss = shell.legacy_sat_count * loss_fraction
            new_loss = shell.new_sat_count * loss_fraction
            total_satellite_loss = legacy_loss + new_loss
            shell.collision_lifespan_damage_years = 0.0
            satellite_mass_per_object = (
                shell.satellite_mass_kg / shell.sat_count
                if shell.sat_count > 0 else SATELLITE_MASS_KG
            )
            lost_satellite_mass = total_satellite_loss * satellite_mass_per_object
            shell.legacy_sat_count -= legacy_loss
            shell.new_sat_count -= new_loss
            shell.satellite_mass_kg -= lost_satellite_mass

            # PMD Failure Logic: Failed sats explicitly become massive Derelicts
            failed_sats = total_satellite_loss * self.pmd_failure_rate
            failed_satellite_mass = lost_satellite_mass * self.pmd_failure_rate
            shell.derelict_count += failed_sats
            shell.derelict_mass_kg += failed_satellite_mass
            self.reentered_mass_kg += lost_satellite_mass - failed_satellite_mass

        # 5. Inject launches across the orbital regimes.
        regime_shares = list(LAUNCH_REGIME_SHARES)
        regime_shares[1] -= diverted_share
        regime_shares[2] += diverted_share
        shells_by_altitude = {shell.altitude_km: shell for shell in self.shells}
        launch_allocations = []
        total_available_share = 0.0
        for altitudes, share in zip(LAUNCH_REGIMES, regime_shares):
            target_shells = [
                shells_by_altitude[altitude]
                for altitude in altitudes
                if altitude in shells_by_altitude
            ]
            if not target_shells:
                continue
            total_available_share += share
            launch_allocations.extend(
                (shell, share / len(target_shells)) for shell in target_shells
            )

        if total_available_share > 0.0:
            for shell, share in launch_allocations:
                per_shell = launch_count * share / total_available_share
                shell.new_sat_count += per_shell
                launched_mass = per_shell * SATELLITE_MASS_KG
                shell.satellite_mass_kg += launched_mass
                self.launched_mass_kg += launched_mass

        self.current_year += dt_years
        for shell in self.shells:
            shell.update_mass_count_snapshot()