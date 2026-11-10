from dataclasses import dataclass, field
from typing import Dict, List, Optional
import math
import sys

from src.physics.breakup import (
    HEAVY_CHUNK_CROSS_SECTION_M2,
    HEAVY_CHUNK_MASS_KG,
    SHRAPNEL_CROSS_SECTION_M2,
    SHRAPNEL_MASS_KG,
    WHIPPLE_SHIELD_LETHAL_ENERGY_J,
    EJECTA_SHELL_SPREAD_WEIGHTS,
    calculate_kinetic_energy_j,
    validate_ejecta_shell_spread_weights,
)

MU = 398600.4418
EARTH_RADIUS_KM = 6378.137
SATELLITE_MASS_KG = 250.0
SATELLITE_CROSS_SECTION_M2 = 10.0
DEFAULT_CAMS_EFFICIENCY = 0.90
DEFAULT_CAM_BURNS_PER_SATELLITE = 75.0
DEFAULT_CAM_SCREENING_RADIUS_KM = 0.75
FRAGMENT_MASS_KG = SHRAPNEL_MASS_KG
MAX_DRAG_DISPLACEMENT_KM = 10.0
MAX_COLLISION_EXPOSURE = 0.1
MAX_SUBSTEPS_PER_YEAR = 25
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
    _occupancy_cache: Dict[float, float] = field(
        default_factory=dict, init=False, repr=False
    )

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
        cached_fraction = self._occupancy_cache.get(shell_altitude_km)
        if cached_fraction is not None:
            return cached_fraction

        perigee_altitude_km = 200.0
        if (
            shell_altitude_km < 400.0
            or not math.isfinite(shell_altitude_km)
            or not math.isfinite(self.origin_altitude_km)
        ):
            self._occupancy_cache[shell_altitude_km] = 0.0
            return 0.0

        perigee_radius_km = EARTH_RADIUS_KM + perigee_altitude_km
        apogee_radius_km = EARTH_RADIUS_KM + self.origin_altitude_km
        if apogee_radius_km <= perigee_radius_km:
            self._occupancy_cache[shell_altitude_km] = 0.0
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
            self._occupancy_cache[shell_altitude_km] = 0.0
            return 0.0

        semi_major_axis_km = (perigee_radius_km + apogee_radius_km) / 2.0
        eccentricity = (
            (apogee_radius_km - perigee_radius_km)
            / (apogee_radius_km + perigee_radius_km)
        )
        if eccentricity <= 0.0:
            self._occupancy_cache[shell_altitude_km] = 0.0
            return 0.0

        def mean_anomaly(radius_km: float) -> float:
            cos_eccentric_anomaly = (
                1.0 - radius_km / semi_major_axis_km
            ) / eccentricity
            eccentric_anomaly = math.acos(
                max(-1.0, min(1.0, cos_eccentric_anomaly))
            )
            return eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly)

        occupancy = (
            mean_anomaly(shell_outer_radius_km)
            - mean_anomaly(shell_inner_radius_km)
        ) / math.pi
        self._occupancy_cache[shell_altitude_km] = occupancy
        return occupancy


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
        self.cam_burns_per_satellite = DEFAULT_CAM_BURNS_PER_SATELLITE
        self.cam_fueled_sat_count = self.sat_count
        self.cam_burns_remaining = (
            self.sat_count * self.cam_burns_per_satellite
        )
        self._cam_sat_count_snapshot = self.sat_count
        self._mass_count_snapshot = self._get_count_snapshot()

        # Shell volume: V = 4/3 * pi * (r_outer^3 - r_inner^3)
        self.inner_radius_km = EARTH_RADIUS_KM + self.altitude_km - 25.0
        self.outer_radius_km = EARTH_RADIUS_KM + self.altitude_km + 25.0
        self.volume_km3 = (4.0 / 3.0) * math.pi * (
            self.outer_radius_km**3 - self.inner_radius_km**3
        )
        self.orbital_radius_km = EARTH_RADIUS_KM + self.altitude_km
        self.circular_orbital_velocity_km_s = math.sqrt(
            MU / self.orbital_radius_km
        )
        from src.physics.collisions import calculate_relative_velocity_km_s
        self.relative_velocity_km_s = calculate_relative_velocity_km_s(
            self.altitude_km
        )

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
        self.reconcile_cam_satellite_count(self._cam_sat_count_snapshot)
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
        self._cam_sat_count_snapshot = self.sat_count

    @property
    def cam_depleted_sat_count(self) -> float:
        return max(0.0, self.sat_count - self.cam_fueled_sat_count)

    def reconcile_cam_satellite_count(self, previous_count: float) -> None:
        current_count = self.sat_count
        if current_count > previous_count:
            added = current_count - previous_count
            if self.cam_burns_per_satellite > 0.0:
                self.cam_fueled_sat_count += added
                self.cam_burns_remaining += (
                    added * self.cam_burns_per_satellite
                )
        elif current_count < previous_count and previous_count > 0.0:
            retention_fraction = current_count / previous_count
            self.cam_fueled_sat_count *= retention_fraction
            self.cam_burns_remaining *= retention_fraction

        self.cam_fueled_sat_count = min(
            current_count, max(0.0, self.cam_fueled_sat_count)
        )
        if self.cam_burns_remaining <= 1e-12:
            self.cam_burns_remaining = 0.0
            self.cam_fueled_sat_count = 0.0
        self._cam_sat_count_snapshot = current_count

    def add_cam_fueled_satellites(self, count: float) -> None:
        if not math.isfinite(count) or count < 0.0:
            raise ValueError("Added CAM cohort must be finite and nonnegative.")
        if self.cam_burns_per_satellite > 0.0:
            self.cam_fueled_sat_count += count
            self.cam_burns_remaining += count * self.cam_burns_per_satellite
        self._cam_sat_count_snapshot = self.sat_count


class KesslerSimulation:
    def __init__(self, data_path: str = "data/celestrak_shells.json"):
        self.shells: List[OrbitalShell] = []
        self.transiting_batches: List[TransitingDebrisBatch] = []
        
        # Operational & mitigation settings
        self.base_annual_launches: float = 2000.0
        self.launch_growth_rate: float = 0.05
        self.sat_lifespan_years: float = 5.0
        self.pmd_failure_rate: float = 0.10    # 10% fail to de-orbit
        self.cams_efficiency: float = DEFAULT_CAMS_EFFICIENCY
        self.cam_burns_per_satellite: float = DEFAULT_CAM_BURNS_PER_SATELLITE
        self.cam_screening_radius_km: float = DEFAULT_CAM_SCREENING_RADIUS_KM
        self.whipple_shield_lethal_energy_j = WHIPPLE_SHIELD_LETHAL_ENERGY_J
        self.shield_coverage_fraction: float = 0.3
        self.satellite_lifespan_damage_years_per_absorbed_impact = 0.1
        self.laser_removals_per_year: float = 0.0
        self.laser_transit_duration_years: float = 2.0
        self.enable_dynamic_highway_routing: bool = True
        self.launches_diverted_this_step: float = 0.0
        
        self.current_year: float = 0.0
        self._load_telemetry(data_path)
        self.set_cam_burns_per_satellite(self.cam_burns_per_satellite)
        self.initial_mass_kg = sum(shell.total_mass_kg for shell in self.shells)
        self.launched_mass_kg = 0.0
        self.reentered_mass_kg = 0.0

    def set_cam_burns_per_satellite(self, burns: float) -> None:
        if not math.isfinite(burns) or burns < 0.0:
            raise ValueError("CAM burns per satellite must be finite and nonnegative.")
        self.cam_burns_per_satellite = burns
        for shell in self.shells:
            shell.cam_burns_per_satellite = burns
            shell.cam_burns_remaining = shell.cam_fueled_sat_count * burns
            if burns == 0.0:
                shell.cam_fueled_sat_count = 0.0

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
            "cams_efficiency": self.cams_efficiency,
            "cam_burns_per_satellite": self.cam_burns_per_satellite,
            "cam_screening_radius_km": self.cam_screening_radius_km,
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
        if not 0.0 <= self.cams_efficiency <= 1.0:
            raise ValueError("cams_efficiency must be between zero and one.")
        if self.cam_burns_per_satellite < 0.0:
            raise ValueError("cam_burns_per_satellite must not be negative.")
        if self.cam_screening_radius_km <= 0.0:
            raise ValueError("cam_screening_radius_km must be greater than zero.")
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
        for shell in self.shells:
            if (
                not math.isfinite(shell.cam_fueled_sat_count)
                or not math.isfinite(shell.cam_burns_remaining)
                or shell.cam_fueled_sat_count < 0.0
                or shell.cam_burns_remaining < 0.0
            ):
                raise ValueError(
                    "CAM cohort reserves must be finite and nonnegative."
                )
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

        from src.physics.collisions import calculate_impact_metrics

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
            highway_shell.relative_velocity_km_s,
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
        if projectile_type not in {"all", "shrapnel", "heavy", "derelict"}:
            raise ValueError(
                "projectile_type must be all, shrapnel, heavy, or derelict."
            )
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
        if projectile_type == "derelict":
            derelict_consumed = min(
                remaining,
                shell.derelict_count,
                shell.derelict_mass_kg / SATELLITE_MASS_KG,
            )
            shell.derelict_count -= derelict_consumed
            derelict_mass_consumed = derelict_consumed * SATELLITE_MASS_KG
            shell.derelict_mass_kg -= derelict_mass_consumed
            remaining -= derelict_consumed
            mass_consumed += derelict_mass_consumed
        return amount - remaining, mass_consumed

    def _cam_screening_rate_per_satellite(
        self, shell: OrbitalShell, relative_velocity_km_s: Optional[float] = None
    ) -> float:
        from src.physics.collisions import calculate_impact_metrics

        if relative_velocity_km_s is None:
            relative_velocity_km_s = shell.relative_velocity_km_s
        screening_area_m2 = math.pi * (
            self.cam_screening_radius_km * 1000.0
        ) ** 2
        macro_density = (
            shell.heavy_debris_count + shell.derelict_count
        ) / shell.volume_km3
        return calculate_impact_metrics(
            macro_density,
            screening_area_m2,
            1.0,
            relative_velocity_km_s,
        )["rate_per_year"]

    def _consume_cam_screenings(
        self, shell: OrbitalShell, dt_years: float, relative_velocity_km_s: float
    ) -> float:
        """Spend reserve on screened conjunctions and return effective CAM coverage."""
        active_satellites = shell.sat_count
        fueled_satellites = min(
            active_satellites, shell.cam_fueled_sat_count
        )
        if (
            active_satellites <= 0.0
            or fueled_satellites <= 0.0
            or shell.cam_burns_remaining <= 0.0
            or self.cams_efficiency <= 0.0
        ):
            return 0.0

        screening_rate = self._cam_screening_rate_per_satellite(
            shell, relative_velocity_km_s
        )
        requested_burns = (
            screening_rate * dt_years * fueled_satellites
        )
        if requested_burns <= 0.0:
            return 0.0

        covered_fraction = min(
            1.0, shell.cam_burns_remaining / requested_burns
        )
        effective_coverage = (
            self.cams_efficiency
            * fueled_satellites
            / active_satellites
            * covered_fraction
        )
        shell.cam_burns_remaining = max(
            0.0, shell.cam_burns_remaining - requested_burns
        )
        if shell.cam_burns_remaining <= 1e-12:
            shell.cam_burns_remaining = 0.0
            shell.cam_fueled_sat_count = 0.0
        return effective_coverage

    def _deposit_collision_ejecta(
        self, impact_shell: OrbitalShell, events: float, breakup_data: dict
    ) -> None:
        """Spread collision products into neighboring shells or the reentry sink."""
        if not math.isfinite(events) or events < 0.0:
            raise ValueError("Collision event count must be finite and nonnegative.")
        validate_ejecta_shell_spread_weights()
        try:
            impact_index = self.shells.index(impact_shell)
        except ValueError as exc:
            raise ValueError("Impact shell must belong to this simulation.") from exc

        shrapnel_count = events * breakup_data["shrapnel_count"]
        shrapnel_mass = events * (
            breakup_data["shrapnel_mass_kg"] + breakup_data["untracked_mass_kg"]
        )
        heavy_count = events * breakup_data["heavy_chunk_count"]
        heavy_mass = events * breakup_data["heavy_chunk_mass_kg"]

        for shell_offset, weight in EJECTA_SHELL_SPREAD_WEIGHTS:
            destination_index = impact_index + shell_offset
            if destination_index < 0:
                self.reentered_mass_kg += weight * (
                    shrapnel_mass + heavy_mass
                )
                continue

            destination = self.shells[min(destination_index, len(self.shells) - 1)]
            destination.debris_count += weight * shrapnel_count
            destination.debris_mass_kg += weight * shrapnel_mass
            destination.heavy_debris_count += weight * heavy_count
            destination.heavy_debris_mass_kg += weight * heavy_mass

    def _adaptive_substep_years(
        self, remaining_years: float, start_year: Optional[float] = None
    ) -> float:
        """Choose a step that resolves shell-to-shell drag and collision exposure."""
        from src.physics.collisions import calculate_impact_metrics
        from src.physics.drag import (
            calculate_orbital_decay_rate,
            maximum_solar_cycle_density_multiplier,
        )

        max_substep_years = remaining_years
        max_collision_rate = 0.0
        if start_year is None:
            start_year = self.current_year
        maximum_density_multiplier = maximum_solar_cycle_density_multiplier(
            start_year, remaining_years
        )

        for shell in self.shells:
            transiting_in_shell = self.get_transiting_debris_in_shell(
                shell.altitude_km
            )
            if shell.debris_count + transiting_in_shell > 0:
                fragment_rate = calculate_orbital_decay_rate(
                    shell.altitude_km,
                    mass_kg=FRAGMENT_MASS_KG,
                    cross_section_m2=SHRAPNEL_CROSS_SECTION_M2,
                    density_multiplier=maximum_density_multiplier,
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
                    density_multiplier=maximum_density_multiplier,
                )["decay_rate_km_year"]
                if heavy_rate > 0:
                    max_substep_years = min(
                        max_substep_years,
                        MAX_DRAG_DISPLACEMENT_KM / heavy_rate,
                    )
            if shell.derelict_count > 0:
                derelict_rate = calculate_orbital_decay_rate(
                    shell.altitude_km,
                    mass_kg=250.0,
                    cross_section_m2=5.0,
                    density_multiplier=maximum_density_multiplier,
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
            impact_metrics = calculate_impact_metrics(
                effective_debris / shell.volume_km3,
                SATELLITE_CROSS_SECTION_M2,
                remaining_years,
                shell.relative_velocity_km_s,
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

    def _integrate_collision_and_drag_substep(
        self, dt_years: float, simulation_year: Optional[float] = None
    ):
        from src.physics.collisions import calculate_impact_metrics
        from src.physics.drag import (
            calculate_orbital_decay_rate,
            solar_cycle_density_multiplier,
        )
        from src.physics.breakup import calculate_breakup_fragments

        if simulation_year is None:
            simulation_year = self.current_year + dt_years / 2.0
        density_multiplier = solar_cycle_density_multiplier(simulation_year)

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
            relative_velocity_km_s = shell.relative_velocity_km_s
            shrapnel_metrics = calculate_impact_metrics(
                (shell.debris_count + transiting_in_shell) / shell.volume_km3,
                SATELLITE_CROSS_SECTION_M2,
                dt_years,
                relative_velocity_km_s,
            )
            heavy_metrics = calculate_impact_metrics(
                shell.heavy_debris_count / shell.volume_km3,
                SATELLITE_CROSS_SECTION_M2,
                dt_years,
                relative_velocity_km_s,
            )
            derelict_metrics = calculate_impact_metrics(
                shell.derelict_count / shell.volume_km3,
                SATELLITE_CROSS_SECTION_M2,
                dt_years,
                relative_velocity_km_s,
            )

            active_targets = shell.sat_count
            derelict_targets = shell.derelict_count
            active_target_mass = shell.satellite_mass_kg
            derelict_target_mass = shell.derelict_mass_kg
            target_capacity = active_targets + derelict_targets
            cam_coverage = self._consume_cam_screenings(
                shell, dt_years, relative_velocity_km_s
            )
            shrapnel_probability = -math.expm1(
                -shrapnel_metrics["rate_per_year"] * dt_years
            )
            heavy_probability = -math.expm1(
                -heavy_metrics["rate_per_year"] * dt_years
            )
            derelict_probability = -math.expm1(
                -derelict_metrics["rate_per_year"] * dt_years
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
            derelict_capacity = min(
                shell.derelict_count,
                shell.derelict_mass_kg / SATELLITE_MASS_KG,
            )
            shrapnel_events = min(
                target_capacity * shrapnel_probability, shrapnel_capacity
            )
            active_fraction = (
                active_targets / target_capacity if target_capacity > 0 else 0.0
            )
            active_shrapnel_events = shrapnel_events * active_fraction
            derelict_shrapnel_events = (
                shrapnel_events - active_shrapnel_events
            )
            active_heavy_events = (
                active_targets
                * heavy_probability
                * (1.0 - cam_coverage)
            )
            derelict_heavy_events = derelict_targets * heavy_probability
            active_derelict_events = (
                active_targets
                * derelict_probability
                * (1.0 - cam_coverage)
            )
            derelict_pair_events = (
                0.5 * derelict_targets * derelict_probability
            )

            active_target_events = (
                active_shrapnel_events
                + active_heavy_events
                + active_derelict_events
            )
            if active_target_events > active_targets and active_target_events > 0:
                active_scale = active_targets / active_target_events
                active_shrapnel_events *= active_scale
                active_heavy_events *= active_scale
                active_derelict_events *= active_scale

            derelict_pool_use = (
                derelict_shrapnel_events
                + derelict_heavy_events
                + active_derelict_events
                + 2.0 * derelict_pair_events
            )
            if derelict_pool_use > derelict_capacity and derelict_pool_use > 0:
                derelict_scale = derelict_capacity / derelict_pool_use
                derelict_shrapnel_events *= derelict_scale
                derelict_heavy_events *= derelict_scale
                active_derelict_events *= derelict_scale
                derelict_pair_events *= derelict_scale

            heavy_events = active_heavy_events + derelict_heavy_events
            if heavy_events > heavy_capacity and heavy_events > 0:
                heavy_scale = heavy_capacity / heavy_events
                active_heavy_events *= heavy_scale
                derelict_heavy_events *= heavy_scale
            shrapnel_events = (
                active_shrapnel_events + derelict_shrapnel_events
            )

            shrapnel_events, shrapnel_mass_consumed = self._consume_projectiles(
                shell, shrapnel_events, projectile_type="shrapnel"
            )
            heavy_events, heavy_mass_consumed = self._consume_projectiles(
                shell, heavy_events, projectile_type="heavy"
            )
            if heavy_events < active_heavy_events + derelict_heavy_events:
                heavy_scale = heavy_events / (
                    active_heavy_events + derelict_heavy_events
                )
                active_heavy_events *= heavy_scale
                derelict_heavy_events *= heavy_scale
            derelict_events = active_derelict_events + derelict_pair_events
            derelict_events, derelict_mass_consumed = self._consume_projectiles(
                shell, derelict_events, projectile_type="derelict"
            )
            if derelict_events < active_derelict_events + derelict_pair_events:
                derelict_scale = (
                    derelict_events
                    / (active_derelict_events + derelict_pair_events)
                )
                active_derelict_events *= derelict_scale
                derelict_pair_events *= derelict_scale

            if target_capacity > 0:
                def apply_collision(
                    target_count: float,
                    target_mass: float,
                    events: float,
                    projectile_mass: float,
                    projectile_is_macroscopic: bool,
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
                        and not projectile_is_macroscopic
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
                        force_catastrophic=projectile_is_macroscopic,
                    )
                    self._deposit_collision_ejecta(
                        shell, events, breakup_data
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
                derelict_mass_per_projectile = (
                    derelict_mass_consumed / derelict_events
                    if derelict_events > 0 else 0.0
                )
                active_destroyed = 0.0
                active_mass_removed = 0.0
                active_mission_kills = 0.0
                active_mission_kill_mass = 0.0
                derelict_destroyed = 0.0
                derelict_mass_removed = 0.0

                for (
                    events,
                    projectile_mass,
                    projectile_is_macroscopic,
                    targets_are_active,
                ) in (
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
                    (
                        active_derelict_events,
                        derelict_mass_per_projectile,
                        True,
                        True,
                    ),
                    (
                        derelict_pair_events,
                        derelict_mass_per_projectile,
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
                        projectile_is_macroscopic,
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
                shell.reconcile_cam_satellite_count(active_targets)

            fragment_rate = calculate_orbital_decay_rate(
                shell.altitude_km,
                mass_kg=FRAGMENT_MASS_KG,
                cross_section_m2=SHRAPNEL_CROSS_SECTION_M2,
                density_multiplier=density_multiplier,
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
                density_multiplier=density_multiplier,
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
                shell.altitude_km,
                mass_kg=250.0,
                cross_section_m2=5.0,
                density_multiplier=density_multiplier,
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
        max_substeps = max(
            1, math.floor(dt_years * MAX_SUBSTEPS_PER_YEAR + 1e-12)
        )
        substeps_taken = 0
        while elapsed_years < dt_years:
            remaining_years = dt_years - elapsed_years
            substep_start_year = self.current_year + elapsed_years
            adaptive_step_years = self._adaptive_substep_years(
                remaining_years, substep_start_year
            )
            remaining_slots = max_substeps - substeps_taken
            minimum_step_years = remaining_years / remaining_slots
            substep_years = min(
                remaining_years,
                max(adaptive_step_years, minimum_step_years),
            )
            if elapsed_years + substep_years == elapsed_years:
                raise ValueError("Adaptive integration step is below time precision.")
            self._integrate_collision_and_drag_substep(
                substep_years,
                substep_start_year + substep_years / 2.0,
            )
            if substep_years == remaining_years:
                elapsed_years = dt_years
            else:
                elapsed_years += substep_years
            substeps_taken += 1

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
            fueled_cohort_before_loss = shell.cam_fueled_sat_count
            depleted_cohort_before_loss = shell.cam_depleted_sat_count
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
            fueled_retirements = fueled_cohort_before_loss * loss_fraction
            depleted_retirements = depleted_cohort_before_loss * loss_fraction
            failed_sats = (
                fueled_retirements * self.pmd_failure_rate
                + depleted_retirements
            )
            failed_satellite_mass = (
                lost_satellite_mass * failed_sats / total_satellite_loss
                if total_satellite_loss > 0.0 else 0.0
            )
            shell.cam_fueled_sat_count = max(
                0.0, shell.cam_fueled_sat_count - fueled_retirements
            )
            shell.cam_burns_remaining *= 1.0 - loss_fraction
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
                shell.add_cam_fueled_satellites(per_shell)
                self.launched_mass_kg += launched_mass

        self.current_year += dt_years
        for shell in self.shells:
            shell.update_mass_count_snapshot()