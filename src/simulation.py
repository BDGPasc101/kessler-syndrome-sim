from dataclasses import dataclass
from typing import List, Optional
import math
import sys

MU = 398600.4418
EARTH_RADIUS_KM = 6378.137
SATELLITE_MASS_KG = 250.0
FRAGMENT_MASS_KG = 1.0
MAX_DRAG_DISPLACEMENT_KM = 10.0
MAX_COLLISION_EXPOSURE = 0.1

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
        if 400.0 <= shell_altitude_km <= self.origin_altitude_km:
            span = self.origin_altitude_km - 200.0  # Perigee at ~200km
            return 50.0 / span if span > 0 else 0.0
        return 0.0


class OrbitalShell:
    def __init__(self, altitude_km: float, debris_count: float = 0.0, legacy_sats: float = 0.0):
        self.altitude_km = float(altitude_km)
        self.debris_count = float(debris_count)
        self.derelict_count = 0.0    # NEW: Track dead satellite buses
        self.legacy_sat_count = float(legacy_sats)
        self.new_sat_count = 0.0
        self.satellite_mass_kg = self.sat_count * SATELLITE_MASS_KG
        self.derelict_mass_kg = self.derelict_count * SATELLITE_MASS_KG
        self.debris_mass_kg = self.debris_count * FRAGMENT_MASS_KG
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
        return self.satellite_mass_kg + self.derelict_mass_kg + self.debris_mass_kg

    def _get_count_snapshot(self) -> tuple:
        return (
            self.legacy_sat_count,
            self.new_sat_count,
            self.derelict_count,
            self.debris_count,
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
        self.laser_removals_per_year: float = 0.0
        self.laser_transit_duration_years: float = 2.0
        
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
            for alt in range(400, 1050, 50):
                s_alt = str(alt)
                sats = data.get(s_alt, {}).get("satellites", 0.0)
                deb = data.get(s_alt, {}).get("debris", 0.0)
                self.shells.append(OrbitalShell(altitude_km=alt, debris_count=deb, legacy_sats=sats))
        else:
            for alt in range(400, 1050, 50):
                self.shells.append(OrbitalShell(altitude_km=alt))

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
        return sum(s.debris_count for s in self.shells) + self.total_active_transiting_debris

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

    def _consume_projectiles(self, shell: OrbitalShell, amount: float) -> float:
        """Remove consumed 1 kg projectiles from this shell and its transit cohorts."""
        remaining = amount
        consumed = min(
            remaining,
            shell.debris_count,
            shell.debris_mass_kg / FRAGMENT_MASS_KG,
        )
        shell.debris_count -= consumed
        shell.debris_mass_kg -= consumed * FRAGMENT_MASS_KG
        remaining -= consumed

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
        return amount - remaining

    def _adaptive_substep_years(self, remaining_years: float) -> float:
        """Choose a step that resolves shell-to-shell drag and collision exposure."""
        from src.physics.collisions import calculate_impact_metrics
        from src.physics.drag import calculate_orbital_decay_rate

        max_substep_years = remaining_years
        max_collision_rate = 0.0

        for shell in self.shells:
            transiting_in_shell = self.get_transiting_debris_in_shell(
                shell.altitude_km
            )
            if shell.debris_count + transiting_in_shell > 0:
                fragment_rate = calculate_orbital_decay_rate(
                    shell.altitude_km, mass_kg=1.0, cross_section_m2=0.01
                )["decay_rate_km_year"]
                max_substep_years = min(
                    max_substep_years,
                    MAX_DRAG_DISPLACEMENT_KM / fragment_rate,
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
                shell.debris_count + transiting_in_shell
            )
            impact_metrics = calculate_impact_metrics(
                effective_debris / shell.volume_km3,
                5.0,
                remaining_years,
                10.0,
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
        from src.physics.collisions import calculate_impact_metrics
        from src.physics.drag import calculate_orbital_decay_rate
        from src.physics.breakup import calculate_breakup_fragments

        downward_debris_flux = 0.0
        downward_debris_mass_flux = 0.0
        downward_derelict_flux = 0.0
        downward_derelict_mass_flux = 0.0

        for i in reversed(range(len(self.shells))):
            shell = self.shells[i]
            transiting_in_shell = self.get_transiting_debris_in_shell(
                shell.altitude_km
            )
            effective_debris = shell.debris_count + transiting_in_shell
            spatial_density = effective_debris / shell.volume_km3
            impact_metrics = calculate_impact_metrics(
                spatial_density, 5.0, dt_years, 10.0
            )

            active_targets = shell.sat_count
            derelict_targets = shell.derelict_count
            collision_probability = -math.expm1(
                -impact_metrics["rate_per_year"] * dt_years
            )
            raw_active_collisions = collision_probability * active_targets
            raw_derelict_collisions = collision_probability * derelict_targets
            raw_total_collisions = raw_active_collisions + raw_derelict_collisions

            transit_projectiles = sum(
                min(
                    batch.count * batch.occupancy_fraction(shell.altitude_km),
                    batch.mass_kg / FRAGMENT_MASS_KG,
                )
                for batch in self.transiting_batches
            )
            projectile_capacity = min(
                shell.debris_count,
                shell.debris_mass_kg / FRAGMENT_MASS_KG,
            ) + transit_projectiles
            target_capacity = active_targets + derelict_targets
            total_shatter_events = min(
                raw_total_collisions,
                target_capacity,
                projectile_capacity,
            )

            if total_shatter_events > 0 and raw_total_collisions > 0:
                projectile_count_consumed = self._consume_projectiles(
                    shell, total_shatter_events
                )
                total_shatter_events = projectile_count_consumed
                active_events = total_shatter_events * (
                    raw_active_collisions / raw_total_collisions
                )
                derelict_events = total_shatter_events - active_events

                def apply_collision(
                    target_count: float, target_mass: float, events: float
                ):
                    if events <= 0 or target_count <= 0:
                        return 0.0, 0.0

                    average_target_mass = target_mass / target_count
                    breakup_data = calculate_breakup_fragments(
                        target_mass_kg=average_target_mass,
                        projectile_mass_kg=FRAGMENT_MASS_KG,
                        impact_velocity_km_s=10.0,
                        fragment_mass_kg=FRAGMENT_MASS_KG,
                    )
                    shell.debris_count += (
                        events * breakup_data["fragments_generated"]
                    )
                    shell.debris_mass_kg += (
                        events * breakup_data["ejecta_mass_kg"]
                    )

                    if breakup_data["is_catastrophic"]:
                        destroyed_count = min(events, target_count)
                        destroyed_mass = destroyed_count * average_target_mass
                        return destroyed_count, destroyed_mass

                    retained_mass_change = (
                        FRAGMENT_MASS_KG - breakup_data["ejecta_mass_kg"]
                    ) * events
                    return 0.0, -retained_mass_change

                active_destroyed, active_mass_removed = apply_collision(
                    active_targets, shell.satellite_mass_kg, active_events
                )
                derelict_destroyed, derelict_mass_removed = apply_collision(
                    derelict_targets, shell.derelict_mass_kg, derelict_events
                )

                if active_destroyed > 0:
                    legacy_ratio = shell.legacy_sat_count / active_targets
                    legacy_destroyed = active_destroyed * legacy_ratio
                    shell.legacy_sat_count -= legacy_destroyed
                    shell.new_sat_count -= active_destroyed - legacy_destroyed
                shell.satellite_mass_kg -= active_mass_removed
                shell.derelict_count -= derelict_destroyed
                shell.derelict_mass_kg -= derelict_mass_removed

            fragment_rate = calculate_orbital_decay_rate(
                shell.altitude_km, mass_kg=1.0, cross_section_m2=0.01
            )["decay_rate_km_year"]
            decay_frac_frag = -math.expm1(
                -(fragment_rate * dt_years) / 50.0
            )
            debris_leaving = shell.debris_count * decay_frac_frag
            debris_mass_leaving = shell.debris_mass_kg * decay_frac_frag
            shell.debris_count -= debris_leaving
            shell.debris_mass_kg -= debris_mass_leaving

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
                shell.derelict_count += downward_derelict_flux
                shell.derelict_mass_kg += downward_derelict_mass_flux

            downward_debris_flux = debris_leaving
            downward_debris_mass_flux = debris_mass_leaving
            downward_derelict_flux = derelicts_leaving
            downward_derelict_mass_flux = derelict_mass_leaving

            if i == 0:
                self.reentered_mass_kg += (
                    debris_mass_leaving + derelict_mass_leaving
                )

    def step(self, dt_years: float = 1.0):
        self._validate_step_inputs(dt_years)
        if dt_years == 0:
            return
        launch_count = self._project_launches(dt_years)

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
            loss_fraction = min(1.0, dt_years / self.sat_lifespan_years)
            legacy_loss = shell.legacy_sat_count * loss_fraction
            new_loss = shell.new_sat_count * loss_fraction
            total_satellite_loss = legacy_loss + new_loss
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

        # 5. Inject new launches
        target_shells = [s for s in self.shells if s.altitude_km in [450, 500, 550]]
        if target_shells:
            per_shell = launch_count / len(target_shells)
            for s in target_shells:
                s.new_sat_count += per_shell
                launched_mass = per_shell * SATELLITE_MASS_KG
                s.satellite_mass_kg += launched_mass
                self.launched_mass_kg += launched_mass

        self.current_year += dt_years
        for shell in self.shells:
            shell.update_mass_count_snapshot()