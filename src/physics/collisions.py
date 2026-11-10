import math
from typing import Optional

EARTH_MU_KM3_S2 = 3.986e5
EARTH_RADIUS_KM = 6371.0
AVERAGE_COLLISION_ANGLE_DEG = 90.0


def calculate_relative_velocity_km_s(altitude_km: float) -> float:
    """Estimate relative speed for circular orbits at a shared altitude."""
    if not math.isfinite(altitude_km) or altitude_km <= -EARTH_RADIUS_KM:
        raise ValueError("altitude_km must be finite and above Earth's center.")

    orbital_radius_km = EARTH_RADIUS_KM + altitude_km
    orbital_velocity_km_s = math.sqrt(EARTH_MU_KM3_S2 / orbital_radius_km)
    collision_angle_rad = math.radians(AVERAGE_COLLISION_ANGLE_DEG)
    relative_velocity_squared = 2.0 * orbital_velocity_km_s**2 * (
        1.0 - math.cos(collision_angle_rad)
    )
    return math.sqrt(relative_velocity_squared)


def calculate_impact_metrics(
    spatial_density_km3: float,
    cross_section_m2: float,
    time_horizon_years: float,
    relative_velocity_km_s: float = 10.0,
    lethal_spatial_density_km3: Optional[float] = None,
) -> dict:
    """
    Calculates the risk of a satellite being hit by debris in a specific orbital shell.
    
    :param spatial_density_km3: Number of debris objects per cubic kilometer in this shell.
    :param cross_section_m2: The frontal area of the satellite in square meters.
    :param time_horizon_years: How long the satellite plans to operate in this shell.
    :param relative_velocity_km_s: Average impact speed in km/s.
    :param lethal_spatial_density_km3: Density of projectiles treated as lethal.
        Defaults to spatial_density_km3 for backward compatibility.
    """
    if lethal_spatial_density_km3 is None:
        lethal_spatial_density_km3 = spatial_density_km3
    inputs = (
        spatial_density_km3,
        lethal_spatial_density_km3,
        cross_section_m2,
        time_horizon_years,
        relative_velocity_km_s,
    )
    if not all(math.isfinite(value) for value in inputs):
        raise ValueError("Impact metric inputs must be finite.")
    if any(value < 0 for value in inputs):
        raise ValueError("Impact metric inputs must be nonnegative.")

    cross_section_km2 = cross_section_m2 / 1_000_000.0
    seconds_per_year = 365.25 * 24 * 3600

    rate_per_year = (
        spatial_density_km3
        * cross_section_km2
        * relative_velocity_km_s
        * seconds_per_year
    )
    lethal_rate_per_year = (
        lethal_spatial_density_km3
        * cross_section_km2
        * relative_velocity_km_s
        * seconds_per_year
    )
    mtbi_years = 1.0 / rate_per_year if rate_per_year > 0 else float('inf')
    probability_of_impact = 1.0 - math.exp(-rate_per_year * time_horizon_years)
    probability_of_catastrophic_breach = (
        1.0 - math.exp(-lethal_rate_per_year * time_horizon_years)
    )
    
    return {
        "rate_per_year": rate_per_year,
        "mtbi_years": mtbi_years,
        "probability_of_impact": probability_of_impact,
        "lethal_rate_per_year": lethal_rate_per_year,
        "probability_of_catastrophic_breach": probability_of_catastrophic_breach,
    }

# --- Quick Validation Test ---
if __name__ == "__main__":
    # Example: A large satellite (10m^2) in a dense LEO shell (800km altitude)
    # Dense LEO typically has around 1e-8 objects per km^3 > 1cm in size.
    test_metrics = calculate_impact_metrics(
        spatial_density_km3=1e-8,
        cross_section_m2=10.0,
        time_horizon_years=5.0
    )
    
    print(f"Impacts per year: {test_metrics['rate_per_year']:.4f}")
    print(f"Mean Time Between Impacts: {test_metrics['mtbi_years']:.1f} years")
    print(f"Probability of impact in 5 years: {test_metrics['probability_of_impact'] * 100:.2f}%")