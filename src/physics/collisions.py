import math

def calculate_impact_metrics(
    spatial_density_km3: float,
    cross_section_m2: float,
    time_horizon_years: float,
    relative_velocity_km_s: float = 10.0
) -> dict:
    """
    Calculates the risk of a satellite being hit by debris in a specific orbital shell.
    
    :param spatial_density_km3: Number of debris objects per cubic kilometer in this shell.
    :param cross_section_m2: The frontal area of the satellite in square meters.
    :param time_horizon_years: How long the satellite plans to operate in this shell.
    :param relative_velocity_km_s: Average impact speed (LEO average is ~10 km/s).
    """
    # 1. Convert cross-section to square kilometers to match density units
    cross_section_km2 = cross_section_m2 / 1_000_000.0
    
    # 2. Calculate volumetric sweep rate (Collisions per second)
    # Rate = Density * Area * Velocity
    rate_per_second = spatial_density_km3 * cross_section_km2 * relative_velocity_km_s
    
    # 3. Convert to annualized metrics
    seconds_per_year = 365.25 * 24 * 3600
    rate_per_year = rate_per_second * seconds_per_year
    
    # 4. Calculate Mean Time Between Impacts (MTBI)
    mtbi_years = 1.0 / rate_per_year if rate_per_year > 0 else float('inf')
    
    # 5. Poisson Probability of AT LEAST ONE collision over the time horizon
    # P(>0) = 1 - e^(-rate * time)
    probability_of_impact = 1.0 - math.exp(-rate_per_year * time_horizon_years)
    
    return {
        "rate_per_year": rate_per_year,
        "mtbi_years": mtbi_years,
        "probability_of_impact": probability_of_impact
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