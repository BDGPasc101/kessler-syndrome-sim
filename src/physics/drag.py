import math

def calculate_atmospheric_density(altitude_km: float) -> float:
    """
    Approximates Earth's atmospheric density using a simple exponential model.
    Values are calibrated for a mean solar cycle in Low Earth Orbit (LEO).
    """
    if altitude_km > 1200:
        return 0.0  # Absolute vacuum for our simulation purposes
        
    # Reference values for 600 km altitude
    h0 = 600.0
    rho0 = 1.37e-13  # Reference density in kg/m^3
    H = 50.0         # Scale height in km (how fast density drops off)
    
    return rho0 * math.exp(-(altitude_km - h0) / H)

def calculate_orbital_decay_rate(
    altitude_km: float,
    mass_kg: float,
    cross_section_m2: float,
    drag_coefficient: float = 2.2
) -> dict:
    """
    Calculates how fast an object is losing altitude due to atmospheric drag.
    
    :param altitude_km: Current orbital altitude.
    :param mass_kg: Mass of the debris/satellite.
    :param cross_section_m2: Frontal area facing the "wind".
    :param drag_coefficient: Aerodynamic shape factor (2.2 is standard for satellites).
    """
    density_kg_m3 = calculate_atmospheric_density(altitude_km)
    
    # 1. Ballistic Coefficient (kg/m^2) = M / (Cd * A)
    # A higher BC means the object punches through the air easier (e.g., a bowling ball).
    # A lower BC means it gets dragged down faster (e.g., a feather).
    if cross_section_m2 == 0 or drag_coefficient == 0:
        ballistic_coefficient = float('inf')
    else:
        ballistic_coefficient = mass_kg / (drag_coefficient * cross_section_m2)
        
    # 2. Orbital mechanics constants
    earth_mass_kg = 5.972e24
    G = 6.67430e-11
    earth_radius_m = 6371000.0
    
    orbital_radius_m = earth_radius_m + (altitude_km * 1000.0)
    
    # Standard gravitational parameter (mu)
    mu = G * earth_mass_kg
    
    # 3. Rate of altitude loss (dr/dt)
    # Derived from orbital energy dissipation: dr/dt = - (rho * sqrt(mu * r)) / BC
    if ballistic_coefficient == float('inf') or density_kg_m3 == 0:
        decay_rate_m_s = 0.0
    else:
        decay_rate_m_s = (density_kg_m3 * math.sqrt(mu * orbital_radius_m)) / ballistic_coefficient
        
    # Convert meters per second to kilometers per year
    seconds_per_year = 365.25 * 24 * 3600.0
    decay_rate_km_year = (decay_rate_m_s / 1000.0) * seconds_per_year
    
    return {
        "density_kg_m3": density_kg_m3,
        "ballistic_coefficient": ballistic_coefficient,
        "decay_rate_km_year": decay_rate_km_year
    }

# --- Quick Validation Test ---
if __name__ == "__main__":
    # Example 1: The ISS (Huge area, heavy, at 400km)
    iss_decay = calculate_orbital_decay_rate(
        altitude_km=400.0,
        mass_kg=420_000.0,
        cross_section_m2=2500.0
    )
    print("--- INTERNATIONAL SPACE STATION (400km) ---")
    print(f"Atmospheric Density: {iss_decay['density_kg_m3']:.2e} kg/m^3")
    print(f"Decay Rate: {iss_decay['decay_rate_km_year']:.1f} km per year")
    
    # Example 2: A 10cm chunk of shattered debris at 800km (Kessler danger zone)
    debris_decay = calculate_orbital_decay_rate(
        altitude_km=800.0,
        mass_kg=1.0,         # 1 kg aluminum chunk
        cross_section_m2=0.01 # 10x10cm area
    )
    print("\n--- DEBRIS CHUNK (800km) ---")
    print(f"Atmospheric Density: {debris_decay['density_kg_m3']:.2e} kg/m^3")
    print(f"Decay Rate: {debris_decay['decay_rate_km_year']:.4f} km per year")