import math
import sys
from typing import Optional

MATERIAL_DENSITY_KG_M3 = 2700.0
SHRAPNEL_DIAMETER_M = 0.01
SHRAPNEL_MASS_KG = (
    MATERIAL_DENSITY_KG_M3 * (math.pi / 6.0) * SHRAPNEL_DIAMETER_M**3
)
SHRAPNEL_CROSS_SECTION_M2 = math.pi * (SHRAPNEL_DIAMETER_M / 2.0) ** 2
HEAVY_CHUNK_MASS_KG = 10.0
HEAVY_CHUNK_RADIUS_M = (
    3.0 * HEAVY_CHUNK_MASS_KG
    / (4.0 * math.pi * MATERIAL_DENSITY_KG_M3)
) ** (1.0 / 3.0)
HEAVY_CHUNK_CROSS_SECTION_M2 = math.pi * HEAVY_CHUNK_RADIUS_M**2
HEAVY_MASS_FRACTION = 0.8
WHIPPLE_SHIELD_LETHAL_ENERGY_J = 100_000.0


def calculate_kinetic_energy_j(
    projectile_mass_kg: float, impact_velocity_km_s: float
) -> float:
    if not math.isfinite(projectile_mass_kg) or projectile_mass_kg < 0:
        raise ValueError("projectile_mass_kg must be finite and nonnegative.")
    if not math.isfinite(impact_velocity_km_s) or impact_velocity_km_s < 0:
        raise ValueError("impact_velocity_km_s must be finite and nonnegative.")
    if projectile_mass_kg == 0.0 or impact_velocity_km_s == 0.0:
        return 0.0

    log_kinetic_energy = (
        math.log(0.5)
        + math.log(projectile_mass_kg)
        + 2.0 * (math.log(impact_velocity_km_s) + math.log(1000.0))
    )
    if log_kinetic_energy > math.log(sys.float_info.max):
        raise ValueError("Impact kinetic energy must be finite.")
    kinetic_energy_j = math.exp(log_kinetic_energy)
    if not math.isfinite(kinetic_energy_j):
        raise ValueError("Impact kinetic energy must be finite.")
    return kinetic_energy_j


def calculate_breakup_fragments(
    target_mass_kg: float,
    projectile_mass_kg: float,
    impact_velocity_km_s: float,
    min_fragment_size_m: float = SHRAPNEL_DIAMETER_M,
    fragment_mass_kg: Optional[float] = None,
    force_catastrophic: bool = False,
) -> dict:
    """
    Splits collision ejecta between 1 cm shrapnel and 10 kg heavy chunks.
    
    :param target_mass_kg: Mass of the satellite being hit.
    :param projectile_mass_kg: Mass of the incoming debris.
    :param impact_velocity_km_s: Relative velocity of the impact.
    :param min_fragment_size_m: Minimum characteristic length of shrapnel to track.
    :param fragment_mass_kg: Mass assigned to each tracked shrapnel piece.
    """
    inputs = {
        "target_mass_kg": target_mass_kg,
        "projectile_mass_kg": projectile_mass_kg,
        "impact_velocity_km_s": impact_velocity_km_s,
        "min_fragment_size_m": min_fragment_size_m,
        "fragment_mass_kg": fragment_mass_kg
    }
    if not all(
        math.isfinite(value)
        for value in inputs.values()
        if value is not None
    ):
        raise ValueError("All inputs must be finite numbers.")
    if target_mass_kg <= 0:
        raise ValueError("target_mass_kg must be greater than zero.")
    if projectile_mass_kg < 0:
        raise ValueError("projectile_mass_kg must not be negative.")
    if impact_velocity_km_s < 0:
        raise ValueError("impact_velocity_km_s must not be negative.")
    if min_fragment_size_m <= 0:
        raise ValueError("min_fragment_size_m must be greater than zero.")
    if fragment_mass_kg is not None and fragment_mass_kg <= 0:
        raise ValueError("fragment_mass_kg must be greater than zero.")

    shrapnel_unit_mass_kg = (
        fragment_mass_kg
        if fragment_mass_kg is not None
        else MATERIAL_DENSITY_KG_M3 * (math.pi / 6.0) * min_fragment_size_m**3
    )
    
    kinetic_energy_j = calculate_kinetic_energy_j(
        projectile_mass_kg, impact_velocity_km_s
    )
    
    # 2. Specific energy (Joules per kg of the target)
    specific_energy_j_kg = kinetic_energy_j / target_mass_kg
    if not math.isfinite(specific_energy_j_kg):
        raise ValueError("Specific impact energy must be finite.")
    
    # 3. Determine if the collision is catastrophic (Threshold: 40 J/g or 40,000 J/kg)
    catastrophic_threshold = 40000.0
    is_catastrophic = (
        force_catastrophic or specific_energy_j_kg > catastrophic_threshold
    )
    
    # 4. Determine total mass involved in fragment generation
    if is_catastrophic:
        # The entire mass of both objects shatters
        ejecta_mass_kg = target_mass_kg + projectile_mass_kg
    else:
        # Cratering: Only a localized chunk is blown off (NASA approximation)
        ejecta_mass_kg = projectile_mass_kg * impact_velocity_km_s

    ejecta_mass_kg = min(
        max(ejecta_mass_kg, 0.0),
        target_mass_kg + projectile_mass_kg,
    )
    if not math.isfinite(ejecta_mass_kg):
        raise ValueError("Ejecta mass must be finite.")
    
    # Reserve most ejecta for intact structural pieces, then distribute the
    # remainder into shrapnel at the mass of a 1 cm aluminum sphere.
    heavy_chunk_count = math.floor(
        ejecta_mass_kg * HEAVY_MASS_FRACTION / HEAVY_CHUNK_MASS_KG
    )
    heavy_chunk_mass = heavy_chunk_count * HEAVY_CHUNK_MASS_KG
    remaining_mass = ejecta_mass_kg - heavy_chunk_mass
    shrapnel_count = remaining_mass / shrapnel_unit_mass_kg
    if not math.isfinite(shrapnel_count):
        raise ValueError("Fragment count exceeds the supported numeric range.")
    shrapnel_mass = remaining_mass
    untracked_mass_kg = max(0.0, ejecta_mass_kg - heavy_chunk_mass - shrapnel_mass)
    
    return {
        "is_catastrophic": is_catastrophic,
        "specific_energy_j_kg": specific_energy_j_kg,
        "ejecta_mass_kg": ejecta_mass_kg,
        "fragments_generated": int(shrapnel_count),
        "shrapnel_count": shrapnel_count,
        "shrapnel_unit_mass_kg": shrapnel_unit_mass_kg,
        "shrapnel_mass_kg": shrapnel_mass,
        "heavy_chunk_count": int(heavy_chunk_count),
        "heavy_chunk_mass_kg": heavy_chunk_mass,
        "untracked_mass_kg": untracked_mass_kg,
    }

# --- Quick Validation Test ---
if __name__ == "__main__":
    # Example 1: Iridium-Cosmos style catastrophic collision
    # 500kg target hit by a 1kg piece of debris at 10 km/s
    result_catastrophic = calculate_breakup_fragments(
        target_mass_kg=500.0,
        projectile_mass_kg=1.0,
        impact_velocity_km_s=10.0,
    )
    
    print("--- CATASTROPHIC IMPACT (1kg at 10km/s vs 500kg) ---")
    print(f"Catastrophic? {result_catastrophic['is_catastrophic']}")
    print(f"Shrapnel pieces (>1cm): {result_catastrophic['shrapnel_count']:,}")
    print(f"Heavy chunks: {result_catastrophic['heavy_chunk_count']:,}")
    
    # Example 2: Non-catastrophic cratering
    # 500kg target hit by a tiny 1-gram (0.001kg) fleck of paint at 10 km/s
    result_cratering = calculate_breakup_fragments(
        target_mass_kg=500.0,
        projectile_mass_kg=0.001,
        impact_velocity_km_s=10.0
    )
    
    print("\n--- CRATERING IMPACT (1g at 10km/s vs 500kg) ---")
    print(f"Catastrophic? {result_cratering['is_catastrophic']}")
    print(f"Shrapnel pieces (>1cm): {result_cratering['shrapnel_count']:,}")
    print(f"Heavy chunks: {result_cratering['heavy_chunk_count']:,}")