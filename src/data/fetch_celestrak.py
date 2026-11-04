import requests
import math
import json
import os

# Orbital Mechanics Constants
MU_EARTH = 398600.4418  # Earth's standard gravitational parameter (km^3/s^2)
EARTH_RADIUS = 6371.0   # Earth's mean radius in km

def calculate_altitude(mean_motion_rev_per_day: float) -> float:
    """
    Converts Mean Motion (Revolutions per Day) into Mean Altitude (km)
    using Kepler's Third Law of Planetary Motion.
    """
    if mean_motion_rev_per_day <= 0:
        return 0.0
        
    mean_motion_rad_s = mean_motion_rev_per_day * (2.0 * math.pi / 86400.0)
    semi_major_axis_km = (MU_EARTH / (mean_motion_rad_s ** 2)) ** (1.0 / 3.0)
    altitude_km = semi_major_axis_km - EARTH_RADIUS
    return altitude_km

def fetch_and_process():
    # CelesTrak General Perturbation (GP) JSON endpoints
    endpoints = {
        "active_satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=json",
        "fengyun_1c_debris": "https://celestrak.org/NORAD/elements/gp.php?GROUP=fengyun-1c-debris&FORMAT=json",
        "iridium_33_debris": "https://celestrak.org/NORAD/elements/gp.php?GROUP=iridium-33-debris&FORMAT=json",
        "cosmos_2251_debris": "https://celestrak.org/NORAD/elements/gp.php?GROUP=cosmos-2251-debris&FORMAT=json"
    }

    # Prepare our 50km altitude buckets (400km to 1000km)
    shells = {alt: {"satellites": 0, "debris": 0} for alt in range(400, 1050, 50)}
    
    for category, url in endpoints.items():
        print(f"Fetching {category}...")
        response = requests.get(url)
        
        if response.status_code != 200:
            print(f"  -> Failed: HTTP {response.status_code}")
            continue
            
        # SAFETY NET: Gracefully handle non-JSON responses
        try:
            data = response.json()
        except requests.exceptions.JSONDecodeError:
            print(f"  -> Failed: CelesTrak returned invalid JSON for {category}. Skipping.")
            continue
        
        for obj in data:
            mean_motion = obj.get("MEAN_MOTION")
            obj_type = obj.get("OBJECT_TYPE", "UNKNOWN")
            
            if not mean_motion:
                continue
                
            alt = calculate_altitude(mean_motion)
            
            if 400 <= alt < 1050:
                bucket = int(alt // 50) * 50
                if obj_type == "PAYLOAD":
                    shells[bucket]["satellites"] += 1
                else:
                    shells[bucket]["debris"] += 1

    # Apply the 1cm dark debris empirical multiplier (50x)
    for alt in shells:
        shells[alt]["debris"] *= 50

    os.makedirs("data", exist_ok=True)
    output_path = "data/celestrak_shells.json"
    with open(output_path, "w") as f:
        json.dump(shells, f, indent=4)
        
    print(f"\nSuccess! Categorized data saved to {output_path}")
    print("\nSample Data (400km - 1000km grid):")
    for alt, counts in shells.items():
        print(f"{alt}km Shell -> Sats: {counts['satellites']:<4} | Untracked Debris Estimate (>1cm): {counts['debris']:,}")

if __name__ == "__main__":
    fetch_and_process()