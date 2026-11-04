import math
import sys
import os
import json

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.physics.collisions import calculate_impact_metrics
from src.physics.breakup import calculate_breakup_fragments
from src.physics.drag import calculate_orbital_decay_rate

class OrbitalShell:
    def __init__(self, altitude_km: float, thickness_km: float, initial_sats: int, initial_debris: int):
        self.altitude_km = altitude_km
        self.thickness_km = thickness_km
        # Track legacy satellites (already in orbit at year 0) separately from newly launched ones
        self.legacy_sat_count = float(initial_sats) 
        self.new_sat_count = 0.0
        self.debris_count = float(initial_debris)
        
        earth_r = 6371.0
        r_inner = earth_r + altitude_km - (thickness_km / 2.0)
        r_outer = earth_r + altitude_km + (thickness_km / 2.0)
        self.volume_km3 = (4.0 / 3.0) * math.pi * (r_outer**3 - r_inner**3)
        
    @property
    def sat_count(self):
        """Total active satellites in this shell."""
        return self.legacy_sat_count + self.new_sat_count

class KesslerSimulation:
    def __init__(self):
        self.SAT_MASS = 500.0       
        self.SAT_AREA = 5.0         
        self.DEBRIS_MASS = 1.0      
        self.DEBRIS_AREA = 0.01     
        self.IMPACT_VEL = 10.0      
        
        self.current_year = 0.0
        
        self.base_annual_launches = 2500  
        self.launch_growth_rate = 0.08    
        self.sat_lifespan_years = 5.0     
        self.legacy_sat_lifespan_years = 5.0 # Assuming existing satellites also have a 5-year remaining life on average
        
        data_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data', 'celestrak_shells.json'))
        if not os.path.exists(data_path):
            raise FileNotFoundError(f"Data file not found at {data_path}. Run fetch_celestrak.py first!")
            
        with open(data_path, 'r') as f:
            real_data = json.load(f)
            
        self.shells = []
        for alt in range(400, 1050, 50):
            str_alt = str(alt)
            starting_sats = real_data.get(str_alt, {}).get("satellites", 0)
            starting_debris = real_data.get(str_alt, {}).get("debris", 0)
            
            self.shells.append(OrbitalShell(
                float(alt), 50.0, starting_sats, starting_debris
            ))

    def get_annual_deployments(self) -> list:
        total_launches = self.base_annual_launches * ((1.0 + self.launch_growth_rate) ** self.current_year)
        shell_deployments = [0.0] * len(self.shells)
        
        for i, shell in enumerate(self.shells):
            if shell.altitude_km == 550:
                shell_deployments[i] = total_launches * 0.70
            elif shell.altitude_km == 600:
                shell_deployments[i] = total_launches * 0.15
            else:
                shell_deployments[i] = total_launches * (0.15 / (len(self.shells) - 2))
                
        return shell_deployments

    def step(self, dt_years: float):
        new_debris = [0.0] * len(self.shells)
        falling_debris = [0.0] * len(self.shells)
        
        deployments = self.get_annual_deployments()
        
        for i, shell in enumerate(self.shells):
            # 1. Update Intact Satellite Population
            shell.new_sat_count += deployments[i] * dt_years
            
            # Active de-orbit for NEW satellites
            new_retirements = shell.new_sat_count * (1.0 / self.sat_lifespan_years) * dt_years
            shell.new_sat_count = max(0.0, shell.new_sat_count - new_retirements)
            
            # Active de-orbit for LEGACY satellites
            # They decay over their estimated remaining lifespan
            if shell.legacy_sat_count > 0:
                legacy_retirements = shell.legacy_sat_count * (1.0 / self.legacy_sat_lifespan_years) * dt_years
                shell.legacy_sat_count = max(0.0, shell.legacy_sat_count - legacy_retirements)

            if shell.debris_count <= 0:
                continue
                
            spatial_density = shell.debris_count / shell.volume_km3
            
            # 2. Collisions (using total sat_count property)
            sat_metrics = calculate_impact_metrics(spatial_density, self.SAT_AREA, dt_years, self.IMPACT_VEL)
            sat_hits = shell.sat_count * sat_metrics['rate_per_year'] * dt_years
            
            actual_sat_hits = min(shell.sat_count, sat_hits)
            
            deb_metrics = calculate_impact_metrics(spatial_density, self.DEBRIS_AREA, dt_years, self.IMPACT_VEL)
            deb_hits = 0.5 * shell.debris_count * deb_metrics['rate_per_year'] * dt_years
            
            # 3. Process Satellite Destruction
            if actual_sat_hits > 0:
                sat_frag = calculate_breakup_fragments(self.SAT_MASS, self.DEBRIS_MASS, self.IMPACT_VEL)['fragments_generated']
                new_debris[i] += actual_sat_hits * sat_frag
                
                # Deduct destroyed satellites proportionally from new and legacy pools
                if shell.sat_count > 0:
                    legacy_ratio = shell.legacy_sat_count / shell.sat_count
                    new_ratio = shell.new_sat_count / shell.sat_count
                    
                    shell.legacy_sat_count = max(0.0, shell.legacy_sat_count - (actual_sat_hits * legacy_ratio))
                    shell.new_sat_count = max(0.0, shell.new_sat_count - (actual_sat_hits * new_ratio))
                
            # 4. Process Debris-on-Debris Grinding
            if deb_hits > 0:
                deb_frag = calculate_breakup_fragments(self.DEBRIS_MASS, self.DEBRIS_MASS, self.IMPACT_VEL)['fragments_generated']
                new_debris[i] += deb_hits * deb_frag
                
            # 5. Atmospheric Drag
            decay = calculate_orbital_decay_rate(shell.altitude_km, self.DEBRIS_MASS, self.DEBRIS_AREA)['decay_rate_km_year']
            fraction_lost = min(1.0, (decay * dt_years) / shell.thickness_km)
            falling_debris[i] = shell.debris_count * fraction_lost
            
        # Apply mass transfers
        for i in range(len(self.shells)):
            self.shells[i].debris_count += new_debris[i] - falling_debris[i]
            if i > 0:
                self.shells[i-1].debris_count += falling_debris[i]
                
        # 6. Active Debris Removal (Laser Broom)
        if hasattr(self, 'laser_removals_per_year') and self.laser_removals_per_year > 0:
            removals_remaining = self.laser_removals_per_year * dt_years
            congested_shells = sorted(self.shells, key=lambda s: s.debris_count, reverse=True)
            for shell in congested_shells:
                if removals_remaining <= 0:
                    break
                if shell.debris_count > 0:
                    removed = min(shell.debris_count, removals_remaining)
                    shell.debris_count -= removed
                    removals_remaining -= removed
                
        self.current_year += dt_years