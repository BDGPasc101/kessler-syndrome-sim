# Orbital Debris Cascade Simulator: Engineering & Physics Learnings

## 1. Physics Engine Architecture

The simulation models Low Earth Orbit (LEO) as a discrete set of 13 concentric shells spanning altitudes from **400 km to 1,000 km** in **50 km increments**. Dynamics are evaluated as coupled ordinary differential equations (ODEs) updated via forward numerical integration ($\Delta t = 1.0\text{ year}$).

### 1.1 Shell Geometry & Spatial Density
Each shell is defined as a spherical shell around Earth:

$$V_{\text{shell}} = \frac{4}{3} \pi \left[ (R_E + h + \Delta h / 2)^3 - (R_E + h - \Delta h / 2)^3 \right]$$

* $R_E = 6,378.137\text{ km}$ (WGS-84 equatorial reference radius)
* $h$: Shell center altitude ($400, 450, \dots, 1000\text{ km}$)
* $\Delta h = 50\text{ km}$ (Shell thickness)

Spatial density is computed dynamically at every time step:
$$\rho = \frac{N_{\text{debris}}}{V_{\text{shell}}} \quad \left[\text{objects} / \text{km}^3\right]$$

### 1.2 Collision Mechanics (Kinetic Gas Theory)
Collisions follow kinetic gas collision theory adapted for orbital cross-sections:
$$\text{Flux} = \rho \cdot v_{\text{rel}}$$
$$\text{Collision Rate} = \rho \cdot \sigma \cdot v_{\text{rel}} \cdot \Delta t$$

* $\sigma$: Combined collision cross-section ($\text{m}^2$)
* $v_{\text{rel}}$: Mean relative velocity in LEO, calibrated to $10.0\text{ km/s}$

The probability of impact over a mission lifetime $T$ uses Poisson statistics:
$$P(\text{impact}) = 1 - e^{-\text{Rate} \cdot T}$$
$$\text{Mean Time Between Impacts (MTBI)} = \frac{1}{\text{Rate}}$$

### 1.3 Fragmentation: NASA Standard Breakup Model
When collisions occur, fragment generation depends on whether the specific kinetic energy surpasses the catastrophic threshold ($E_p / M_t \ge 40\text{ kJ/kg}$):

* **Catastrophic Breakup ($E_p / M_t \ge 40\text{ kJ/kg}$):** Complete destruction of both target and projectile. Total mass involved is $M = M_t + M_p$.
* **Non-Catastrophic / Crater Impact ($E_p / M_t < 40\text{ kJ/kg}$):** The target remains intact; only projectile mass ejects local ejecta fragments.

For fragments larger than characteristic length $L_c \ge 0.1\text{ m}$ (and scaling down to $1\text{ cm}$):
$$N(L > L_c) = 0.1 \cdot (M)^{0.75} \cdot L_c^{-1.71}$$

### 1.4 Atmospheric Drag & Downward Flux
Debris loses orbital energy through thermospheric drag. The circular orbital decay rate is parameterized by:

$$\frac{da}{dt} = -2\pi \left(\frac{C_D \cdot A}{m}\right) \rho_{\text{atm}}(h) \cdot a^2$$

To ensure computational performance across large sweeps without integrating every trajectory individually, the downward transfer between shells is modeled as a mass flux:
$$\text{Fraction Lost} = \min\left(1.0, \frac{\dot{a}_{\text{decay}} \cdot \Delta t}{\Delta h}\right)$$
Mass purged from shell $i$ cascades down into shell $i-1$. Debris leaving the $400\text{ km}$ shell burns up in the dense lower mesosphere and is permanently excised from the system.

### 1.5 Active Debris Removal (Laser Broom)
Ground- or space-based pulsed laser ablation imparts a small retrograde velocity vector ($\Delta v$) onto debris perimeters, lowering their perigee into dense atmosphere for accelerated burn-up.
* The system prioritizes the shell with the highest absolute debris count.
* Removals are bounded by the facility's annual firing capacity ($N_{\text{shots}} \cdot \Delta t$).

---

## 2. Key Engineering Decisions & Edge Cases

### 2.1 The "Legacy Satellite" Pool vs. New Launches
* **Problem:** In early iterations, setting `base_annual_launches = 0` caused the entire active satellite population to immediately plummet to zero, despite ~12,000+ active satellites orbiting Earth in reality.
* **Resolution:** Separated the shell satellite population into:
  * `legacy_sat_count`: Seeded directly from current live CelesTrak ephemeris data.
  * `new_sat_count`: Injected yearly based on deployment rates and constellation lifespan.
* Both decay over their active lifespans, but legacy satellites naturally taper off across the first decade without needing artificial replenishment.

### 2.2 Keplerian Altitude Binning & Constellation Geography
Converting Mean Motion ($n$, revs/day) to geometric altitude above WGS-84 revealed why ~9,800 satellites cluster into the $450\text{ km}$ and $500\text{ km}$ shells:
* $a = \left(\frac{\mu}{n^2}\right)^{1/3}$
* $h = a - 6,378.137\text{ km}$
* Starlink V2 Mini / Direct-to-Cell fleets operate across $480\text{ km}$–$530\text{ km}$ rather than older $550\text{ km}$ shells. With $50\text{ km}$ bin widths centered at $450\text{ km}$ (`[425, 475)`) and $500\text{ km}$ (`[475, 525)`), satellites stationed at $521\text{ km}$ fall squarely into the $500\text{ km}$ bin.

### 2.3 CelesTrak Ingestion: OMM JSON vs. CSV Streams
* **Bug Encountered:** Large JSON queries (`GROUP=active` or `GROUP=starlink`) triggered `HTTP 403 Forbidden` or `HTTP 500 Internal Server Error` due to back-end PHP memory exhaustion when generating dynamic JSON for 11,000+ records.
* **Fix:** Switched from `FORMAT=json` to `FORMAT=csv`. CSV payloads are a fraction of the size, parse sequentially via `csv.DictReader`, and bypass server-side JSON serialization limits.
* **Resilience:** Implemented an empirical baseline fallback so the simulation remains functional offline or during CelesTrak downtime.

### 2.4 Economic Modeling: Upfront vs. Recurring Costs
* Fixed capital expenditure (CAPEX) for laser facility construction is modeled as an upfront cost ($1B at $t = 0$).
* Annual operating expenditure (OPEX) scales with fixed operational upkeep ($50M/yr) plus variable costs per shot ($2,500/shot).
* Lifetime program costs ($T = 100\text{ years}$) reach $\approx \$9.75\text{ Billion}$.

---

## 3. Test Suite Documentation (`tests/test_simulation.py`)

A five-point validation suite guarantees physical and data integrity:

| Test Name | Primary Validation Target | Key Physics / Logic Verified |
|---|---|---|
| `test_initial_sat_count_nonzero_and_matches_data` | Data Ingestion | Confirms `celestrak_shells.json` exists, contains non-zero active satellite data, and initializes the simulation shells correctly. |
| `test_atmospheric_drag_monotonic_downward_flux` | Mass Conservation & Decay | Verifies that under zero launches and no collisions, debris at $800\text{ km}$ drains downwards into $750\text{ km}$ monotonically via atmospheric drag. |
| `test_nasa_breakup_energy_and_fragment_scaling` | Breakup Dynamics | Ensures catastrophic collision fragment counts scale with projectile and target mass according to NASA power-law distributions. |
| `test_collision_probability_density_scaling` | Kinetic Gas Theory | Asserts that collision rates scale linearly when spatial density doubles ($\text{Rate} \propto \rho$). |
| `test_laser_broom_targets_most_congested_shells_first` | Remediation Routing | Validates that active laser removal prioritizes the densest shell and does not exceed configured engagement quotas. |


My Questions:
"In what way is this model too simplistic - what assumptions do we make? For example, is the laser model too simple as it assumes every laser hit will cause a de-orbit? The material may remain an active hazard even if its velocity is reduce, as it drops to lower orbital shells"

R: My point "highlights the central compromise of this simulation: it is a 1D thermodynamic "particle-in-a-box" model, not an orbital mechanics propagator."

Steps to move it closer to an orbital mechanics propagator

## 4. Model Limitations & Physical Simplifications

While the 1D discrete-shell kinetic engine captures macro-scale cascade dynamics, it relies on several simplifying assumptions that abstract away real-world orbital and operational complexities:

1. **Laser Broom De-Orbit Dynamics & Transiting Hazards:**
   * **Instant Annihilation vs. Incremental $\Delta v$:** The model treats laser ablation as single-shot removals rather than thousands of low-impulse burns required to lower perigee.
   * **Crossing Hazards:** Retarding an 800 km orbit transforms it into an eccentric transfer ellipse that slices through lower operational shells (400–750 km). In reality, de-orbiting debris temporarily increases lower-shell collision risk before atmospheric burn-up.
   * **Secondary Fragmentation:** Real laser ablation carries risk of thermal shock or irregular torque, which can shatter composite debris into un-trackable millimeter fragments.

2. **Isotropic Spatial Distribution vs. Orbital Plane Geometry:**
   * **Homogeneous Gas Model:** The model distributes objects uniformly within each 50 km spherical shell.
   * **Pinch Points & Convergence:** Real orbits concentrate around specific inclination bands (e.g., Sun-Synchronous Orbits at $97^\circ\text{–}98^\circ$ and polar paths), creating high-density congestion zones over Earth's poles rather than uniform global dispersion.
   * **Neglected Orbital Eccentricity:** Rocket upper stages and transfer orbits (e.g., GTO) regularly sweep across multiple concentric shells at high relative velocities.

3. **In-Shell Fragment Confinement & Uniform Area-to-Mass ($A/m$):**
   * **Single-Shell Deposition:** All fragments generated in a collision are assigned back to the impact altitude. In reality, hypervelocity ejecta velocity vectors ($\Delta v_{\text{ejecta}}$) fling fragments hundreds of kilometers above and below the point of impact.
   * **Monolithic Ballistic Properties:** Debris is assigned fixed mass (1 kg) and area ($0.01\text{ m}^2$). Real debris follows a broad log-normal $A/m$ distribution, where high-drag foils decay in weeks while dense metallic fragments linger for centuries.

4. **Static Atmosphere vs. the 11-Year Solar Cycle:**
   * **Constant Thermospheric Density:** The atmospheric drag module uses a static density-altitude curve.
   * **Solar Flux Modulation:** Real thermospheric density fluctuates by a factor of 10 to 100 across the 11-year solar cycle ($F_{10.7}$ index). Solar maximums periodically scrub low-altitude debris, while solar minimums allow debris populations to stagnate.

5. **Operational Behavior & Numerical Discretization:**
   * **Zero Collision Avoidance (CAMs):** Active satellites are modeled as passive targets, ignoring real-world propulsive avoidance maneuvers triggered by Conjunction Assessment Risk Analysis (CARA).
   * **100% Post-Mission Disposal (PMD):** Retired satellites are assumed to de-orbit reliably at end-of-life, whereas actual industry PMD compliance hovers between 60% and 85%, leaving a steady fraction of derelict hulks in orbit.
   * **Time-Step Smoothing ($\Delta t = 1.0\text{ year}$):** Annual forward numerical integration averages out rapid, short-term collision cascade chain reactions that unfold over days or weeks.