# Orbital Debris Cascade Simulator: Engineering & Physics Learnings

## 1. Physics Engine Architecture: How the Core Works

To build this model, we didn't track every single satellite's specific orbit—that would require a supercomputer. Instead, we used a "particle-in-a-box" approach, treating Earth's orbit like a series of traffic lanes.

### 1.1 Slicing the Sky (Shell Geometry)
We divided Low Earth Orbit (LEO) into 13 distinct "shells" or traffic lanes. Starting at 400 km and going up to 1,000 km, each shell is exactly 50 km thick. 

* **Spatial Density:** Every simulated year, the engine calculates the volume of each shell (in cubic kilometers) and counts the objects inside it. This gives us the "Spatial Density" (objects per cubic kilometer)—the single most important metric for determining collision risk.

### 1.2 Bumping in the Dark (Collision Mechanics)
Rather than simulating exactly when two objects cross paths, we use **Kinetic Gas Theory**—the same math used to predict how gas molecules bounce off each other in a balloon.

* If a shell has a high density of objects, and those objects have a large physical size (cross-section), and they are moving incredibly fast (relative velocity), they will hit each other more often.
* **Math Note:** We calculate average collision speeds based on an assumed 90° intersection angle, recognizing that objects rarely hit perfectly head-on or perfectly rear-end each other.

### 1.3 The Fireworks (NASA Breakup Model)
When things collide in space, they don't just crack in half—they vaporize and shatter. We use the official **NASA Standard Breakup Model** to determine what happens.

* **Heavy Hits (Catastrophic):** If the impact energy is massive, both objects completely shatter. A 250 kg satellite might instantly turn into 15,000 pieces of shrapnel.
* **Light Hits (Crater Impact):** If a tiny piece of shrapnel hits a big satellite, it just chips off a few pieces rather than destroying the whole thing.

### 1.4 The Vacuum Cleaner (Atmospheric Drag)
Satellites don't stay in orbit forever. The Earth's atmosphere extends out into space, creating friction (drag) that slowly pulls objects down.

* In our simulation, we model this as a "waterfall." Debris from the 800 km shell slowly leaks down into the 750 km shell, which leaks into the 700 km shell.
* Anything that falls below 400 km burns up completely in the dense lower atmosphere and is permanently deleted from the simulation.

---

## 2. Model Limitations (Where We Simplified)

Because this is a 1D thermodynamic model and not a full 3D orbital propagator, we had to make a few intentional compromises:

1. **Traffic Jams at the Poles:** We assume objects are spread perfectly evenly throughout a 50 km shell. In reality, orbits often cross over the North and South poles, creating high-density "pinch points" where collisions are much more likely.
2. **Instant De-Orbiting:** Real orbital mechanics dictate that if you slow down a piece of debris, it becomes an oval (elliptical) orbit that swoops down and back up. Early versions of our model just instantly deleted debris, ignoring the danger it poses while falling.
3. **No Swerving:** Our active satellites are "dumb." They do not use their thrusters to dodge incoming debris, even though real-world operators (like SpaceX) maneuver their satellites thousands of times a month to avoid collisions.

---

## 3. Phase 2 Upgrades: Fixing the Math and the Code

As the simulation became more complex, stress-testing it with extreme inputs revealed fascinating intersections between numerical mathematics and orbital physics.

### 3.1 Debugging Infinity: The ODE Singularity Bug
* **The Problem:** In a massive runaway Kessler scenario, the model crashed, reporting 1.49 × 10³² fragments—more than the number of atoms in Earth's atmosphere!
* **The Fix:** The math was evaluating total collision probability without a cap. If the risk was high enough, it was mathematically having a single satellite explode 10,000 times in one year. We clamped the math: a satellite can only be completely destroyed *once*.

### 3.2 Linear Lasers vs. Exponential Explosions
* **The Finding:** We maxed out our "Laser Broom" to fire 50,000 times a year to clean up the debris. Surprisingly, a century later, it had barely made a dent.
* **The Lesson:** A laser is a *linear* tool (you shoot one piece of debris at a time). The Kessler Syndrome is an *exponential* process (one collision makes 10,000 fragments, which cause 10,000 more collisions). You cannot "shoot your way out" of an exponential explosion. Mitigation must focus on preventing the cascade from starting.

---

## 4. Phase 3 Upgrades: Engineering Supreme Realism

To push the model beyond simple mathematics and turn it into a true "reactive" physical environment, we implemented four massive layer-by-layer physics upgrades.

### 4.1 The "Whipple Shield" Sponges
* **The Reality:** Satellites aren't totally defenseless. They are wrapped in thin, spaced-out aluminum and Kevlar layers called Whipple Shields that vaporize tiny shrapnel impacts.
* **The Code:** We separated our debris into "Heavy Chunks" and "Shrapnel". We introduced a Kinetic Energy Threshold. Now, if shrapnel hits an active satellite, it is safely absorbed. This turned active satellites into giant "sponges" that actually help sweep the orbit clean of tiny debris, vastly stabilizing the simulation!

### 4.2 Dynamic Intersection Speeds
* **The Reality:** Objects at 400 km travel faster than objects at 800 km.
* **The Code:** We stopped hardcoding the collision speed to 10 km/s. The engine now calculates the exact Keplerian orbital velocity for every specific shell, applies a mathematically rigorous 90° intersection geometry, and feeds that unique speed into the collision probability formulas.

### 4.3 Kepler's Second Law & Transiting Debris
* **The Reality:** When we knock debris out of orbit, it falls in an elliptical path. Objects in elliptical orbits spend far more time lingering at the top (apogee) than whipping around the bottom (perigee).
* **The Code:** Using the "eccentric anomaly" from orbital mechanics, the simulation now calculates exactly what fraction of an orbital period a falling piece of debris spends in each 50 km shell. It perfectly simulates debris dangerously lingering in the upper shells while quickly passing through the lower ones.

### 4.4 Dynamic Highway Routing (A Reactive Economy)
* **The Reality:** 60% of all real-world launches go to the 450 km "Starlink Highway." But if that orbit becomes a Kessler-choked death trap, aerospace companies won't keep launching there—they will move their satellites to a higher, safer orbit.
* **The Code:** We built a reactive traffic controller. The engine constantly calculates the lethal-impact risk at the 450 km shell. If survivability drops below 95%, it automatically reroutes the excess traffic up into the 500–600 km "Expansion Shelf." You can literally watch the space industry flee upward as the lower orbits become unusable, physically accelerating the destruction of the higher orbits!

### 4.5 Hypervelocity Ejecta Cross-Shell Bleed
* **The Reality:** Explosions in orbit do not stay confined to a single altitude. Hypervelocity impact kinetic energy imparts a delta-v dispersion that casts fragments into eccentric orbits spanning above and below the point of impact.
* **The Code:** Implemented a 5-shell discrete Gaussian kernel (60% parent shell, 15% at ±50 km, 5% at ±100 km). Mass ejected below the 300 km boundary is routed directly to the atmospheric reentry sink, while ejecta crossing the 1,000 km ceiling is clamped to the upper boundary, maintaining strict mass conservation across 14 altitude bands.

### 4.6 Conjunction Screening, CAMs, and Propellant Exhaustion
* **The Reality:** Satellites do not passively wait to be struck; operators perform Collision Avoidance Manoeuvres (CAMs) whenever an object enters an operational safety bubble (~0.75 km). However, onboard propellant is finite.
* **The Code:** Decoupled the conjunction screening cross-section ($\sigma_{\text{screen}} = \pi R_{\text{screen}}^2$) from physical satellite collision destruction ($\sigma_{\text{phys}} = 10\text{ m}^2$). Active cohorts reduce macroscopic lethal risks by 90% until their maneuver budget is exhausted. Propellant-depleted cohorts lose maneuver capability and suffer 100% Post-Mission Disposal (PMD) failure at retirement, creating a non-linear tipping point that naturally drives the runaway cascade.
