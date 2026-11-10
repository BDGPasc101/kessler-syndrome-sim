# Kessler Syndrome Simulator
Stochastic modeling of orbital debris cascade and satellite lifespan prediction.

Catastrophic breakup ejecta is tracked as two debris populations: 1 cm shrapnel
and 10 kg structural chunks. Their sizes, masses, and cross-sectional areas
drive distinct atmospheric-drag rates; the dashboard reports heavy debris
separately in its inventory metrics and population area chart.

Active satellites use a configurable 100 kJ Whipple-shield energy cutoff and
30% shield coverage by default. Sub-threshold shrapnel hitting shielded areas
is absorbed into satellite mass and adds configurable lifespan wear (0.1 years
per hit by default); hits to exposed areas cause a mission kill and convert
the satellite to a derelict without a catastrophic breakup. Lifespan wear is
represented as fractional cohort retirements because individual satellite
ages are not tracked. Heavy debris remains catastrophic. Dashboard breach risk
uses heavy debris and derelicts only; shrapnel is excluded from that estimate.

Annual launches are distributed across 50 km shells in four LEO regimes:
5% to 300-400 km, 60% to the 450 km Starlink Highway, 25% to 500-600 km,
and 10% to the 650-850 km Legacy Grid. The Legacy Grid starts with 200
historical derelicts per shell at 700, 750, 800, and 850 km. Dynamic Highway
routing is enabled by default: if lethal-impact survival at 450 km over the
configured constellation lifespan plus five disposal years falls below 95%,
the excess risk above the 5% threshold is deducted from the Highway share and
reallocated evenly across the Expansion Shelf. Risk uses heavy-debris and
derelict spatial density, a 10 m^2 reference cross-section, and the altitude-
dependent relative-velocity model. The dashboard reports launches diverted in
the final simulated year; its viability-chart altitude selector is independent
of launch routing.

Atmospheric density and orbital decay vary with an illustrative 11-year solar
cycle. The mean density profile is scaled by a sinusoidal multiplier ranging
from 0.5 at solar minimum to 1.5 at solar maximum, with solar maximum at
simulation year 0 and every 11 years afterward. Drag is evaluated at each
adaptive substep's midpoint; the step-size bound accounts for the maximum
solar-cycle density over the candidate interval. Solar maxima are highlighted
on the debris inventory and mission-viability charts. Debris mass transported
downward remains accounted for, with material leaving the lowest modeled shell
recorded as reentered mass.

Collision ejecta is spread across the impact shell and its two neighbors in
each direction using a discrete Gaussian-inspired kernel: 60% remains in the
impact shell, 15% goes to each adjacent shell, and 5% goes to each shell two
steps away. Ejecta below the lowest modeled shell is recorded as reentered;
ejecta above the highest shell is retained in the top shell. Shrapnel and heavy
chunks use the same normalized weights so the spread preserves mass.

Active satellites use collision-avoidance maneuvers (CAMs) against tracked
heavy debris and derelicts, but not against 1 cm shrapnel. CAM warnings are
screened within a configurable radius (0.75 km by default), distinct from the
10 m^2 physical collision area. The default 90% avoidance efficiency uses a
finite reserve of 75 burns per satellite. Exhausted satellites remain active
but cannot maneuver; at retirement they fail PMD and become derelicts. Reserve
and screening are represented at aggregate shell-cohort level, and the burn
counter does not change modeled satellite mass. Derelict-derelict collisions
use a one-half symmetry factor to avoid counting identical-species pairs twice.
