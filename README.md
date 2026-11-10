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
