# CSG FW-CADIS smoke test

This is the Path-A isolation test from the 2026-08-21 FW-CADIS handoff. It
uses OpenMC's canonical three-region CSG random-ray model with two detector
tallies and local FW-CADIS targets. It deliberately excludes DAGMC, ParaStell,
continuous-energy MGXS generation, and the stellarator source.

A pass establishes that the selected Ginsburg build can run the native
multigroup random-ray forward/adjoint solve and write a valid
`weight_windows.h5`. A failure therefore belongs to the core build/API; a pass
moves diagnosis to the DAGMC/CE-to-MG integration in Path B.

The batch job is staged and submitted through the AdvancedAgents Ginsburg
protocol. Outputs remain in the configured scratch root.
