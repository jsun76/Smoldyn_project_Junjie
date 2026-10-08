# Filament excluded volume implementation plan

Date: 2026-10-02. Status: design for review; collision code has not been implemented or benchmarked.

## Objective

Give filament segments finite physical volume, transmit contact loads through branched networks, and retain useful simulation speed. Establish correctness and convergence before interpreting network height or force as biological predictions. Excluded volume alone is not expected to correct the existing branching kinetics, load-independent elongation, or boundary control.

The recommended first implementation is capsule geometry, spatial bins, a cached neighbor list, and a finite harmonic repulsion integrated with the existing overdamped dynamics. Decide whether to replace the mechanical integrator from measured overlap, convergence, and runtime results. Do not assume that a generic rigid-body engine or a larger spring constant solves this problem.

## What the current checkout does

- `source/Smoldyn/smoldyn.h`: boxes store molecules and surface panels, but no filament segments. There is an existing `filradius` field, which must be audited for units, parsing, and compatibility before reuse. Drawing thickness is a separate parameter.
- `source/Smoldyn/smolboxes.c`: `boxesupdatelists` chooses a dense Cartesian grid from molecule count or requested box size. With no molecules and automatic sizing, each dimension has one box. `line2nextbox` offers segment traversal machinery; its face/edge/corner and boundary behavior needs targeted validation.
- `source/Smoldyn/smolfilament.c`: `filSegmentXFilament` already computes segment distances in an exhaustive search, but it is not part of the mechanical force path. `filComputeForces` adds stretching, bending, thermal, junction, and surface forces without steric forces.
- `filEulerDynamics` computes and moves one filament at a time. Pair forces must instead be evaluated at one common configuration, accumulated, and followed by movement of all participating filaments.
- `filPinBranches` translates entire daughters onto their branch points after motion. This can introduce overlaps and does not provide a coupled positional constraint that transfers arbitrary contact loads to the mother.
- `filAddJunctionForces` includes reactions for the polar angle spring, evaluated sequentially, but its azimuthal spring is explicitly one-way. Audit this when validating network force and torque balance.
- `filImplicitDynamics` solves a system for one filament and uses a local numerical Jacobian; it does not already provide the cross-filament blocks needed for implicit contacts. The older matrix solvers omit junction and confinement forces and ignore node mobility.
- Thermal forces are cached per simulation time and use the outer `sim->dt`. Mechanical substepping needs an explicit substep/noise interface; repeatedly calling the existing routine is insufficient.

## Physical and geometric model

Represent each straight segment as a capsule: its centerline plus a radius. Start with an explicitly documented actin diameter of 7 nm (`radius = 0.0035` in micrometre units), and sweep a plausible effective diameter rather than treating one value as exact. Associated proteins and electrostatics may justify an effective diameter later; do not inflate it just to obtain the desired network height.

For segments A and B, compute closest centerline points p and q, their interpolation fractions s and t, and distance d. Define penetration

    delta = max(0, rA + rB - d)
    U = 0.5 * k_contact * delta^2
    F = k_contact * delta * normal

Apply +F to A and -F to B, distributed as `(1-s, s)` and `(1-t, t)` to their endpoints. Verify force and torque balance and agreement with a numerical energy gradient. The existing `Geo_NearestSeg2SegDist` returns only distance; add a robust closest-point variant with fractions, scale-aware tolerances, and explicit treatment of degenerate and nearly parallel segments.

Exact coincident centerlines require a deterministic, documented normal or an invalid-initial-configuration diagnostic. Silently dividing by d or choosing a fresh random direction is unacceptable. Nearly parallel overlapping segments require testing for sufficient contact coverage; a single closest point can be inadequate for long side-by-side contacts. Avoid counting shared-endpoint contacts repeatedly and ensure results converge under segment refinement.

The first physical model is repulsive and frictionless. Friction, adhesion, bundling, crosslinks, and hydrodynamic interactions are additional assumptions requiring separate support.

### Topology and boundaries

Exclude a segment against itself and the intended overlap between consecutive bonded segments. Check nonlocal self-contact. Limit mother-daughter exemptions to the geometry near their actual shared junction; their distant portions must still collide. Test nearby distinct branches rather than exempting a whole branch family. Validate the junction exemption under changes in segment spacing and branch angle.

Make confinement radius-aware. For a planar floor or ceiling, centerlines must remain at least one radius inside the allowed region, up to the chosen penalty tolerance. Extend this to finite panel edges with capsule-panel geometry; endpoints alone are insufficient for arbitrary panels. The existing reflective domain walls for molecules are not a substitute for filament boundary handling. Initially support the project's nonperiodic planar chamber explicitly; unsupported steric/periodic combinations should produce a clear diagnostic until periodic images are implemented and tested.

Growth and branch creation must query occupied space before introducing a full new segment. The current 10 nm growth increment can jump across a 7 nm obstacle. Use a swept geometry check or smaller physical growth increments. Define blocked growth semantics explicitly, including the treatment of accumulated elongation length. Rejecting an obstructed branch changes the realized branching rate: record attempted, accepted, and blocked events. Do not silently choose repeated random directions until every event succeeds.

## Spatial indexing and neighbor lists

Follow Steven Andrews's ownership proposal: spatial cells hold filament segment references; segments need not hold box memberships. Use stable filament identity plus segment identity/index, with explicit invalidation after allocation, deletion, copying, branching, elongation, treadmilling, or any other topology change. Avoid stale pointers when arrays expand or segment indices move.

Choose and document one complete broad phase:

1. Store centerlines in all crossed cells and search a cell stencil expanded by both radii plus skin; or
2. Insert radius/skin-expanded segment bounds in all intersected cells and consider shared-cell pairs.

Both can work. Centerline-only insertion followed by same-cell comparisons is incomplete. The search extent must account for actual cell widths and mixed radii, rather than assuming 26 neighbors always suffice. Validate line traversal when a segment lies on a face or crosses an edge/corner.

Assess existing box occupancy first. Integrate segment lists into existing boxes when their resolution is suitable. If they are too coarse, use a sparse filament subdivision associated with occupied boxes, or a separate sparse bin index sharing the same geometry conventions. Keep molecule grid tuning independent from the filament contact cutoff. The current 4 x 4 x 10 micrometre domain would need about 20 million cells at 20 nm spacing if made dense; allocating full box structures at that resolution is unsuitable.

Start a cell-width sweep around 10, 20, and 40 nm for the 10 nm discretization. These are benchmark choices, not correctness requirements or fixed defaults. Keep storage contiguous where practical, reuse capacity, clear only occupied cells, and avoid allocation inside distance/force loops.

Generate a deterministic unique pair list within `rA + rB + skin`. Cells may discover a pair multiple times; sort/unique stable pair identifiers or use an equivalent proven ownership rule. Perform inexpensive bounding-box rejection before exact segment distance. Recompute forces every mechanical step even when the neighbor list is reused.

Rebuild when any endpoint displacement exceeds half the skin, or immediately after topology/box/panel/constraint changes that can invalidate coverage. Monitor every endpoint, including motion caused by junction corrections. Endpoint displacement bounds the displacement of all points on a straight segment while topology is fixed. Start with a 2-5 nm skin sweep and measure the rebuild/candidate tradeoff. Do not use a fixed rebuild interval without a displacement guarantee.

Expected cost is O(number of segments + number of candidates), at bounded local density. Dense packing can increase candidate count substantially; this is not an unconditional O(N) guarantee. At 14,766 segments an exhaustive comparison would require 109,009,995 unordered pairs each evaluation.

## Mechanical solver and stability

The current mechanics is overdamped, approximately

    dx = M F(x) dt + sqrt(2 kT M) dW.

There are no inertial collision impulses. Velocity Verlet, rigid-body bounce, or a conventional semi-implicit velocity update is not a direct replacement.

First make the enabled steric path global: topology update, common geometry snapshot, force initialization for all filaments, internal/junction/wall forces, pair forces once, consistent attachment constraints, all-node update, geometry refresh, and contact diagnostics. Respect immobile nodes and mixed node/type mobilities. Avoid force clearing after pair contributions have been accumulated. Force stages for RK and implicit integration also require a common network state; enable only validated integrators and reject unsupported combinations.

Branch attachment is a prerequisite for credible load transmission. Evaluate either shared physical junction coordinates or a mobility-weighted coupled constraint with reaction forces. Account consistently for drag and thermal noise at constrained/merged coordinates; assigning multiple independent Brownian kicks to one physical shared coordinate is incorrect. Validate the positional attachment separately from the existing angular springs. Do not preserve whole-daughter translation as the final steric correction mechanism.

For an isolated relaxing scalar spring, forward Euler requires `mu*k*dt < 2`. For two movable point contacts, relative relaxation involves `(muA+muB)*k`, not just one mobility. For a network, the relevant bound involves the largest eigenvalue of the mobility-weighted stiffness matrix, including stretch, bend, contacts, junctions, and walls. Contact count and node interpolation weights matter. The file's scalar `mobility*force_length*dt` comment is not a complete network guarantee.

At the present `mu=1`, `dt=1e-5`, a trial `k_contact=4000` gives the isolated equal-mobility contact factor 0.08. This is a feasibility estimate, not a stability proof. The thermal penetration scale is approximately `sqrt(kT/k_contact)`, or 0.5 nm for `kT=0.001`. Load-induced penetration is approximately `F/k_contact` and must also meet the tolerance. Sweep 1000, 4000, and 10000 in model stiffness units and compare measured penetration distributions.

Reducing mobility reduces deterministic relaxation speed and diffusion (`D=kT*mu`). It can preserve equilibrium distributions in an appropriately converged thermal model but changes mechanical kinetics relative to growth and capping. The claim that mechanical relaxation remains faster than chemistry must be measured for relevant network modes, especially long collective modes. The comment that the physical mobility is about 10000 is not itself a calibration: derive node drag from viscosity, filament dimensions, discretization, and boundary assumptions.

At the current settings the unconstrained Brownian displacement scale per coordinate is approximately 0.14 nm per step. Increasing mobility by 10000 at unchanged dt increases this to about 14 nm, while multiplying deterministic stiffness factors by 10000. This is strong reason to redesign/tune integration before increasing mobility; it does not establish that mobility 1 is adequate biologically.

### Escalation path

If a finite penalty meets penetration and convergence targets without excessive runtime, retain it. If not, prototype a network-level linearly implicit or constrained overdamped solver, following the Cytosim approach. Include cross-filament contact/junction blocks and use sparse or matrix-free iterative solves with residual diagnostics. Reuse structure/preconditioners when appropriate. Compare simulated physical time per wall-clock second, not cost per step alone.

Mechanical substeps can separate fast mechanics from slower chemistry, but must evaluate every stiff mechanical force at the substep. Specify substep Brownian increments consistently, and use Brownian bridges if stochastic steps are rejected/subdivided. A hard distance projection/XPBD approach is a comparison candidate; its thermal statistics and physical contact force estimates require validation. Ordinary positional projection or a converged implicit solve does not automatically produce accurate Brownian sampling or prevent geometric tunneling.

## Staged work and acceptance gates

1. **Baseline and diagnostics.** Record revision/build/configuration hashes. Measure fixed-topology sparse and dense snapshots separately from growing runs. Add timing for index, pair construction, exact distances, force evaluation, constraints, and integration; record candidates, contacts, rebuilds, memory, penetration, wall violations, junction error, segment strain, and energies. Quantify relaxation relative to chemistry. Use multiple random seeds for statistical results.
2. **Geometry and indexing.** Implement capsule queries and box segment indexing in diagnostic mode. Compare every candidate/contact with a brute-force oracle on small randomized systems and adversarial cell boundaries. Demand no missed contacts and no duplicated forces. Verify allocation/free and topology invalidation. Diagnose pre-existing overlaps without silently changing initial state.
3. **Coupled mechanics.** Implement the common-state force/update path and consistent branch position constraints. Test a force applied to a daughter, its transmission to the mother, force/torque balance, attachment error, free diffusion, equilibrium fluctuations, and immobile nodes. Audit azimuthal reaction behavior.
4. **Steric forces, boundaries, and creation.** Enable finite penalties with explicit physical radii. Test skew, parallel, endpoint, coincident, and nonlocal self-contact; local branch exemptions; radius-aware floor/ceiling; growth into an obstacle; and creation in crowded regions. Test contacts after attachment correction, not only before it.
5. **Convergence and optimization.** Sweep dt, stiffness, skin, cell width, radius, and spacing. Recompute bend stiffness and mobility consistently when spacing changes. Compare fixed-topology timings for 1k, 5k, 15k, and 50k segments at controlled density, then representative growing networks. Optimize the measured bottleneck and decide whether a different solver is necessary.
6. **Biological evaluation.** Compare max and percentile height, z profiles, footprint, contact/packing distributions, force during relaxation/hold, realized growth/branch/cap rates, and network relaxation across seeds. Sterics is not validated solely because height increases. Address NPF-limited branching, load-dependent chemistry, and force-controlled boundaries as separate work.

Provisional targets to discuss with the lab:

- No broad-phase omissions relative to the small-system oracle.
- Deterministic force/energy checks within numerical tolerance; no systematic net contact force or torque.
- For nonbonded contacts, penetration P99 below 1 nm in representative loaded runs, with excursions above 2 nm investigated. Report temperature, load, stiffness, and discretization alongside this target.
- Halving dt and tightening solver tolerance changes key ensemble observables by less than about 5%, or less than justified sampling uncertainty. Assess force and equilibrium fluctuation bias in addition to geometry.
- Aim for no more than 25% runtime overhead on the representative workload, measured at fixed topology with matched output; this is an optimization target, not a promise. Publish dense-case overhead and memory separately. If accurate sterics exceeds the budget, choose a validated solver/coarser discretization rather than suppressing contacts.
- Disabled collision mode avoids allocating its data and preserves the old path. Contact-enabled results are expected to differ, including realized growth and filament count.

## Proposed code organization

- `source/Smoldyn/smoldyn.h`: physical radius/contact configuration, box segment references, neighbor/cache and constraint work structures.
- `source/Smoldyn/smolboxes.c`: segment insertion, traversal, occupancy management, and query helpers.
- `source/libSteve/Geometry.c` and `.h`: closest segment points/fractions and capsule-panel support with independent geometry verification.
- `source/Smoldyn/smolfilament.c`: parameter parsing, lifecycle invalidation, global mechanics staging, growth/branch collision hooks, and diagnostics.
- Optional `source/Smoldyn/smolfilamentsteric.c`: contact pair building and force accumulation if this keeps `smolfilament.c` manageable; register it in the relevant CMake source lists.
- `examples/S13_filaments/`: small demonstrators and reproducible physical cases. Add targeted automated geometry/mechanics checks to an appropriate test target.

Proposed configuration names (not implemented): `steric_radius`, `steric_stiffness`, `steric_skin`, and a neighbor/index setting. Define pair mixing, units, solver compatibility, and output diagnostics explicitly before exposing the interface. Use the project's unit convention; physical radius must not inherit the default `filradius=1` inadvertently.

## Sources and how they inform the design

- [PhysX geometry](https://nvidia-omniverse.github.io/PhysX/physx/5.5.1/docs/Geometry.html): capsule primitives. Useful geometry reference; its inertial dynamics is not the target model.
- [LAMMPS neighbor lists](https://docs.lammps.org/Developer_par_neigh.html): spatial bins, unique pair lists, cutoff plus skin, displacement-triggered reuse. Adapt point indexing to extended moving segments.
- [Nedelec and Foethke, Collective Langevin Dynamics of Flexible Cytoskeletal Fibers](https://arxiv.org/abs/0903.5178), New Journal of Physics 9 (2007), 427: overdamped connected fibers, implicit integration, and length constraints. The closest solver reference for this application; convergence and noise checks remain necessary.
- [Macklin et al., XPBD](https://doi.org/10.1145/2994258.2994272): compliant constraint formulation that addresses conventional PBD timestep/iteration-dependent stiffness. A comparison option requiring adaptation to overdamped thermal mechanics.
- [Letort et al., 2015](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004245): distinguishes the approximately 7 nm actin diameter from enlarged effective interaction diameters used in some simulations. Do not copy such enlarged diameters without their physical assumptions.

No evidence in this checkout establishes whether prior mobility/solver compromises were comprehensively benchmarked. Existing comments mention particular measurements, and example files exercise filament features, but neither proves convergence of the present dense growing network. Judge these choices from reproducible sweeps, not the identity of the coding assistant.
