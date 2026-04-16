# xtal Module — Multi-Agent Workflow

## Dependency Graph

```
                    ┌─────────────┐
                    │  element.rs  │ (already exists)
                    └──────┬──────┘
                           │
    ┌──────────────────────┼──────────────────────┐
    │                      │                      │
    ▼                      ▼                      ▼
┌──────────┐      ┌──────────────┐      ┌───────────┐      ┌───────────┐
│ types.rs │      │form_factors.rs│     │ bessel.rs │      │ fft_cpu.rs│
│          │      │              │      │           │      │           │
│ UnitCell │      │ IT92 table   │      │ I0, I1/I0 │      │ 3D FFT   │
│ Symop    │      │ vdW radii    │      │ log(I0)   │      │ grid dims│
│ SpaceGrp │      │              │      │           │      │           │
│ Reflect  │      └──────┬───────┘      └─────┬─────┘      └─────┬─────┘
│ DensGrid │             │                    │                  │
└────┬─────┘      ┌──────┴────────────────────┴──────────────────┘
     │            │
     ├────────────┼──────────────────────────────┐
     │            │                              │
     ▼            ▼                              ▼
┌──────────┐  ┌──────────────┐  ┌──────────┐  ┌───────────┐
│density.rs│  │solvent_mask.rs│  │scaling.rs│  │sigma_a.rs │
│          │  │              │  │          │  │           │
│ splat    │  │ F-H EDT      │  │ Lev-Marq │  │ binned D  │
│ symm_sum │  │ mask+shrink  │  │ Wilson   │  │ sigma²    │
│ deblur   │  │              │  │          │  │           │
└────┬─────┘  └──────┬───────┘  └────┬─────┘  └─────┬─────┘
     │               │               │              │
     └───────────────┴───────┬───────┴──────────────┘
                             │
                    ┌────────┴────────┐
                    │                 │
                    ▼                 ▼
            ┌──────────────┐  ┌───────────┐
            │map_coeffs.rs │  │ targets.rs│
            │              │  │           │
            │ 2mFo-DFc     │  │ ML target │
            │ mFo-DFc      │  │ B gradient│
            │ FOM          │  │ R-factors │
            └──────┬───────┘  └─────┬─────┘
                   │                │
                   └────────┬───────┘
                            │
                    ┌───────▼───────┐
                    │    mod.rs     │
                    │               │
                    │ XtalRefine    │
                    │ compute_map   │
                    │ refine_bfac   │
                    │ macro-cycles  │
                    └───────────────┘
```

---

## Stage 0: Setup (sequential, 1 agent)

**Agent: setup**

Update Cargo.toml, create `src/xtal/mod.rs` skeleton with feature gate, add module
declaration to `lib.rs`.

```toml
# Add to Cargo.toml [features]
xtal = ["dep:rustfft", "dep:argmin", "dep:argmin-math"]

# Add to [dependencies]
rustfft = { version = "6", optional = true }
argmin = { version = "0.10", optional = true }
argmin-math = { version = "0.4", features = ["vec"], optional = true }
```

```rust
// Add to lib.rs
#[cfg(feature = "xtal")]
pub mod xtal;
```

Create empty files for all submodules:
- `src/xtal/mod.rs`
- `src/xtal/types.rs`
- `src/xtal/form_factors.rs`
- `src/xtal/bessel.rs`
- `src/xtal/fft_cpu.rs`
- `src/xtal/density.rs`
- `src/xtal/solvent_mask.rs`
- `src/xtal/scaling.rs`
- `src/xtal/sigma_a.rs`
- `src/xtal/map_coefficients.rs`
- `src/xtal/targets.rs`

**Input:** Cargo.toml, lib.rs
**Output:** Feature-gated xtal module compiles (empty)
**Validation:** `cargo check --features xtal`

---

## Stage 1: Foundation (4 agents in parallel)

These four files have ZERO interdependencies. Each agent reads only `xtal-spec.md`
and the relevant research output.

### Agent 1A: `types.rs`

**Implement:** UnitCell, Symop, GroupOps, SpaceGroup, CrystalSystem, Reflection,
DensityGrid, epsilon_factor, is_centric, is_systematically_absent, grid_factors,
requires_equal_uv, round_up_to_smooth, has_small_factorization.

**Data to include:**
- All 10 space group symop arrays (from symops research agent output)
- Centering vectors for P, C, I lattices
- Grid factor match statement
- UnitCell::calculate_properties() with all matrix formulas
- UnitCell::d_star_sq(h,k,l)
- UnitCell::orthogonalize() / fractionalize()

**Spec section:** xtal-spec.md §1, §4 (grid dims)

**Lines estimate:** ~800-1000

**Tests to include:**
- UnitCell for orthorhombic cell (alpha=beta=gamma=90): verify orth matrix is diagonal
- UnitCell for monoclinic (beta≠90): verify orth[0][2] ≠ 0
- UnitCell volume for known cell: e.g., lysozyme P43212 a=b=79.1, c=37.9 → V=237,000 ų
- d_star_sq for (1,0,0) = (1/a)²
- epsilon_factor for general reflection = num_cen_ops
- epsilon_factor for special reflection (0,0,l) in P43212 > num_cen_ops
- is_systematically_absent for centering violations (h+k odd in C lattice)
- Symop::apply_to_frac roundtrip: apply identity → same coords
- Symop::apply_to_hkl: verify P212121 maps (1,0,0)→(-1,0,0) via second op
- round_up_to_smooth: 7→8, 11→12, 13→15, 31→32

### Agent 1B: `form_factors.rs`

**Implement:** FormFactor struct, FORM_FACTORS table (18 elements), VDW_RADII table,
lookup functions `form_factor(Element) -> &FormFactor`, `vdw_radius(Element) -> f64`.

**Data to include:**
- IT92 4-Gaussian coefficients for all 18 elements (from form factors research output)
- Bondi/Mantina vdW radii for all 18 elements

**Spec section:** xtal-spec.md §2

**Lines estimate:** ~150

**Tests to include:**
- sum(a) + c ≈ Z for each element (f(0) = atomic number)
- f(s=0.5) < f(s=0) for all elements (form factor decreases with s)
- Specific known values: C at s=0 should be ~6.0
- vdw_radius(Element::C) == 1.70

### Agent 1C: `bessel.rs`

**Implement:** bessel_i0(x), bessel_i1_over_i0(x) (the "sim" function),
log_bessel_i0(x) (the "sim_integ" function).

**Spec section:** xtal-spec.md §8

**Lines estimate:** ~120

**Tests to include:**
- I0(0) = 1.0
- I0(1.0) ≈ 1.2660658... (known value)
- I0(3.75) tests both branches of piecewise approximation
- I1(x)/I0(x) at x=0 → 0.0
- I1(x)/I0(x) at large x → approaches 1.0
- I1(x)/I0(x) at x=1.0 ≈ 0.4400506... (known)
- log(I0(0)) = 0.0
- log(I0(x)) ≈ x for large x (asymptotic)
- Relative error < 1e-6 across full range

### Agent 1D: `fft_cpu.rs`

**Implement:** fft_3d_forward (real grid → complex), fft_3d_inverse (complex → real grid),
using rustfft. Handle stride-based 1D passes along each axis.

**Spec section:** xtal-spec.md §4

**Lines estimate:** ~200

**Tests to include:**
- Round-trip: forward then inverse recovers original (within tolerance)
- Known 1D DFT: delta function → constant, constant → delta
- 3D: single point at origin → uniform in reciprocal space
- Grid dimension handling: works for smooth numbers (96, 108, 120)
- Normalization convention: document whether IFFT divides by N or not (rustfft does NOT)

---

## Stage 2: Core Algorithms (4 agents in parallel)

Each depends on Stage 1 types but NOT on each other.

### Agent 2A: `density.rs`

**Implement:** Precalculate per-atom Gaussian terms, cutoff_radius(), splat_density()
(atom loop + grid accumulation), symmetrize_sum(), compute_blur(), deblur_fc().

**Dependencies:** types.rs (UnitCell, SpaceGroup, DensityGrid), form_factors.rs (FormFactor)

**Spec section:** xtal-spec.md §3

**Lines estimate:** ~350

**Tests to include:**
- Single carbon atom at origin in P1, large B: density at origin > 0, decays with distance
- Symmetrize_sum in P212121: 4 equivalent positions have equal density
- Deblur: apply blur then deblur = identity (within tolerance)
- Total integrated density ≈ sum of atomic numbers (conservation of electrons)
- Cutoff radius increases with B-factor

### Agent 2B: `solvent_mask.rs`

**Implement:** edt_1d(), edt_1d_periodic(), edt_3d(), solvent_mask() pipeline
(mark protein → EDT → shrink threshold).

**Dependencies:** types.rs (UnitCell, DensityGrid), form_factors.rs (vdw_radius)

**Spec section:** xtal-spec.md §5, grid factors + solvent mask research output

**Lines estimate:** ~300

**Tests to include:**
- EDT of [0, INF, INF, INF, 0] → [0, 1, 4, 1, 0] (squared distances)
- EDT periodic: [0, INF, INF] → [0, 1, 1] (wraps around)
- 3D EDT: single point source, verify distance at known offsets
- Solvent mask: single atom in large box → spherical protein region
- Solvent fraction for typical protein ≈ 40-55%
- Shrink step reduces protein volume

### Agent 2C: `scaling.rs`

**Implement:** Levenberg-Marquardt optimizer (small, self-contained), Wilson plot initial
estimates, analytical derivatives (dy/dk_overall, dy/dk_sol, dy/dB_sol, dy/dB*),
B* symmetry constraints per crystal system, fit_scaling() pipeline.

**Dependencies:** types.rs (UnitCell, Reflection, CrystalSystem)

**Spec section:** xtal-spec.md §6

**Lines estimate:** ~500

**Tests to include:**
- LM solver on simple quadratic: converges to known minimum
- Wilson plot: synthetic data with known k and B → recovers them
- Constraint basis: orthorhombic has 3 free params, cubic has 1
- B* roundtrip: B_cart → B* → B_cart via orth/frac matrices
- Scaling synthetic data: known k_sol=0.35, B_sol=50 → fit recovers them within 10%
- Derivatives: numerical gradient check (finite difference vs analytical)

### Agent 2D: `sigma_a.rs`

**Implement:** SigmaAResult struct, estimate_sigma_a() (binning, D estimation,
sigma² computation, merging, smoothing, monotonicity, interpolation), r_free().

**Dependencies:** types.rs (Reflection), bessel.rs (for figure_of_merit used downstream)

**Spec section:** xtal-spec.md §7, sigmaA research output

**Lines estimate:** ~250

**Tests to include:**
- Perfect model (Fc = Fo): D ≈ 1.0 in all bins
- Random Fc: D ≈ 0.0
- Monotonicity: D[i] <= D[i-1] after enforcement
- Empty bins get interpolated, not NaN
- sigma² is always positive
- R-free computation: known synthetic data → expected R value

---

## Stage 3: Integration (2 agents in parallel)

### Agent 3A: `map_coefficients.rs`

**Implement:** figure_of_merit() (acentric/centric), compute_map_coefficients()
(2mFo-DFc and mFo-DFc), apply inverse FFT to get density grids.

**Dependencies:** types.rs, bessel.rs, sigma_a.rs (SigmaAResult), fft_cpu.rs

**Spec section:** xtal-spec.md §9

**Lines estimate:** ~200

**Tests to include:**
- FOM at x=0 → 0, FOM at large x → 1
- With D=1, sigma²→0: 2mFo-DFc ≈ 2Fo-Fc (classical coefficients)
- Missing Fobs: F_best = D*Fc, F_diff = 0
- Map coefficient computation is per-reflection (no inter-reflection dependencies)
- Inverse FFT of map coefficients produces real-valued grid

### Agent 3B: `targets.rs`

**Implement:** ml_target_value() (negative log-likelihood, acentric + centric),
r_work(), r_free(), BFactorGradient struct (compute B-gradient map coefficients,
IFFT, trilinear interpolation at atom positions), BFactorProblem (argmin CostFunction
+ Gradient traits with RefCell caching, sigmoid reparameterization for bounds).

**Dependencies:** types.rs, bessel.rs, sigma_a.rs, fft_cpu.rs, density.rs,
solvent_mask.rs, scaling.rs, map_coefficients.rs

**Spec section:** xtal-spec.md §10-11, B-factor gradient research output

**Lines estimate:** ~450

**Tests to include:**
- ML target decreases when Fc is moved toward Fobs
- R-work = 0 when Fc == Fo exactly
- B-gradient: finite difference check (perturb B by epsilon, recompute target)
- Sigmoid reparameterization: B stays in [2, 300]
- argmin integration: SteepestDescent on simple quadratic converges
- Trilinear interpolation: known linear field → exact interpolation

---

## Stage 4: Module Wiring (1 agent, sequential)

### Agent 4: `mod.rs`

**Implement:** XtalRefinement struct (holds all state), public API methods:
- `XtalRefinement::new(atoms, reflections, unit_cell, spacegroup)` — initialize
- `compute_map(atoms) -> DensityGrid` — full pipeline: splat → FFT → mask → scale → sigma_a → map_coeffs → IFFT
- `refine_b_factors(atoms) -> Vec<f32>` — one macro-cycle of B-factor refinement
- `fit_scaling(atoms)` — fit k_overall, B_aniso, k_sol, B_sol
- `update_sigma_a(atoms)` — re-estimate D and sigma²
- `r_factors(atoms) -> (f32, f32)` — R-work and R-free
- `refine(atoms, n_macro_cycles)` — full refinement loop

Also: module-level documentation, re-exports of key types.

**Dependencies:** ALL previous files

**Spec section:** xtal-gpu-plan.md §8 (API), §9 (macro-cycle structure)

**Lines estimate:** ~300

**Tests to include:**
- End-to-end: load a small test structure (manually constructed or from test data),
  compute map, verify R-free < 0.5 (better than random)
- Macro-cycle: R-free decreases over 3 cycles
- compute_map produces a non-zero grid

---

## Stage 5: Parallel Linting + Compliance (4 agents in parallel)

Each agent runs `cargo check --features xtal` and fixes issues iteratively.

### Agent 5A: clippy compliance
```
cargo clippy --features xtal -- -D warnings
```
Fix all clippy lints. The crate has strict clippy config (deny all, pedantic, nursery).

### Agent 5B: missing docs
```
cargo doc --features xtal --no-deps
```
Fix all `missing_docs` lint errors. Every public item needs a doc comment.

### Agent 5C: unused results + dead code
```
cargo check --features xtal 2>&1 | grep -E "warning|error"
```
Fix `unused_results`, `unused_qualifications`, etc.

### Agent 5D: test runner
```
cargo test --features xtal
```
Run all tests, fix any failures.

---

## Stage 6: Documentation (2 agents in parallel)

### Agent 6A: Module-level docs
Add `//!` module documentation to each file explaining what it does, the algorithm,
and key references.

### Agent 6B: Integration test
Create `tests/xtal_integration.rs` with a synthetic test case:
- Construct a small P1 structure (5-10 atoms) with known positions and B-factors
- Generate synthetic Fobs from a slightly perturbed model
- Run the full refinement pipeline
- Assert R-free improves

---

## Execution Summary

| Stage | Agents | Parallel? | Depends on | Est. total lines |
|-------|--------|-----------|------------|------------------|
| 0: Setup | 1 | no | — | ~30 |
| 1: Foundation | 4 | **yes** | Stage 0 | ~1300 |
| 2: Core | 4 | **yes** | Stage 1 | ~1400 |
| 3: Integration | 2 | **yes** | Stage 2 | ~650 |
| 4: Wiring | 1 | no | Stage 3 | ~300 |
| 5: Lint/Test | 4 | **yes** | Stage 4 | fixes only |
| 6: Docs | 2 | **yes** | Stage 5 | ~200 |
| **Total** | **18 agents** | | | **~3900 lines** |

## Agent Prompt Template

Each code-generation agent receives:
1. This workflow doc (their stage + dependencies)
2. `xtal-spec.md` (the relevant section)
3. The research output for their specific topic (if applicable)
4. The existing molex crate conventions (strict clippy, `missing_docs = "deny"`, etc.)
5. Instruction: "Write the file. Include tests in a `#[cfg(test)] mod tests` block.
   Follow molex's lint config. Do NOT add dependencies not listed in Stage 0."

## Files Produced

```
src/xtal/
  mod.rs              (Stage 4)
  types.rs            (Stage 1A)
  form_factors.rs     (Stage 1B)
  bessel.rs           (Stage 1C)
  fft_cpu.rs          (Stage 1D)
  density.rs          (Stage 2A)
  solvent_mask.rs     (Stage 2B)
  scaling.rs          (Stage 2C)
  sigma_a.rs          (Stage 2D)
  map_coefficients.rs (Stage 3A)
  targets.rs          (Stage 3B)
tests/
  xtal_integration.rs (Stage 6B)
```
