# GPU-Accelerated Crystallographic Refinement — Implementation Plan

## Purpose

This document captures every technical decision, architectural insight, and implementation detail from a deep research session on adding GPU-accelerated electron density map computation and ML crystallographic refinement to the molex/viso ecosystem. It is designed to be fed into a Claude Code session so no context is lost.

---

## 1. What We're Building

A feature in **molex** (the molecular IO/computation crate) that computes electron density maps from X-ray crystallographic data (Fobs + model coordinates) using the Maximum Likelihood (ML) refinement target function, with GPU acceleration via CubeCL's wgpu backend. The output is a real-space density grid that viso (the molecular renderer) consumes for isosurface extraction via marching cubes.

### Why This Matters

- Currently, electron density map computation is done by tools like Phenix/Refmac (CPU, Python/C++/Fortran, offline batch processing).
- Nobody has built an interactive-speed GPU-accelerated refinement pipeline.
- Amber 2023 proved the ML target function works on GPU (their `xray` module in `pmemd.cuda`), but it's designed for MD trajectories (hours-long runs), not interactive visualization. Their code is also restricted-license ($25k commercial).
- We're building the interactive version: student releases mouse → map updates with refined B-factors before they blink.

### What the Pipeline Does

Given: a set of observed structure factor amplitudes (|Fobs|) from an X-ray diffraction experiment, and a 3D atomic model (coordinates, B-factors, occupancies):

1. **Density splatting**: Place Gaussian atomic density at each atom position onto a 3D grid, using scattering form factors from the International Tables for Crystallography (ITC).
2. **FFT → Fc**: Transform real-space atomic density to reciprocal space to get calculated structure factors (Fc).
3. **Solvent mask**: Compute a binary mask of the bulk solvent region (Jiang & Brünger method — probe sphere + shrink), FFT it to get Fmask.
4. **Scaling**: Fit bulk solvent parameters (k_sol, B_sol) and overall scale (k_overall, B_aniso) by least-squares against Fobs. This gives Ftotal = k · exp(-B·s²) · (Fc + k_sol · exp(-B_sol·s²) · Fmask).
5. **σA estimation**: Estimate resolution-dependent D and σ parameters from Fobs/Fc agreement in ~15 resolution bins, using free reflections to avoid bias.
6. **Map coefficients**: Compute 2mFo-DFc and mFo-DFc weighted map coefficients using σA weights. These are the "best" map coefficients (Read 1986).
7. **Inverse FFT → density**: Transform map coefficients back to real space to get the electron density grid.
8. **Marching cubes**: Extract an isosurface mesh from the density grid (viso's job, not molex's).

For ML B-factor refinement, steps 1-7 are repeated iteratively (10-15 cycles), with B-factors updated each cycle using the gradient of the ML target function sampled at atomic positions.

---

## 2. The ML Refinement Target Function

The ML target minimizes:

```
-log L = Σ_h [ log(2|Fobs|/ε σ²) + (|Fobs|² + D²|Fc|²)/(2ε σ²) - log I₀(D|Fobs||Fc|/(ε σ²)) ]
```

Where:
- h indexes reflections (Miller indices)
- D and σ are resolution-dependent σA parameters (estimated in bins from free reflections)
- ε is the expected intensity factor for the reflection's symmetry
- I₀ is the modified Bessel function of the first kind, order zero

The gradient uses the Bessel function ratio I₁(x)/I₀(x), which is smooth and can be approximated by rational polynomials (~20 FLOPs per reflection in a GPU shader).

The gradient of the ML target IS the mFo-DFc difference map sampled at atomic positions. Per-atom B-factor gradients come from sampling this gradient map and multiplying by the atom's radial derivative ∂f/∂B.

### Key References for the Math

- Murshudov 1997 (ML refinement in Refmac): pubmed.ncbi.nlm.nih.gov/15299926/
- Pannu & Read 1996 (ML target function)
- Brünger 1997 (bulk solvent + scaling)
- Read 1986 (σA weighting for map coefficients)

---

## 3. Algorithmic Reference: Gemmi

**Gemmi** (C++14, MPL 2.0, by Marcin Wojdyr / Global Phasing + CCP4) has the full pipeline we need to reimplement:

- `DensityCalculatorX`: atomic density splatting with B-factors and scattering form factors
- `SolventMasker`: Jiang & Brünger method with rprobe/rshrink parameters
- `Scaling` class: fits k_overall, B_aniso, k_sol, B_sol via least-squares against Fobs
- FFT via PocketFFT
- Map coefficient computation

Gemmi compiles to WebAssembly (used by UglyMol, a web-based crystallography viewer). It's the best algorithmic reference because it's clean C++14 without the massive dependency graphs of CCTBX/Phenix or the CCP4 entanglement of Clipper.

**Important**: We are reimplementing the algorithms in Rust, using gemmi's C++ source as an algorithmic reference. We are NOT binding to gemmi. Attribute the algorithmic origin appropriately.

Gemmi docs for the relevant pipeline: gemmi.readthedocs.io/en/latest/scattering.html

### What Rosetta's Existing Code (ElectronDensity.cc) Does vs. What It's Missing

Rosetta's existing code computes ρC via Gaussian splatting with one-Gaussian form factors, FFTs to get Fc, extracts phases, combines Fobs amplitudes with model phases, and inverse FFTs. What it's MISSING:

- Bulk solvent correction (Fmask)
- Scale/B optimization (k_sol, B_sol fitting)
- σA weighting (2mFo-DFc map coefficients)

These are what separate a toy electron density viewer from a real crystallographic refinement tool.

---

## 4. Architecture

### Where Everything Lives

```
gpu-fft                    ← existing crate (by Eugene Hauptmann), dependency
                             MIT licensed, CubeCL-based, batched 1D FFT
                             github.com/eugenehp/gpu-fft

molex (feature = "xtal")   ← CPU crystallographic refinement (pure Rust)
molex (feature = "xtal-gpu") ← GPU-accelerated path
  ├── cubecl-wgpu          ← ALL GPU compute: FFT, splatting, mask, ML kernels
  └── xtal/
      ├── mod.rs
      ├── form_factors.rs    ← ITC scattering factor tables (public data, 5-Gaussian approximation)
      ├── density.rs         ← atomic density splatting onto grid (CubeCL kernel for GPU, loop for CPU)
      ├── solvent_mask.rs    ← bulk solvent mask, Jiang & Brünger method (CubeCL kernel / CPU)
      ├── scaling.rs         ← ML scaling: k_overall, B_aniso, k_sol, B_sol (least-squares fit)
      ├── sigma_a.rs         ← σA estimation in resolution bins (from free reflections)
      ├── map_coefficients.rs ← 2mFo-DFc, mFo-DFc computation (CubeCL kernel / CPU)
      ├── fft_cpu.rs         ← 3D FFT wrapper using rustfft (CPU fallback)
      ├── targets.rs         ← ML target function + gradients (∂L/∂B, ∂L/∂xyz)
      └── gpu/
          ├── fft.rs           ← CubeCL 3D FFT: butterfly kernel + batched 1D composition
          ├── density.rs       ← CubeCL density splatting kernel
          ├── solvent_mask.rs  ← CubeCL mask kernel
          ├── map_coeffs.rs    ← CubeCL ML weight + map coefficient kernel
          └── gradient.rs      ← CubeCL B-factor gradient sampling kernel

viso
  ├── molex                ← calls molex compute, receives buffer handles
  ├── wgpu                 ← render pipeline, marching cubes (raw WGSL stays here)
  └── compute/
      ├── prefix_sum.wgsl    ← reusable primitive (MC, SES, electrostatics)
      ├── stream_compact.wgsl
      ├── mc_classify.wgsl   ← marching cubes classify + count
      └── mc_vertices.wgsl   ← marching cubes vertex generation
```

### Feature Gating

```toml
# molex/Cargo.toml
[features]
default = []
xtal = ["rustfft", "argmin"]          # CPU path, pure Rust, no GPU deps
xtal-gpu = ["xtal", "cubecl-wgpu"]    # GPU-accelerated path
```

### Why GPU Compute Lives in molex, Not viso

The crystallographic computation is molex's domain. Viso is the renderer — it takes a density grid and renders it. Viso shouldn't know anything about structure factors, B-factors, or bulk solvent. The handoff is just a `wgpu::Buffer` containing a density grid.

### Device Sharing Between CubeCL and viso

CubeCL has a `WgpuDevice::Existing(u32)` variant, initialized via `init_device()`, that accepts an external wgpu device/queue setup. The docs explicitly state: "Useful when you want to share a device between CubeCL and other wgpu-dependent libraries." (Examples: egui, Bevy.)

This means:
- viso owns the `wgpu::Device` and `wgpu::Queue` (it creates them at startup)
- viso passes its device/queue to molex
- molex initializes CubeCL with `WgpuDevice::Existing`, sharing the same device
- CubeCL compute dispatches and viso render/compute passes share the same device
- Buffers produced by CubeCL kernels are accessible to viso's wgpu pipeline
- Zero CPU readback for the common case

**Critical verification needed**: confirm that CubeCL buffer handles can be unwrapped back to raw `wgpu::Buffer` for use in viso's bind groups. The `WgpuDevice::Existing` design implies this should work, but it needs to be tested.

### Why CubeCL for molex's Compute, but NOT for viso's Shaders

CubeCL is a compute-only abstraction. Viso's render pipeline (impostor shaders, deferred shading, SSAO, backbone geometry, pick buffers) uses vertex/fragment shaders — CubeCL can't express these. Viso's existing compute shaders (marching cubes, prefix sum) are tightly integrated with render passes and stay as raw WGSL.

CubeCL is appropriate for molex because molex's workload is pure compute with no rendering. The advantage: CubeCL kernels written in Rust via `#[cube]` macros compile to WGSL (for wgpu), CUDA (for NVIDIA), and HIP (for AMD) — molex gets native CUDA performance on NVIDIA hardware for free.

---

## 5. The FFT Strategy

### Custom CubeCL FFT Kernels (Buffer-Resident, Zero Readback)

**gpu-fft was evaluated and rejected.** The `gpu-fft` crate (github.com/eugenehp/gpu-fft) provides batched 1D FFT on CubeCL, but its API forces host round-trips (`Vec<f32>` in, `Vec<f32>` out). For a pipeline targeting <10ms per cycle, four FFTs × ~4ms memcpy overhead each = 16ms in transfers alone — this blows the entire budget before any computation. It also lacks complex-to-complex support and is power-of-two only.

**Instead, we write our own FFT butterfly kernel in CubeCL**, using gpu-fft's Cooley-Tukey radix-2/4 implementation as algorithmic reference. The key difference is the API surface: our kernel operates on CubeCL tensor handles (GPU-resident buffers) and returns CubeCL tensor handles. No host round-trips. The output buffer of one FFT feeds directly into the next compute kernel.

The actual kernel code is ~100-150 lines of CubeCL (the butterfly operation, twiddle factor computation, bit-reversal permutation). The 3D composition — three batched 1D passes with axis transposition — is ~100 lines of Rust orchestration.

A 3D FFT decomposes into three batched 1D FFT passes:

1. Batch of Ny×Nz 1D FFTs along axis 0
2. Transpose/restride (CubeCL kernel)
3. Batch of Nx×Nz 1D FFTs along axis 1
4. Transpose/restride (CubeCL kernel)
5. Batch of Nx×Ny 1D FFTs along axis 2

All passes stay on GPU. No data touches the CPU between passes.

### Power-of-Two Constraint

The initial CubeCL FFT implementation will be radix-2/4 (power-of-two only). For crystallographic grids with non-power-of-two dimensions (e.g., 135, 150), bias `findSampling` to prefer powers of two when the GPU path is active. This wastes some compute (padding 150 to 256) but keeps the kernel simple for v1. Mixed-radix (3, 5, 7) can be added later as a performance optimization.

### CPU Fallback

For web/WASM targets and systems without GPU compute support, use `rustfft` (pure Rust, well-optimized, handles mixed radix) for the FFT steps. Everything else in the pipeline has trivial CPU implementations (loops over atoms/grid points).

---

## 6. Performance Estimates

### Performance Target

The real comparison is **ISOLDE** (Coot-based interactive crystallographic refinement), which does live map recalculation at ~1.8 updates/second (~550ms each) on a 16-core CPU for ~3300 atoms. Our target: **<10ms per refinement cycle** for a 1000-residue protein, enabling 60fps-capable updates with breathing room. This is a ~55x speedup over ISOLDE.

We are NOT competing with Phenix/Refmac on final model quality after convergence — we converge to the same ML target, just faster. The goal is interactivity: student moves an atom → map updates before they lift their finger.

### Prior Art: SFCalculator (Hekstra Lab, 2025)

SFCalculator — a differentiable structure factor calculator in PyTorch/JAX — demonstrated 50-200x speedup over Phenix for structure factor calculation on an A100 GPU. This validates the core thesis that these computations are massively GPU-parallelizable. Different stack (PyTorch vs CubeCL) but identical math.

### Common Case: 1000-Residue Protein at 2.0 Å Resolution

Grid ~128×128×256 (power-of-two padded), ~80k reflections, ~8000 heavy atoms. All buffers GPU-resident, zero host round-trips.

| Step | CPU (Rust) | GPU (CubeCL, buffer-resident) |
|---|---|---|
| Density splatting | 15-30ms | ~0.3ms |
| FFT → Fc (3 batched 1D passes) | 10-20ms | ~0.5ms |
| Solvent mask | 10-20ms | ~0.2ms |
| FFT mask → Fmask | 10-20ms | ~0.5ms |
| Scaling fit (CPU, ~10 params) | 1-5ms | ~0.3ms (stays CPU) |
| σA update (CPU, binned stats) | <1ms | <0.1ms (stays CPU) |
| Map coefficients | <1ms | ~0.1ms |
| Inverse FFT → density | 10-20ms | ~0.5ms |
| Gradient sampling | 5-10ms | ~0.2ms |
| CubeCL dispatch overhead (~8 kernels × 65μs) | — | ~0.5ms |
| **Subtotal (raw compute)** | | **~3.2ms** |
| **With 2-4x wgpu efficiency pessimism** | | **~6-10ms** |
| **Total per cycle (realistic)** | **~60-120ms** | **~8ms** |

The wgpu efficiency gap (vs native CUDA) is real — WebGPU compute reaches ~17-25% of theoretical peak on Apple Silicon vs ~75% for CUDA. The 2-4x pessimism factor accounts for this. On NVIDIA hardware with CubeCL's CUDA backend, expect the raw ~3ms number.

Full ML B-factor refinement (5 macro-cycles × 3 inner iterations): ~120ms on GPU. Still interactive.

### Important: FFT Is Not The Bottleneck

Phenix's own documentation notes that FFT is <50% of refinement cycle runtime. The rest is minimization, geometry restraints, scaling, and overhead. Our GPU pipeline eliminates the FFT bottleneck, but the scaling/σA steps stay on CPU (~0.4ms) and the argmin minimizer step adds per-iteration overhead. The per-cycle budget is dominated by the density splatting + 3 FFTs + gradient sampling, with dispatch overhead as a fixed tax.

### Worst Case: Large Ribosomal Subunit (~3000 Residues)

Unit cell dimensions up to ~200×400×600 Å. At 3× oversampling at 2.5 Å: grid ~240×480×720 ≈ 83M points.

Each complex grid (2×f32 per point) = ~640MB. Full pipeline needs ~4-5GB simultaneously (density grid, mask grid, Fc, Fmask, Ftotal, map coefficients, gradient map).

GPU FFT at this grid size: ~30-60ms per transform, ~120-240ms for the 4 FFTs per refinement cycle. Full convergence (10 cycles): ~1.5-4 seconds on GPU vs ~30-60 seconds on CPU.

### Design Target

Architect for grids up to **512×512×1024** (~268M points). This covers every ribosome and large complex ever solved by X-ray crystallography.

The largest structure ever solved by X-ray crystallography is the Bluetongue virus core (PDB 2BTV) at 700 Å diameter, ~1000 protein components, 3.5 Å resolution. Virus capsids exploit 60-fold NCS — refinement programs work with the asymmetric unit, not the full virus. For the full unit cell grid of a virus crystal (potentially 1000³), fall back to CPU or require 8+ GB VRAM.

### Buffer Size Considerations

wgpu's default `max_buffer_size` is 256MB, `max_storage_buffer_binding_size` is 128MB. These are defaults, not hard caps. Desktop GPUs typically support 1-4GB per buffer when elevated limits are requested:

```rust
let limits = wgpu::Limits {
    max_buffer_size: 1 << 30,                  // 1 GB
    max_storage_buffer_binding_size: 1 << 30,   // 1 GB
    ..wgpu::Limits::default()
};
```

For the common 1000-residue case, each buffer is ~10-20MB — well within default limits. Elevated limits only needed for ribosome-scale structures. Check adapter limits at startup and gracefully degrade.

Use **multiple separate buffers** (density grid, mask grid, Fc grid, gradient map) rather than one monolithic allocation, keeping each individual buffer under the adapter's reported limits.

---

## 7. Build Order

### Phase 1: CPU Implementation (Get the Math Right)

1. **`molex::xtal::form_factors`** — ITC 5-Gaussian scattering factor tables. Pure data, no computation. These are published constants from the International Tables for Crystallography.

2. **`molex::xtal::density`** — Atomic density splatting onto a 3D grid. For each atom: evaluate the Gaussian form factor at each nearby grid point (typically within a cutoff of ~2.5 Å / d_min), accounting for the atom's B-factor and occupancy. CPU version is a triple-nested loop.

3. **`molex::xtal::fft_cpu`** — 3D FFT wrapper using `rustfft`. Compose from 1D transforms along each axis. This is the reference implementation for correctness testing.

4. **`molex::xtal::solvent_mask`** — Jiang & Brünger bulk solvent mask. For each grid point: check if it's within rprobe of any atom (outside = solvent). Then shrink the mask by rshrink. CPU version: distance checks per grid point.

5. **`molex::xtal::scaling`** — Least-squares fitting of k_overall, B_aniso (6-parameter anisotropic), k_sol, B_sol against Fobs. This is a small optimization problem (~10 parameters) — Gauss-Newton or Levenberg-Marquardt, runs on CPU even in the GPU path because it's tiny.

6. **`molex::xtal::sigma_a`** — Estimate σA parameters (D and σ) in ~15 resolution bins from the agreement between Fobs and Fc, using only the free reflections (R-free set). This is binned statistics, runs on CPU.

7. **`molex::xtal::map_coefficients`** — Compute weighted map coefficients: 2mFo-DFc = 2m·|Fobs|·exp(iφc) - D·Fc, where m = I₁(x)/I₀(x) and x = 2|Fobs|·D·|Fc|/(ε·σ²). Per-reflection arithmetic.

8. **`molex::xtal::targets`** — ML target function value and gradients. The gradient with respect to B-factors is sampled from the inverse FFT of the gradient map coefficients at each atom position.

**Validation**: Compare output maps against Refmac/Phenix output for the same structure + data. Use a well-characterized test case (e.g., lysozyme, PDB 1LYZ or similar with deposited structure factors).

### Phase 2: GPU Acceleration

9. **CubeCL FFT butterfly kernel** — Implement Cooley-Tukey radix-2/4 butterfly in CubeCL, using gpu-fft's source (github.com/eugenehp/gpu-fft) as algorithmic reference. This is ~100-150 lines of kernel code. Key design: operates on CubeCL tensor handles, NOT `Vec<f32>` — data stays GPU-resident. Compose 3D FFT from three batched 1D passes with CubeCL transposition kernels between passes. **This is the highest-risk implementation task and should be spiked first.**

10. **CubeCL density splatting kernel** — One workgroup per atom (or per atom neighborhood), atomic adds to grid buffer. The main challenge is handling the atomic float addition (WGSL doesn't have native atomic float add — use atomic compare-exchange loop or partition the grid).

11. **CubeCL solvent mask kernel** — One thread per grid point, distance check against atom positions. Straightforward parallel map.

12. **CubeCL ML weights + map coefficients** — Per-reflection arithmetic, trivial parallel map. Bessel function I₁(x)/I₀(x) approximated by rational polynomial.

13. **CubeCL B-factor gradient kernel** — Sample the gradient map at each atom position (trilinear interpolation), multiply by radial derivative.

### Phase 3: Integration with viso

14. **Density grid → marching cubes handoff** — The density grid buffer from molex's GPU pipeline feeds directly into viso's marching cubes classify pass. This is the key zero-copy integration point. Verify CubeCL buffer → wgpu buffer interop works.

15. **GPU marching cubes migration** — Move viso's current CPU marching cubes to a GPU compute shader pipeline: classify+count → prefix sum → stream compaction → vertex generation. Reference: Will Usher's WebGPU marching cubes (github.com/Twinklebear/webgpu-marching-cubes). The lookup tables (EDGE_TABLE, TRI_TABLE) move to a GPU storage buffer.

### Phase 4: Marching Cubes as Reusable Primitive

The GPU prefix sum and stream compaction primitives built for marching cubes are reusable across:
- Molecular surface (SES) construction
- Volumetric electron density ray-marching
- Real-time electrostatics (APBS-style)
- Fluid solvent visualization (MLS-MPM)

---

## 8. API Design

### molex Public API (CPU Path)

```rust
// molex::xtal

pub struct ReflectionData {
    pub h: i32, pub k: i32, pub l: i32,
    pub f_obs: f32,
    pub sigma_f: f32,
    pub free_flag: bool,
}

pub struct DensityRefinement {
    pub unit_cell: UnitCell,
    pub spacegroup: SpaceGroup,
    pub f_obs: Vec<ReflectionData>,
    pub grid_dims: [usize; 3],
    // Refined parameters
    pub b_factors: Vec<f32>,
    pub k_sol: f32,
    pub b_sol: f32,
    pub k_overall: f32,
    pub b_aniso: [f32; 6],
    pub sigma_a_d: Vec<f32>,      // per-resolution-bin D values
    pub sigma_a_sigma: Vec<f32>,  // per-resolution-bin σ values
}

impl DensityRefinement {
    /// Full pipeline: model → Fc → mask → scaling → map coefficients → density grid
    pub fn compute_map(&self, atoms: &[Atom]) -> DensityGrid;

    /// One cycle of ML B-factor refinement, returns updated B-factors
    pub fn refine_b_factors(&mut self, atoms: &[Atom]) -> &[f32];

    /// Fit bulk solvent + overall scale against Fobs
    pub fn fit_scaling(&mut self, atoms: &[Atom]);

    /// Update σA estimates from free reflections
    pub fn update_sigma_a(&mut self, atoms: &[Atom]);
}

pub struct DensityGrid {
    pub data: Vec<f32>,
    pub dims: [usize; 3],
    pub unit_cell: UnitCell,
}
```

### molex GPU API

```rust
// molex::xtal::gpu (behind feature = "xtal-gpu")

pub struct GpuDensityRefinement {
    device: Arc<wgpu::Device>,
    queue: Arc<wgpu::Queue>,
    // GPU buffers
    atom_buffer: wgpu::Buffer,
    grid_buffer: wgpu::Buffer,
    fobs_buffer: wgpu::Buffer,
    // CubeCL runtime
    cubecl_device: WgpuDevice,
    // ... pipelines, etc.
}

impl GpuDensityRefinement {
    pub fn new(device: Arc<wgpu::Device>, queue: Arc<wgpu::Queue>, ...) -> Self;

    /// Same pipeline, all on GPU, returns handle to GPU-resident density grid
    /// that can go directly into marching cubes with no readback
    pub fn compute_map(&self, atom_buffer: &wgpu::Buffer) -> &wgpu::Buffer;

    /// One ML refinement cycle on GPU
    pub fn refine_cycle(&mut self);
}
```

### Handoff to viso

```rust
// In viso's engine, after molex computes the density:
let density_buf = gpu_refinement.compute_map(&atoms_buf);
self.marching_cubes.extract(&density_buf, iso_level);
// density_buf is a wgpu::Buffer on the same device — zero copy
```

---

## 9. Minimization Strategy

### Use argmin — Don't Roll Your Own

Use the `argmin` crate (argmin-rs.org) for all optimization. It's a pure-Rust optimization framework where you implement `CostFunction` and `Gradient` traits for your problem, then plug in any solver. The problem definition stays the same regardless of which algorithm you use.

```rust
use argmin::prelude::*;
use argmin::solver::gradientdescent::SteepestDescent;
use argmin::solver::quasinewton::LBFGS;
use argmin::solver::linesearch::MoreThuenteLineSearch;

struct BFactorRefinement {
    atoms: Vec<Atom>,
    f_obs: Vec<ReflectionData>,
    unit_cell: UnitCell,
    spacegroup: SpaceGroup,
    // ... all the data needed to compute the ML target
}

impl CostFunction for BFactorRefinement {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, b_factors: &Vec<f64>) -> Result<f64, Error> {
        // Run the full pipeline: splat → FFT → mask → scale → map coeffs
        // Return the ML target function value (negative log-likelihood)
        // PLUS restraint terms (bonded B smoothness, Wilson B prior)
    }
}

impl Gradient for BFactorRefinement {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;

    fn gradient(&self, b_factors: &Vec<f64>) -> Result<Vec<f64>, Error> {
        // Run the pipeline, compute mFo-DFc gradient map,
        // sample at atom positions, add restraint gradients
        // Return ∂L/∂B for each atom
    }
}

// Start with gradient descent:
let linesearch = MoreThuenteLineSearch::new();
let solver = SteepestDescent::new(linesearch);

// Upgrade to LBFGS when needed — same problem, one line change:
// let solver = LBFGS::new(linesearch, 7);

let result = Executor::new(problem, solver)
    .configure(|state| state.param(initial_b_factors).max_iters(15))
    .add_observer(SlogLogger::term(), ObserverMode::Every(1))  // free logging
    .run()?;
```

### Why argmin Instead of Hand-Written Gradient Descent

- **Step size / line search**: argmin handles this. More-Thuente line search is robust and well-tested. Hand-rolling step size selection is error-prone and was likely part of the difficulty with past Rosetta minimization attempts.
- **Convergence monitoring**: argmin tracks cost function value, gradient norm, iteration count, and supports custom termination conditions. You get R-free monitoring for free by implementing it as a termination criterion or observer.
- **Solver swapping**: Start with `SteepestDescent`, upgrade to `LBFGS` or `ConjugateGradient` by changing one line. The `CostFunction` and `Gradient` implementations don't change.
- **Checkpointing**: argmin can checkpoint and resume optimization. Useful for long refinement runs.
- **No reinventing wheels**: Line search algorithms, convergence criteria, and Hessian approximations are tricky to get right. Let a tested library handle them.

### Complications Beyond "Just Minimize B-factors"

The minimization is not as simple as a single call to argmin. The full refinement loop has interleaved parameter updates:

```
for macro_cycle in 0..5:
    1. Compute Fc from current model (splat + FFT)
    2. Compute Fmask (solvent mask + FFT)
    3. Fit scaling parameters (k_overall, B_aniso, k_sol, B_sol) — this is a SEPARATE
       small optimization (Gauss-Newton or Levenberg-Marquardt on ~10 parameters,
       also via argmin)
    4. Update σA estimates from free reflections (binned statistics, not optimization)
    5. Run argmin on B-factors for N iterations (inner loop)
       - Each cost/gradient evaluation runs: splat → FFT → combine → map coeffs → sample
       - Restraint terms added to cost and gradient
    6. Check R-free for convergence / overfitting
```

Steps 3 and 4 update global parameters that change the landscape for step 5. This is why it's macro-cycles wrapping micro-cycles, not one big optimization. argmin handles the inner loop (step 5). The outer loop (macro-cycle) is your orchestration code.

### Restraints

B-factor restraints are additional terms added to BOTH the cost function and the gradient inside the argmin problem:

- **Bonded smoothness**: For each bonded atom pair (i,j): cost += w_bond × (B_i - B_j)², gradient_i += 2 × w_bond × (B_i - B_j)
- **Wilson B prior**: For each atom: cost += w_wilson × (B_i - B_wilson)², gradient_i += 2 × w_wilson × (B_i - B_wilson). B_wilson is the overall Wilson B-factor estimated from the data.
- **Bounds**: B-factors should stay in [2, 300] Å². Use argmin's L-BFGS-B (bounded LBFGS) or clamp after each step.

The restraint weights (w_bond, w_wilson) are hyperparameters. Refmac uses defaults that work well for most structures. Start with Refmac's defaults.

### Future: Coordinate Refinement

If you later want to refine atomic coordinates (xyz), the same framework extends naturally:
- Parameter vector becomes [x1,y1,z1, B1, x2,y2,z2, B2, ...]
- Gradient includes ∂L/∂x, ∂L/∂y, ∂L/∂z per atom (sampled from the gradient map, chain-ruled through the density function)
- Add geometry restraints (bond lengths, angles, torsions, Ramachandran, clashes) to cost and gradient
- LBFGS becomes essential at this point — the parameter space is 4N and the landscape is more complex

This is where argmin really pays off. The coordinate refinement problem is exactly what LBFGS was designed for, and you don't have to rewrite any minimization machinery.

---

## 10. Key Dependencies

| Crate | Purpose | License |
|---|---|---|
| `cubecl-wgpu` | GPU compute kernel framework (FFT, splatting, mask, ML kernels) | MIT/Apache-2.0 |
| `rustfft` | CPU FFT fallback | MIT/Apache-2.0 |
| `argmin` | Optimization framework (gradient descent, LBFGS, CG, etc.) | MIT/Apache-2.0 |
| `wgpu` | GPU API (already in viso) | MIT/Apache-2.0 |

### NOT Dependencies

- **gpu-fft**: Evaluated and rejected. API forces host round-trips that blow the <10ms budget. Butterfly kernel logic used as algorithmic reference for our custom CubeCL FFT implementation.
- **gemmi**: Used as algorithmic reference only. We reimplement in Rust, do not link to C++.
- **VkFFT**: Considered and rejected. Would require breaking out of wgpu's abstraction to access raw Vulkan/Metal handles.
- **CubeCL for viso**: CubeCL is compute-only, can't express render passes. Viso stays on raw wgpu + WGSL for its render pipeline.
- **liblbfgs / hand-rolled minimizers**: Use argmin instead. It provides the same algorithms with a consistent, tested framework.

---

## 11. Open Questions / Risks

1. **CubeCL ↔ wgpu buffer interop**: Can CubeCL buffer handles be unwrapped to raw `wgpu::Buffer` for use in viso's bind groups? The `WgpuDevice::Existing` design implies yes, but needs verification. **This is the first thing to spike.**

2. **CubeCL FFT kernel correctness**: Writing a GPU FFT from scratch is the highest-risk task. Use gpu-fft's butterfly logic (github.com/eugenehp/gpu-fft) as reference. Validate against rustfft output for identical inputs before integrating into the pipeline.

3. **CubeCL dispatch overhead**: ~65μs per kernel dispatch. With ~8 kernel dispatches per cycle, that's ~0.5ms of fixed overhead. Acceptable within the 10ms budget but leaves less room for the actual compute. Monitor this — if CubeCL's Metal backend has the known launch overhead issues on macOS, it could be worse.

4. **wgpu efficiency gap**: WebGPU compute reaches ~17-25% of theoretical GPU peak on Apple Silicon, vs ~75% for native CUDA. The 2-4x pessimism factor in our estimates accounts for this. On NVIDIA hardware with CubeCL's CUDA backend, performance will be significantly better. Test on both platforms.

5. **Atomic float addition in WGSL/CubeCL**: The density splatting kernel needs multiple atoms to contribute to the same grid point. WGSL has `atomicAdd` for integers but not floats. Options: (a) fixed-point accumulation (multiply by 1000, use integer atomics, divide at end), (b) atomic compare-exchange loop on f32 bit patterns, (c) privatize grid segments per workgroup and merge. Option (a) is simplest and proven in practice.

6. **Mixed-radix FFT**: Initial implementation is power-of-two only. Crystallographic grids often have dimensions like 135, 150 (products of 2, 3, 5). Bias `findSampling` to prefer powers of two on GPU path. Padding 150→256 wastes ~70% compute on that axis. Mixed-radix (3, 5) kernels are a later optimization.

7. **f32 precision**: All GPU computation is single precision. Crystallographic maps don't need f64 (Refmac and Phenix compute maps in single precision). The scaling/fitting step (least-squares for k_sol, B_sol) might want f64 for numerical stability of the normal equations — do this step on CPU regardless since it's tiny (~10 parameters).

---

## 12. Key References

- **Gemmi** (algorithmic reference for SF calc + bulk solvent + scaling): gemmi.readthedocs.io/en/latest/scattering.html
- **gpu-fft** (algorithmic reference for CubeCL FFT butterfly kernels): github.com/eugenehp/gpu-fft
- **Will Usher WebGPU Marching Cubes**: willusher.io/graphics/2024/04/22/webgpu-marching-cubes/
- **CubeCL**: github.com/tracel-ai/cubecl
- **CubeCL WgpuDevice::Existing**: docs.rs/cubecl-wgpu/latest/cubecl_wgpu/enum.WgpuDevice.html
- **argmin optimization framework**: argmin-rs.org
- **SFCalculator** (GPU SF calc in PyTorch/JAX, validates approach): Hekstra Lab, 2025
- **ISOLDE** (real-time crystallographic refinement, our performance target): ~550ms/update on CPU
- **Amber xray module** (proof ML target works on GPU): doi.org/10.1021/acs.jcim.3c01531
- **VkFFT** (reference for FFT algorithm design): github.com/DTolm/VkFFT
- **Murshudov 1997** (ML refinement in Refmac): pubmed.ncbi.nlm.nih.gov/15299926/
- **Pannu & Read 1996** (ML target function)
- **Brünger 1997** (bulk solvent + scaling)
- **Read 1986** (σA weighting for map coefficients)
