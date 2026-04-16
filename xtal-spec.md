# xtal Module — Complete Implementation Specification

All algorithms, formulas, data tables, and pseudocode needed to implement `molex::xtal`.
Output of 12 parallel research agents reading gemmi, Clipper, argmin, CCTBX, and
crystallographic references. Each section maps to one source file.

---

## 1. `types.rs` — Core Crystallographic Types

### Constants

```rust
pub const DEN: i32 = 24;
```

### Symop

Integer arithmetic with DEN=24. Represents all crystallographic fractional translations
exactly (1/2=12, 1/3=8, 2/3=16, 1/4=6, 3/4=18, 1/6=4, 5/6=20, 1/8=3).

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Symop {
    pub rot: [[i32; 3]; 3],  // scaled by DEN. Identity = [[24,0,0],[0,24,0],[0,0,24]]
    pub tran: [i32; 3],      // each element in [0, DEN)
}
```

**Apply to fractional coords:**
```
x' = (R[0][0]*x + R[0][1]*y + R[0][2]*z) / DEN + t[0] / DEN
y' = (R[1][0]*x + R[1][1]*y + R[1][2]*z) / DEN + t[1] / DEN
z' = (R[2][0]*x + R[2][1]*y + R[2][2]*z) / DEN + t[2] / DEN
```

**Apply to Miller indices (h,k,l) — uses TRANSPOSE:**
```
h' = (R[0][0]*h + R[1][0]*k + R[2][0]*l) / DEN
k' = (R[0][1]*h + R[1][1]*k + R[2][1]*l) / DEN
l' = (R[0][2]*h + R[1][2]*k + R[2][2]*l) / DEN
```

**Phase shift:** `phase = -2*pi/DEN * (h*t[0] + k*t[1] + l*t[2])`

**Combining A*B:**
```
C.rot[i][j] = (A.rot[i][0]*B.rot[0][j] + A.rot[i][1]*B.rot[1][j] + A.rot[i][2]*B.rot[2][j]) / DEN
C.tran[i] = (A.tran[i]*DEN + A.rot[i][0]*B.tran[0] + A.rot[i][1]*B.tran[1] + A.rot[i][2]*B.tran[2]) / DEN
```
Then wrap translations: `t = ((t % DEN) + DEN) % DEN`.

### GroupOps

```rust
pub struct GroupOps {
    pub sym_ops: Vec<Symop>,
    pub cen_ops: Vec<[i32; 3]>,
}
```
Total order = `sym_ops.len() * cen_ops.len()`.

**Centering vectors (DEN-scaled):**
- P: `[[0,0,0]]`
- C: `[[0,0,0], [12,12,0]]`
- I: `[[0,0,0], [12,12,12]]`
- F: `[[0,0,0], [0,12,12], [12,0,12], [12,12,0]]`
- R (hex): `[[0,0,0], [16,8,8], [8,16,16]]`

### Epsilon, Centric, Absences

```rust
fn epsilon_factor(hkl: [i32; 3], group: &GroupOps) -> i32 {
    let denh = [DEN * hkl[0], DEN * hkl[1], DEN * hkl[2]];
    let mut eps = 0;
    for op in &group.sym_ops {
        if op.apply_to_hkl_without_division(hkl) == denh { eps += 1; }
    }
    eps * group.cen_ops.len() as i32
}

fn is_centric(hkl: [i32; 3], ops: &[Symop]) -> bool {
    let neg = [-DEN * hkl[0], -DEN * hkl[1], -DEN * hkl[2]];
    ops.iter().any(|op| op.apply_to_hkl_without_division(hkl) == neg)
}

fn is_systematically_absent(hkl: [i32; 3], group: &GroupOps) -> bool {
    for cen in &group.cen_ops[1..] {
        if (hkl[0]*cen[0] + hkl[1]*cen[1] + hkl[2]*cen[2]) % DEN != 0 { return true; }
    }
    let denh = [DEN*hkl[0], DEN*hkl[1], DEN*hkl[2]];
    for op in &group.sym_ops[1..] {
        if op.apply_to_hkl_without_division(hkl) == denh {
            for cen in &group.cen_ops {
                let t = [op.tran[0]+cen[0], op.tran[1]+cen[1], op.tran[2]+cen[2]];
                if (hkl[0]*t[0] + hkl[1]*t[1] + hkl[2]*t[2]) % DEN != 0 { return true; }
            }
        }
    }
    false
}
```

### Top 10 Space Group Symops (hardcoded, DEN=24)

Convention: `rot[row][col]`, `tran[xyz]`. Source: CCTBX, verified against ITA.

```rust
pub type SymopData = ([[i32; 3]; 3], [i32; 3]);

// ---- P 1 (SG 1) ---- 1 op, P centering
pub const SG1_CEN: &[[i32; 3]] = &[[0, 0, 0]];
pub const SG1_OPS: &[SymopData] = &[
    ([[24,0,0],[0,24,0],[0,0,24]], [0,0,0]),   // x, y, z
];

// ---- P 1 21 1 (SG 4) ---- 2 ops, P centering
pub const SG4_CEN: &[[i32; 3]] = &[[0, 0, 0]];
pub const SG4_OPS: &[SymopData] = &[
    ([[24,0,0],[0,24,0],[0,0,24]], [0,0,0]),     // x, y, z
    ([[-24,0,0],[0,24,0],[0,0,-24]], [0,12,0]),   // -x, y+1/2, -z
];

// ---- C 1 2 1 (SG 5) ---- 2 ops, C centering
pub const SG5_CEN: &[[i32; 3]] = &[[0,0,0], [12,12,0]];
pub const SG5_OPS: &[SymopData] = &[
    ([[24,0,0],[0,24,0],[0,0,24]], [0,0,0]),     // x, y, z
    ([[-24,0,0],[0,24,0],[0,0,-24]], [0,0,0]),   // -x, y, -z
];

// ---- P 21 21 21 (SG 19) ---- 4 ops, P centering
pub const SG19_CEN: &[[i32; 3]] = &[[0, 0, 0]];
pub const SG19_OPS: &[SymopData] = &[
    ([[24,0,0],[0,24,0],[0,0,24]], [0,0,0]),       // x, y, z
    ([[-24,0,0],[0,-24,0],[0,0,24]], [12,0,12]),   // -x+1/2, -y, z+1/2
    ([[24,0,0],[0,-24,0],[0,0,-24]], [12,12,0]),   // x+1/2, -y+1/2, -z
    ([[-24,0,0],[0,24,0],[0,0,-24]], [0,12,12]),   // -x, y+1/2, -z+1/2
];

// ---- I 2 2 2 (SG 23) ---- 4 ops, I centering
pub const SG23_CEN: &[[i32; 3]] = &[[0,0,0], [12,12,12]];
pub const SG23_OPS: &[SymopData] = &[
    ([[24,0,0],[0,24,0],[0,0,24]], [0,0,0]),     // x, y, z
    ([[-24,0,0],[0,-24,0],[0,0,24]], [0,0,0]),   // -x, -y, z
    ([[24,0,0],[0,-24,0],[0,0,-24]], [0,0,0]),   // x, -y, -z
    ([[-24,0,0],[0,24,0],[0,0,-24]], [0,0,0]),   // -x, y, -z
];

// ---- P 41 21 2 (SG 92) ---- 8 ops, P centering
pub const SG92_CEN: &[[i32; 3]] = &[[0, 0, 0]];
pub const SG92_OPS: &[SymopData] = &[
    ([[24,0,0],[0,24,0],[0,0,24]], [0,0,0]),        // x, y, z
    ([[0,-24,0],[24,0,0],[0,0,24]], [12,12,6]),     // -y+1/2, x+1/2, z+1/4
    ([[-24,0,0],[0,-24,0],[0,0,24]], [0,0,12]),     // -x, -y, z+1/2
    ([[0,24,0],[-24,0,0],[0,0,24]], [12,12,18]),    // y+1/2, -x+1/2, z+3/4
    ([[24,0,0],[0,-24,0],[0,0,-24]], [12,12,18]),   // x+1/2, -y+1/2, -z+3/4
    ([[0,24,0],[24,0,0],[0,0,-24]], [0,0,0]),       // y, x, -z
    ([[-24,0,0],[0,24,0],[0,0,-24]], [12,12,6]),    // -x+1/2, y+1/2, -z+1/4
    ([[0,-24,0],[-24,0,0],[0,0,-24]], [0,0,12]),    // -y, -x, -z+1/2
];

// ---- P 43 21 2 (SG 96) ---- 8 ops, P centering
pub const SG96_CEN: &[[i32; 3]] = &[[0, 0, 0]];
pub const SG96_OPS: &[SymopData] = &[
    ([[24,0,0],[0,24,0],[0,0,24]], [0,0,0]),        // x, y, z
    ([[0,-24,0],[24,0,0],[0,0,24]], [12,12,18]),    // -y+1/2, x+1/2, z+3/4
    ([[-24,0,0],[0,-24,0],[0,0,24]], [0,0,12]),     // -x, -y, z+1/2
    ([[0,24,0],[-24,0,0],[0,0,24]], [12,12,6]),     // y+1/2, -x+1/2, z+1/4
    ([[24,0,0],[0,-24,0],[0,0,-24]], [12,12,6]),    // x+1/2, -y+1/2, -z+1/4
    ([[0,24,0],[24,0,0],[0,0,-24]], [0,0,0]),       // y, x, -z
    ([[-24,0,0],[0,24,0],[0,0,-24]], [12,12,18]),   // -x+1/2, y+1/2, -z+3/4
    ([[0,-24,0],[-24,0,0],[0,0,-24]], [0,0,12]),    // -y, -x, -z+1/2
];

// ---- P 31 2 1 (SG 152) ---- 6 ops, P centering
pub const SG152_CEN: &[[i32; 3]] = &[[0, 0, 0]];
pub const SG152_OPS: &[SymopData] = &[
    ([[24,0,0],[0,24,0],[0,0,24]], [0,0,0]),         // x, y, z
    ([[0,-24,0],[24,-24,0],[0,0,24]], [0,0,8]),      // -y, x-y, z+1/3
    ([[-24,24,0],[-24,0,0],[0,0,24]], [0,0,16]),     // -x+y, -x, z+2/3
    ([[0,24,0],[24,0,0],[0,0,-24]], [0,0,0]),        // y, x, -z
    ([[-24,0,0],[-24,24,0],[0,0,-24]], [0,0,8]),     // -x, -x+y, -z+1/3
    ([[24,-24,0],[0,-24,0],[0,0,-24]], [0,0,16]),    // x-y, -y, -z+2/3
];

// ---- P 32 2 1 (SG 154) ---- 6 ops, P centering
pub const SG154_CEN: &[[i32; 3]] = &[[0, 0, 0]];
pub const SG154_OPS: &[SymopData] = &[
    ([[24,0,0],[0,24,0],[0,0,24]], [0,0,0]),         // x, y, z
    ([[0,-24,0],[24,-24,0],[0,0,24]], [0,0,16]),     // -y, x-y, z+2/3
    ([[-24,24,0],[-24,0,0],[0,0,24]], [0,0,8]),      // -x+y, -x, z+1/3
    ([[0,24,0],[24,0,0],[0,0,-24]], [0,0,0]),        // y, x, -z
    ([[-24,0,0],[-24,24,0],[0,0,-24]], [0,0,16]),    // -x, -x+y, -z+2/3
    ([[24,-24,0],[0,-24,0],[0,0,-24]], [0,0,8]),     // x-y, -y, -z+1/3
];

// ---- P 61 2 2 (SG 178) ---- 12 ops, P centering
pub const SG178_CEN: &[[i32; 3]] = &[[0, 0, 0]];
pub const SG178_OPS: &[SymopData] = &[
    ([[24,0,0],[0,24,0],[0,0,24]], [0,0,0]),         // x, y, z
    ([[24,-24,0],[24,0,0],[0,0,24]], [0,0,4]),       // x-y, x, z+1/6
    ([[0,-24,0],[24,-24,0],[0,0,24]], [0,0,8]),      // -y, x-y, z+1/3
    ([[-24,0,0],[0,-24,0],[0,0,24]], [0,0,12]),      // -x, -y, z+1/2
    ([[-24,24,0],[-24,0,0],[0,0,24]], [0,0,16]),     // -x+y, -x, z+2/3
    ([[0,24,0],[-24,24,0],[0,0,24]], [0,0,20]),      // y, -x+y, z+5/6
    ([[0,-24,0],[-24,0,0],[0,0,-24]], [0,0,20]),     // -y, -x, -z+5/6
    ([[24,-24,0],[0,-24,0],[0,0,-24]], [0,0,0]),     // x-y, -y, -z
    ([[24,0,0],[24,-24,0],[0,0,-24]], [0,0,4]),      // x, x-y, -z+1/6
    ([[0,24,0],[24,0,0],[0,0,-24]], [0,0,8]),        // y, x, -z+1/3
    ([[-24,24,0],[0,24,0],[0,0,-24]], [0,0,12]),     // -x+y, y, -z+1/2
    ([[-24,0,0],[-24,24,0],[0,0,-24]], [0,0,16]),    // -x, -x+y, -z+2/3
];
```

### SpaceGroup

```rust
pub struct SpaceGroup {
    pub number: u16,
    pub hm: &'static str,
    pub ops: GroupOps,
    pub crystal_system: CrystalSystem,
}

pub enum CrystalSystem { Triclinic, Monoclinic, Orthorhombic, Tetragonal, Trigonal, Hexagonal, Cubic }
```

Lookup by number (match or HashMap over 10 entries).

### Grid Factors

```rust
fn grid_factors(sg_number: u16) -> [usize; 3] {
    match sg_number {
        1   => [1, 1, 1],  // P 1
        4   => [1, 2, 1],  // P 21: screw along b
        5   => [2, 2, 1],  // C 2: C centering (1/2,1/2,0)
        19  => [2, 2, 2],  // P 212121: three 21 screws
        23  => [2, 2, 2],  // I 222: I centering (1/2,1/2,1/2)
        92  => [2, 2, 4],  // P 41212: 41 screw (1/4 along c)
        96  => [2, 2, 4],  // P 43212: 43 screw (1/4 along c)
        152 => [1, 1, 3],  // P 3121: 31 screw (1/3 along c)
        154 => [1, 1, 3],  // P 3221: 32 screw (1/3 along c)
        178 => [1, 1, 6],  // P 6122: 61 screw (1/6 along c)
        _   => panic!("unsupported space group"),
    }
}

fn requires_equal_uv(sg_number: u16) -> bool {
    matches!(sg_number, 92 | 96 | 152 | 154 | 178)  // tetragonal, trigonal, hexagonal
}
```

### Grid Dimension Selection

```rust
fn has_small_factorization(mut n: i32) -> bool {
    while n % 2 == 0 { n /= 2; }
    while n % 3 == 0 { n /= 3; }
    while n % 5 == 0 { n /= 5; }
    n == 1
}

fn round_up_to_smooth(exact: f64) -> usize {
    let mut n = exact.ceil() as i32;
    while !has_small_factorization(n) { n += 1; }
    n as usize
}
```

For GPU path: round to nearest power of 2 instead.

### UnitCell

```rust
pub struct UnitCell {
    pub a: f64, pub b: f64, pub c: f64,
    pub alpha: f64, pub beta: f64, pub gamma: f64,
    // Precomputed:
    pub volume: f64,
    pub ar: f64, pub br: f64, pub cr: f64,
    pub cos_alphar: f64, pub cos_betar: f64, pub cos_gammar: f64,
    pub orth: [[f64; 3]; 3],
    pub frac: [[f64; 3]; 3],
}
```

**calculate_properties():**
```
1. Trig (exact for 90°): cos_a = (alpha==90)?0:cos(rad), sin_a = (alpha==90)?1:sin(rad), etc.
2. volume = a*b*c*sqrt(1 - cos_a² - cos_b² - cos_g² + 2*cos_a*cos_b*cos_g)
3. ar = b*c*sin_a/V,  br = a*c*sin_b/V,  cr = a*b*sin_g/V
4. cas_b = (cos_b*cos_g - cos_a)/sin_g    // cos(alpha*)*sin(beta)
   cos_alphar = cas_b / sin_b
   cos_betar  = (cos_a*cos_g - cos_b)/(sin_a*sin_g)
   cos_gammar = (cos_a*cos_b - cos_g)/(sin_a*sin_b)
   sin_alphar = sqrt(1 - cos_alphar²)
5. orth[0] = [a, b*cos_g, c*cos_b]
   orth[1] = [0, b*sin_g, -c*cas_b]
   orth[2] = [0, 0, c*sin_b*sin_alphar]    // = V/(a*b*sin_g)
6. frac[0] = [1/a, -cos_g/(sin_g*a), -(cos_g*cas_b + cos_b*sin_g)/(sin_alphar*sin_b*sin_g*a)]
   frac[1] = [0, 1/(b*sin_g), cos_alphar/(sin_alphar*sin_g*b)]
   frac[2] = [0, 0, 1/(c*sin_b*sin_alphar)]
```

**d_star_sq(h,k,l) = 1/d²:**
```
h²*ar² + k²*br² + l²*cr² + 2hk*ar*br*cos_gammar + 2hl*ar*cr*cos_betar + 2kl*br*cr*cos_alphar
```

### Reflection, DensityGrid

```rust
pub struct Reflection {
    pub h: i32, pub k: i32, pub l: i32,
    pub f_obs: f32, pub sigma_f: f32, pub free_flag: bool,
}

pub struct DensityGrid {
    pub data: Vec<f32>,
    pub nu: usize, pub nv: usize, pub nw: usize,
    pub unit_cell: UnitCell,
}
```

---

## 2. `form_factors.rs` — IT92 Scattering Factors + vdW Radii

IT92 (Cromer-Mann) 4-Gaussian + constant. Source: gemmi `it92.hpp`.

```rust
pub struct FormFactor {
    pub a: [f64; 4],
    pub b: [f64; 4],
    pub c: f64,
}
```

`f(s) = c + Σ a_i * exp(-b_i * s²)` where `s = sin(θ)/λ`.

**Coefficients (a1,b1,a2,b2,a3,b3,a4,b4,c):**

| Element | a1 | b1 | a2 | b2 | a3 | b3 | a4 | b4 | c |
|---------|----|----|----|----|----|----|----|----|---|
| H | 0.493002 | 10.5109 | 0.322912 | 26.1257 | 0.140191 | 3.14236 | 0.04081 | 57.7997 | 0.003038 |
| C | 2.31 | 20.8439 | 1.02 | 10.2075 | 1.5886 | 0.5687 | 0.865 | 51.6512 | 0.2156 |
| N | 12.2126 | 0.0057 | 3.1322 | 9.8933 | 2.0125 | 28.9975 | 1.1663 | 0.5826 | -11.529 |
| O | 3.0485 | 13.2771 | 2.2868 | 5.7011 | 1.5463 | 0.3239 | 0.867 | 32.9089 | 0.2508 |
| S | 6.9053 | 1.4679 | 5.2034 | 22.2151 | 1.4379 | 0.2536 | 1.5863 | 56.172 | 0.8669 |
| P | 6.4345 | 1.9067 | 4.1791 | 27.157 | 1.78 | 0.526 | 1.4908 | 68.1645 | 1.1149 |
| Fe | 11.7695 | 4.7611 | 7.3573 | 0.3072 | 3.5222 | 15.3535 | 2.3045 | 76.8805 | 1.0369 |
| Zn | 14.0743 | 3.2655 | 7.0318 | 0.2333 | 5.1652 | 10.3163 | 2.41 | 58.7097 | 1.3041 |
| Ca | 8.6266 | 10.4421 | 7.3873 | 0.6599 | 1.5899 | 85.7484 | 1.0211 | 178.437 | 1.3751 |
| Mg | 5.4204 | 2.8275 | 2.1735 | 79.2611 | 1.2269 | 0.3808 | 2.3073 | 7.1937 | 0.8584 |
| Cl | 11.4604 | 0.0104 | 7.1964 | 1.1662 | 6.2556 | 18.5194 | 1.6455 | 47.7784 | -9.5574 |
| Na | 4.7626 | 3.285 | 3.1736 | 8.8422 | 1.2674 | 0.3136 | 1.1128 | 129.424 | 0.676 |
| K | 8.2186 | 12.7949 | 7.4398 | 0.7748 | 1.0519 | 213.187 | 0.8659 | 41.6841 | 1.4228 |
| Mn | 11.2819 | 5.3409 | 7.3573 | 0.3432 | 3.0193 | 17.8674 | 2.2441 | 83.7543 | 1.0896 |
| Se | 17.0006 | 2.4098 | 5.8196 | 0.2726 | 3.9731 | 15.2372 | 4.3543 | 43.8163 | 2.8409 |
| Cu | 13.338 | 3.5828 | 7.1676 | 0.247 | 5.6158 | 11.3966 | 1.6735 | 64.8126 | 1.191 |
| Co | 12.2841 | 4.2791 | 7.3409 | 0.2784 | 4.0034 | 13.5359 | 2.3488 | 71.1692 | 1.0118 |
| Ni | 12.8376 | 3.8785 | 7.292 | 0.2565 | 4.4438 | 12.1763 | 2.38 | 66.3421 | 1.0341 |

Sanity check: `sum(a) + c ≈ Z` (atomic number). E.g., C: 2.31+1.02+1.5886+0.865+0.2156=5.999 (Z=6).

**Bondi/Mantina vdW Radii (Angstroms):**

| Element | Radius |
|---------|--------|
| H | 1.20 |
| C | 1.70 |
| N | 1.55 |
| O | 1.52 |
| S | 1.80 |
| P | 1.80 |
| Fe | 1.26 |
| Zn | 1.39 |
| Ca | 2.31 |
| Mg | 1.73 |
| Cl | 1.75 |
| Na | 2.27 |
| K | 2.75 |
| Mn | 1.19 |
| Se | 1.90 |
| Cu | 1.40 |
| Co | 1.13 |
| Ni | 1.63 |

---

## 3. `bessel.rs` — Bessel Function Approximations

### I0(x) (Abramowitz & Stegun)

**|x| <= 3.75** (t = (x/3.75)²):
```
I0 = 1 + 3.5156229t + 3.0899424t² + 1.2067492t³ + 0.2659732t⁴ + 0.0360768t⁵ + 0.0045813t⁶
```

**|x| > 3.75** (t = 3.75/|x|):
```
I0 = exp(|x|)/sqrt(|x|) * (0.39894228 + 0.01328592t + 0.00225319t² - 0.00157565t³
     + 0.00916281t⁴ - 0.02057706t⁵ + 0.02635537t⁶ - 0.01647633t⁷ + 0.00392377t⁸)
```

### I1(x)/I0(x) — figure of merit "sim"

- x < 1e-6: return x/2 (Taylor)
- x > 700: return 1 - 0.5/x (asymptotic)
- Otherwise: compute I1(x) and I0(x) separately and divide
  (A&S also has I1 polynomial approximations analogous to I0)

### log(I0(x)) — "sim_integ"

- Small x: x²/4
- Large x: x - 0.5*ln(2πx)

---

## 4. `fft_cpu.rs` — 3D FFT Wrapper

Use `rustfft` crate. 3D FFT = three batched 1D passes along each axis.
`rustfft` does NOT normalize on inverse (multiply output by 1/N yourself if needed,
but for crystallography the `1/V_cell` normalization is applied at the point of use).

Grid dimension selection: see `round_up_to_smooth` and `grid_factors` in types.rs.

---

## 5. `density.rs` — Atomic Density Splatting

### Gaussian precalculation (5 terms: 4 from IT92 + constant)

```
For i in 0..4:
    t_i = 4π / (b_i + B_atom + blur)
    precal_a[i] = a_i * t_i^(3/2)
    precal_b[i] = -t_i * π

For constant (i=4):
    t_c = 4π / (B_atom + blur)
    precal_a[4] = c * t_c^(3/2)
    precal_b[4] = -t_c * π
```

Density at distance² r2: `ρ = occ * Σ precal_a[i] * exp(precal_b[i] * r2)`

### Blur

Single global constant:
```
spacing = d_min / (2 * rate)
blur = max((8π²/1.1) * spacing² - B_min, 0.0)
```

### Cutoff radius

Per-atom: `radius ≈ (8.5 + 0.075*B_eff) / (2.4 + 0.0045*B_eff)` where B_eff=B+blur.
Refine by stepping 0.5A to find where |ρ| < 1e-5.

### Splatting loop

For each ASU atom: fractionalize → find nearest grid point → iterate bounding box →
compute Cartesian distance via `orth_n = orth * diag(1/nu, 1/nv, 1/nw)` → evaluate
Gaussian sum → add to grid with periodic wrapping.

### Symmetrize (real space)

After all atoms splatted, call `symmetrize_sum`:
For each grid point orbit, sum all symmetry-equivalent values, assign sum to all points.

### Deblur (reciprocal space)

After FFT: `Fc(h,k,l) *= exp(blur * d_star_sq(h,k,l) / 4)`

---

## 6. `solvent_mask.rs` — Bulk Solvent Mask via EDT

### Pipeline
1. Initialize grid: 0=protein, 1=solvent. Mark protein where distance to any atom < r_vdw + r_probe (default r_probe=1.0A).
2. Build EDT input: solvent→0.0, protein→INF (want distance from protein points to nearest solvent).
3. Run 3D F-H EDT.
4. Threshold: protein points with distance < r_shrink (default 1.1A) become solvent.
5. FFT the mask → Fmask(h,k,l).

### Felzenszwalb-Huttenlocher 1D EDT

```
fn edt_1d(f: &[f64]) -> Vec<f64> {
    let n = f.len();
    let mut v = vec![0usize; n];      // parabola locations
    let mut z = vec![0.0f64; n + 1];  // boundaries
    z[0] = f64::NEG_INFINITY;
    z[1] = f64::INFINITY;
    let mut k = 0usize;

    for q in 1..n {
        loop {
            let vk = v[k];
            let s = (f[q] + (q*q) as f64 - f[vk] - (vk*vk) as f64)
                    / (2.0 * (q as f64 - vk as f64));
            if s > z[k] {
                k += 1; v[k] = q; z[k] = s; z[k+1] = f64::INFINITY;
                break;
            }
            k -= 1;
        }
    }

    let mut dt = vec![0.0f64; n];
    k = 0;
    for q in 0..n {
        while z[k+1] < q as f64 { k += 1; }
        let d = q as f64 - v[k] as f64;
        dt[q] = d*d + f[v[k]];
    }
    dt
}
```

### Periodic 1D EDT

Pad array by n/2+1 on each side (wrapping), run edt_1d on padded, extract center.

### 3D EDT with anisotropic spacing

Three passes. Pass 1 (along u): scale output by hu². Pass 2 (along v): divide input
by hv², run edt_1d, multiply output by hv². Pass 3 (along w): same with hw².

Spacing: `hu = |column_0(orth)| / nu`, `hv = |column_1(orth)| / nv`, `hw = |column_2(orth)| / nw`.

---

## 7. `scaling.rs` — Levenberg-Marquardt Scaling

**CORRECTION**: All parameters fit simultaneously (NOT alternating). Source: gemmi scaling.hpp.

### Model

`|F_total(h)| = k_overall * exp(-0.25 * h^T B* h) * |Fc(h) + k_sol * exp(-B_sol * stol2) * Fmask(h)|`

where `stol2 = (1/d²)/4`, `h^T B* h = h²B11 + k²B22 + l²B33 + 2hkB12 + 2hlB13 + 2klB23`.

### Parameters

`[k_overall, k_sol, B_sol, d_0..d_N]` where d_j are B* components projected onto
crystal-system constraint basis. N = 6(triclinic), 3(ortho), 2(tetra), 1(cubic).

### Initial estimates (Wilson plot)

For reflections where Fobs >= 1 and Fobs >= sigma:
`y = ln(Fobs/|Fc_total|)`, `x = stol2`. Linear regression → `B_iso = -slope`, `k = exp(intercept)`.
Initial k_sol=0.35, B_sol=46.0.

### Analytical derivatives

```
dy/dk_overall = fe = |Fc_total| * kaniso
dy/dk_sol     = solv_b * Re(Fc_total* . Fmask) / |Fc_total| * k_overall * kaniso
dy/dB_sol     = -stol2 * k_sol * solv_b * Re(Fc_total* . Fmask) / |Fc_total| * k_overall * kaniso
dy/dB*_ij:    du.B11 = -0.25*y*h², du.B22 = -0.25*y*k², du.B33 = -0.25*y*l²
              du.B12 = -0.50*y*h*k, du.B13 = -0.50*y*h*l, du.B23 = -0.50*y*k*l
dy/d(d_j) = dot(constraint_basis[j], du)
```

### LM settings

lambda_start=0.001, lambda_up=10, lambda_down=0.1, lambda_limit=1e15, eval_limit=100,
convergence=1e-5 relative WSSR change (2 consecutive). Marquardt damping: `alpha[j][j] *= (1+lambda)`.

---

## 8. `sigma_a.rs` — Simplified Binned SigmaA (v1)

### Algorithm

1. **Bin**: 20 bins equally spaced in 1/d². Bin_width = (1/d_min²)/20.
2. **Accumulate** (working set only): `sum_fo_fc[b] += |Fo|*|Fc|`, `sum_fc2[b] += |Fc|²`, `count[b]++`.
3. **Merge underpopulated bins** (<3 reflections) into nearest populated neighbor.
4. **Compute D**: `D[b] = clamp(sum_fo_fc[b] / sum_fc2[b], 0.0, 1.0)`. Empty bins → interpolate from neighbors.
5. **Smooth**: 3-point moving average. Edges weighted 2:1.
6. **Monotonicity**: `D[i] = min(D[i], D[i-1])` for i=1..19.
7. **Compute σ²**: `σ²[b] = Σ(|Fo| - D[b]*|Fc|)² / N_bin`. Floor at 1e-6. Smooth same as D (no monotonicity).
8. **Interpolate**: Linear between bin centers for per-reflection values.

Use f64 accumulators to avoid precision loss.

### Output

```rust
pub struct SigmaAResult {
    pub d_bins: [f32; 20],
    pub sigma_sq_bins: [f32; 20],
    pub s2_min: f32,
    pub bin_width: f32,
}
```

---

## 9. `map_coefficients.rs` — Weighted Map Coefficients

### Figure of merit

`x = 2 * |Fo| * D * |Fc| / (2*σ_Fo² + ε*σ²)`, then:
- Acentric: `m = I1(x)/I0(x)`
- Centric: `m = tanh(x)`
- (For v1: treat all as acentric — Sohncke groups have very few centric reflections)

### Map coefficients

- **2mFo-DFc**: `(2m*|Fo| - D*|Fc|) * exp(i*φ_calc)`
- **mFo-DFc**: `(m*|Fo| - D*|Fc|) * exp(i*φ_calc)`
- Missing Fo: `F_best = D*Fc`, `F_diff = 0`

Inverse FFT → real-space density grids.

---

## 10. `targets.rs` — ML Target + B-factor Gradient

### ML target (negative log-likelihood)

Acentric: `-log L = log(d) + (Fo² + D²Fc²)/d - log(I0(2FoDFc/d))`
where `d = 2σ_Fo² + ε*σ²`.

### R-factors

`R = Σ||Fo| - k|Fc|| / Σ|Fo|` over working (R-work) or free (R-free) set.

### B-factor gradient via FFT

1. **B-gradient map coefficients**: `F_Bgrad(h) = F_grad(h) * (-stol2)` where
   `F_grad = (m*|Fo|*exp(iφ) - D*Fc)` and `stol2 = (1/d²)/4`.
2. **Inverse FFT** → real-space B-gradient map.
3. **Trilinear interpolation** at each atom's fractional position (8 corners, periodic wrapping).
4. **Gradient**: `dL/dB_i = occ_i * interpolated_value / V_cell`.

The `-stol2` multiplication converts the position gradient into a B-factor gradient.
This costs one extra per-reflection multiply and one extra inverse FFT.

### argmin Integration

- Implement `CostFunction` + `Gradient` on same struct with `RefCell<Option<CachedEval>>` for caching.
- `gradient()` runs full pipeline, caches both cost and gradient.
- `cost()` checks cache first, returns cached value if params match.
- B-factor bounds [2, 300] via sigmoid reparameterization: `B = 2 + 298*sigmoid(x)`, chain rule `dL/dx = dL/dB * dB/dx`.
- Solver: `LBFGS::new(MoreThuenteLineSearch::new(), 7).with_tolerance_grad(1e-6)`.
- Parameter type: `Vec<f64>`.

### Macro-cycle structure

```
for macro_cycle in 0..5 {
    compute Fc (splat + FFT + deblur)
    compute Fmask (mask + FFT)
    fit scaling (LM, all params simultaneously)
    update sigmaA (binned estimation)
    run argmin on B-factors for ≤15 iterations
    check R-free for convergence
}
```
