# 1D Solver for Respiratory Airflow Simulations with Distal Time-Dependent Boundary Conditions

This solver computes air **pressure, flow and cross-sectional area** in a
deformable 1D network of airways (trachea → terminal bronchi) over the
breathing cycle.

It is a modified version of SimVascular's
[svOneDSolver](https://github.com/SimVascular/svOneDSolver), originally written
for blood flow in arteries. The equations are the same 1D model: conservation
of mass, axial momentum balance and a constitutive law that links pressure to
cross-sectional area, solved with a finite element method. The changes for airways are:

- a new wall law, `AIRWAYS`, for compliant bronchi;
- flow conditions imposed at **every terminal branch** instead of a single inlet flow;
- periodic boundary conditions over a breathing cycle of chosen length.

The full list of changes is in [Changes with respect to svOneDSolver](#changes-with-respect-to-svonedsolver).

> **Citation.** If you use this code, please cite the associated paper
> ([doi:10.46298/mbj.19176](https://doi.org/10.46298/mbj.19176)) and the original SimVascular paper:
> Updegrove A. et al., *SimVascular: An open source pipeline for cardiovascular
> simulation*, Annals of Biomedical Engineering 45:525–541 (2017),
> doi:10.1007/s10439-016-1762-8.

---

## Contents

```
1D-Solver-for-Respiratory-Airflow-Simulations-with-Distal-Time-Dependant-Boundary-Conditions/
├── CMakeLists.txt              build configuration
├── Code/Source/                C++ sources of the solver
│   ├── main.cxx                reads the input file and starts the simulation
│   ├── cvOneDBFSolver.*        time loop, Newton iterations, output files
│   ├── cvOneDMaterialAirways.* airway wall law (added for this work)
│   ├── cvOneDMth*Model.*       finite element equations and boundary conditions
│   ├── cvOneDSubdomain.*       one airway segment (mesh, material, imposed flow)
│   └── sparse/                 optional sparse linear solvers (csparse, superlu)
├── examples/bifurcation/       small test case (trachea + 2 bronchi)
├── CITATION.cff                citation metadata
└── LICENSE                     SimVascular license (MIT-style)
```

## Requirements

| Dependency | Version tested | Notes |
|---|---|---|
| C++ compiler (g++ or clang++) | g++ 11.4 | must support C++11 |
| CMake | 3.22 | 3.5 or newer |
| make | GNU make 4.3 | |
| Python 3 + NumPy | 3.10 | *optional*, only used to check the example |

There are no other dependencies with the default linear solver. The code was
developed and tested on Linux (Ubuntu 22.04). It should also build on macOS.

On Ubuntu/Debian:

```bash
sudo apt-get install build-essential cmake python3-numpy
```

## Building

Build in a separate folder:

```bash
git clone https://github.com/QuentinHerszkowicz/1D-Solver-for-Respiratory-Airflow-Simulations-with-Distal-Time-Dependant-Boundary-Conditions.git
cd 1D-Solver-for-Respiratory-Airflow-Simulations-with-Distal-Time-Dependant-Boundary-Conditions
mkdir build && cd build
cmake ..
make -j4
```

This creates the executable `build/bin/OneDSolver`.

### Linear solver option

The linear system is solved with the built-in **skyline** solver by default.
All the simulations in the paper used this solver. Two other options exist:

| `sparseSolverType` | Description | External library |
|---|---|---|
| `skyline` (default) | serial skyline solver | none |
| `csparse` | serial sparse solver, sources included | none |
| `superlu` | parallel sparse solver (shared memory) | [SuperLU_MT](https://portal.nersc.gov/project/sparse/superlu/), built with PTHREAD |

```bash
cmake .. -DsparseSolverType=csparse
# or
cmake .. -DsparseSolverType=superlu -DSUPERLU_DIR=/path/to/SuperLU_MT_3.x
```

## Running

```bash
/path/to/build/bin/OneDSolver input_file.in
```

Results are written to the **current folder**, so run the solver from the
folder where you want the output files.

### Quick test

```bash
cd examples/bifurcation
../../build/bin/OneDSolver bifurcation.in > log.txt
python3 check_results.py
```

The run takes about 25 s on one core. `check_results.py` should print `OK`. The
expected values are a peak tracheal flow of about 392.7 mL/s (two terminal
branches at 196.3 mL/s each), a flow balance error at the bifurcation below
1e-6 mL/s, and a peak distal pressure of about 7.37 dyn/cm².

The solver prints a lot of diagnostic text to the console (nodes, joints,
segments, convergence at each time step). Redirect it to a file as shown above.

## Input file format

The input is a plain text file. Each line starts with a keyword. Lines starting
with `#` are comments. All quantities are in **CGS units** (cm, g, s; pressure
in dyn/cm², flow in mL/s = cm³/s).

A complete commented example is
[examples/bifurcation/bifurcation.in](examples/bifurcation/bifurcation.in).

### Geometry

```
NODE     name  x  y  z
SEGMENT  name  id  length  nElements  inNode  outNode  inArea  outArea  initialFlow \
         material  lossType  branchAngle  upstreamSeg  branchSeg  boundaryType  dataTable
JOINT    name  node  inletListName  outletListName
JOINTINLET   inletListName   nSegments  seg1 seg2 ...
JOINTOUTLET  outletListName  nSegments  seg1 ...
```

- **Orientation.** In the airway trees used in the paper, every segment goes
  from its **distal** node (`inNode`) to its **proximal** node (`outNode`). At a
  bifurcation, the two daughter branches are listed in `JOINTINLET` and the
  parent branch in `JOINTOUTLET`.
- `inArea`, `outArea`: reference areas (cm²) at the two ends of the segment.
  If they differ, the area changes smoothly along the segment (see below).
- `lossType NONE`, `branchAngle 0.0`, `upstreamSeg 0`, `branchSeg 0`: not used
  for airways. Keep these values.
- `boundaryType`, applied at the free end of a segment:
  - `PRESSURE <table>`: imposed pressure (used at the trachea);
  - `FLOW <table>`: imposed flow (used at every terminal branch);
  - `NOBOUND NONE`: internal segment.

### Boundary condition tables

```
DATATABLE  name  LIST
t0  value0
t1  value1
...
ENDDATATABLE
```

The values are interpolated linearly in time and repeated with the period
`cycleTime` (see `SOLVEROPTIONS`). A `FLOW` table must cover the whole interval
`[0, cycleTime]`.

### Solver options

```
SOLVEROPTIONS  dt  saveFreq  maxSteps  quadPoints  cycleTime  tolerance  useIV  useStab
```

| Field | Meaning | Paper value |
|---|---|---|
| `dt` | time step (s) | 0.0008–0.002, depending on the subject |
| `saveFreq` | save results every `saveFreq` steps | 20 |
| `maxSteps` | total number of time steps | 6000 or 10000 |
| `quadPoints` | Gauss points per element | 2 |
| `cycleTime` | period of the boundary conditions (s) | 6.0 |
| `tolerance` | Newton convergence tolerance | 1e-5 |
| `useIV` | include the viscous term (1 = yes) | 1 |
| `useStab` | use stabilisation (1 = yes) | 1 |

### Material

```
MATERIAL  name  AIRWAYS  density  viscosity  pRef  exponent  k1  k2  k3
```

Values used for air: density 0.0012 g/cm³, dynamic viscosity 1.85e-4 g/(cm·s),
`pRef` 0, `exponent` 2.0 (velocity profile exponent).

The `AIRWAYS` law (`Code/Source/cvOneDMaterialAirways.cxx`) relates area `S`
and pressure `p` at position `z`:

```
S(p, z) = S0(z) * [ 1 + (p - pRef) / E(z) ]^(1/4)
E(z)    = 4/3 * ( K1 * r0(z) + K3 ),    r0 = sqrt(S0 / pi)
```

`S0(z)` is the reference area. It goes from `inArea` to `outArea` along a tanh
profile and is constant when both areas are equal.

**Parameter reading.** The current input parser assigns the **7th field
(`k1`) to both `K1` and `K3`**. `K2` is set to 0, and the `k2` and `k3` fields
of the line are ignored. The parameters actually used are printed at start-up,
for example `Setting material K's 211748 0 211748`.

### Output

```
OUTPUT  TEXT | VTK | BOTH  [vtkType]
```

- `TEXT` writes one file per segment and per field: `<seg>_flow.dat`,
  `<seg>_pressure.dat`, `<seg>_area.dat`, `<seg>_Re.dat` (Reynolds number) and
  `<seg>_wss.dat` (wall shear stress). Each row is a finite element node, from
  `inNode` (first row) to `outNode` (last row). Each column is a saved time step.
- `VTK` writes ParaView files. `vtkType` = 0: 3D, one file per time step;
  1: 3D, a single file; 2: 1D polylines, one `.vtp` file per time step.
- `BOTH` writes both.

## Changes with respect to svOneDSolver

This code is based on svOneDSolver, commit
[`3598cb55`](https://github.com/SimVascular/svOneDSolver/commit/3598cb55) (July 2022).

- **New material `AIRWAYS`** (`cvOneDMaterialAirways.h/.cxx`, registered in
  `cvOneDMaterialManager`, `cvOneDModelManager`, `cvOneDEnums.h` and `main.cxx`).
- **Flow imposed at each terminal segment.** `cvOneDSegment` and
  `cvOneDSubdomain` store their own flow table (`SetFlowRate` /
  `GetFlowRate`). `cvOneDMthModelBase` applies these flows and sums them for
  the mass balance. The single global inlet condition of the original code
  was removed.
- **`SOLVEROPTIONS` format.** The `inletDataTable` and `boundaryType` fields
  were replaced by `cycleTime`.
- **1D VTK output** (`OUTPUT ... 2`) enabled.
- Extra console diagnostics (nodes, joints, segments, cycle number).
  `debugMode` is on by default.

## Reproducibility

The simulations in the paper were run with the default `skyline` solver,
compiled with g++ 11.4 on Ubuntu 22.04 with the flags in `CMakeLists.txt`
(`-O3 -std=c++0x`). The `examples/bifurcation` case gives bit-for-bit
identical results with this repository and with the executable used for the
paper.

## License

The code is distributed under the SimVascular license (MIT-style permissive
license); see [LICENSE](LICENSE). The original copyright belongs to Stanford
University, The Regents of the University of California, and others.
