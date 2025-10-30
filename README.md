# MATLAB Timoshenko Beam FEA Toolkit

This repository contains a modular MATLAB implementation of a two-dimensional Timoshenko beam finite element analysis (FEA) workflow. The toolkit reads a beam definition from an Excel workbook, assembles the global stiffness system, solves for nodal displacements, evaluates internal forces and stresses, and visualizes both the original and deformed beam configurations. The emphasis is on clarity and extensibility so the code can support coursework, research prototypes, or quick structural what-if studies.

## Highlights
- End-to-end MATLAB pipeline from Excel input to publication-ready plots
- Modular source tree that invites experimentation with material laws or elements
- Built-in utilities for post-processing, visualization, and result inspection
- Verification guidance to compare with analytical Timoshenko beam solutions

## Table of Contents
- [Why Timoshenko Theory?](#why-timoshenko-theory)
- [Repository Layout](#repository-layout)
- [Prerequisites](#prerequisites)
- [Preparing the Input Workbook](#preparing-the-input-workbook)
- [Quick Start](#quick-start)
- [Understanding the Outputs](#understanding-the-outputs)
- [Extending or Customizing](#extending-or-customizing)
- [Testing and Verification](#testing-and-verification)
- [Troubleshooting](#troubleshooting)
- [References](#references)
- [Contributing](#contributing)
- [License](#license)

---

## Why Timoshenko Theory?
Classical Euler-Bernoulli beam theory ignores shear deformation, which can lead to large errors for thick beams, composite members, or low shear-stiffness materials. Timoshenko theory remedies this by coupling bending and shear effects. The provided MATLAB scripts implement this richer model while keeping the analysis pipeline approachable.

---

## Repository Layout
- `scripts/` – Entry points runnable from MATLAB. The default driver is `scripts/TimoshenkoBeamFEA.m` which adds the `src` tree to the path and kicks off a full analysis.
- `src/` – Reusable library code split by concern:
  - `analysis/` – Core finite element routines (`performTimoshenkoFEA.m`, `calculateTimoshenkoStiffness.m`, `calculateTimoshenkoStresses.m`).
  - `io/` – Data loading helpers (`readInputData.m`).
  - `properties/` – Cross-section property calculations (`calculateTimoshenkoSectionProperties.m`).
  - `reporting/` – Console summaries (`displayResults.m`).
  - `visualization/` – Plotting utilities (`plotBeamSystem.m`, `plotDeformedShape.m`, `getColorFromValue.m`).
- `data/` – Place analysis workbooks such as `input_data.xlsx`.
- `docs/` – Project documentation, figures, and references.

---

## Prerequisites
- MATLAB R2020a or newer (earlier versions may work but have not been verified).
- No special toolboxes are required beyond base MATLAB.
- Microsoft Excel or any editor capable of creating `.xlsx` files in the expected layout.
- Optional: familiarity with basic structural analysis terminology.

---

## Preparing the Input Workbook
Create an `input_data.xlsx` file under `data/`. Each sheet feeds a different part of the solver; column order matters. Units are assumed to be SI (meters, Newtons, Pascals, kg/m³).

### 1. `Nodes`
| NodeID | X | Y |
|--------|---|---|
| 1 | 0.0 | 0.0 |
| 2 | 2.0 | 0.0 |
| ⋮ | ⋮ | ⋮ |

- `NodeID` must be unique integers.
- Coordinates are given in meters.

### 2. `Elements`
| ElementID | Node1 | Node2 |
|-----------|-------|-------|
| 1 | 1 | 2 |
| 2 | 2 | 3 |

- Elements are straight-line members connecting two existing nodes.
- Orientation is inferred from node coordinates.

### 3. `Supports`
| NodeID | Type |
|--------|------|
| 1 | Fixed |
| 3 | Roller |

- `Type` accepts `Fixed`, `Pinned`, or `Roller` (case-insensitive).
- Multiple supports can be assigned across the model.

### 4. `Forces`
| NodeID | Fx | Fy | Mz |
|--------|----|----|----|
| 2 | 0 | -5000 | 0 |
| 3 | 1000 | 0 | 150 |

- Forces are applied at nodes in Newtons and follow the right-handed sign convention: positive `Fx` acts along +X (right), positive `Fy` acts along +Y (up), and positive `Mz` induces a counter-clockwise rotation about +Z (out of plane).

### 5. `Properties`
At minimum include:

| YoungsModulus | CrossSectionalArea | Density |
|---------------|--------------------|---------|
| 2.1e11 | 8.0e-3 | 7850 |

Optional columns (recommended for high accuracy):

| SectionType | Width | Height | Diameter | ShearModulus | PoissonRatio |
|-------------|-------|--------|----------|--------------|--------------|
| Rectangle | 0.1 | 0.2 | | 8.1e10 | 0.3 |

- `SectionType` supports `Rectangle`, `Square`, or `Circle`.
- Provide either `ShearModulus` directly or `PoissonRatio` so `G` can be derived.
- Leaving geometric columns blank triggers conservative defaults based on area.

> **Tip:** Keep IDs consecutive and start all tables on row 2, with headers on row 1.

---

## Quick Start
1. Open MATLAB and change the working directory to the repository folder.
2. Review `data/input_data.xlsx` (or duplicate your preferred template workbook) and confirm IDs, units, and load directions.
3. Launch the driver script:
  ```matlab
  TimoshenkoBeamFEA
  ```
  Use `TimoshenkoBeamFEA('relative/or/absolute/path.xlsx')` to analyze a different workbook without relocating it.
4. Inspect the MATLAB command window and generated figures to validate the model setup before moving to detailed studies.

---

## Understanding the Outputs
- **Console log**
  - Beam summary (counts, material constants, derived section properties).
  - Max horizontal/vertical displacements and rotations with contextual notes.
  - Peak combined stress, bending moment, shear force, and shear stress per element.
- **Plots**
  - Undeformed layout for quick verification of geometry, support types, and load directions.
  - Deformed shape using an automatic scale factor; legend indicates the magnification.
  - Colorbar reports actual displacement magnitudes (meters).

If you need numerical displacements or stresses for post-processing, capture the outputs directly by calling the lower-level API:
```matlab
workbookPath = fullfile('data', 'input_data.xlsx');
[nodes, elements, supports, forces, E, A, rho, sectionType, w, h, d, G, nu] = readInputData(workbookPath);
[displacements, stresses] = performTimoshenkoFEA(nodes, elements, supports, forces, E, A, rho, sectionType, w, h, d, G, nu);
```

---

## Extending or Customizing
- **Different load cases**: Duplicate and modify `data/input_data.xlsx`, then rerun the driver.
- **Alternate visualization**: Use `displacements` and `stresses` from `performTimoshenkoFEA` to script custom plots or animations.
- **Element refinement**: Add more nodes/elements to capture curvature or localized loads with greater fidelity.
- **Cross-section library**: Extend `calculateTimoshenkoSectionProperties.m` with additional shapes (e.g., I-beams, hollow rectangles) and return the appropriate shear correction factors.
- **Non-uniform properties**: The current implementation assumes constant material/section properties. To vary them per element, pass vectors of `E`, `A`, etc., and adapt the loops accordingly.

---

## Testing and Verification
- Start with a canonical cantilever or simply supported beam problem that has a published Timoshenko solution.
- Compare tip deflections, rotations, and bending moments against textbook results before exploring complex geometries.
- For peer review or coursework submissions, capture load cases, numerical results, and plots in `docs/` for traceability.

---

## Troubleshooting
- **"Error reading input file"** – Ensure `data/input_data.xlsx` exists and sheet headers match the documented format.
- **Ill-conditioned stiffness warning** – Typically indicates missing supports or an unconstrained degree of freedom. Double-check boundary conditions.
- **Zero or tiny displacements** – Units may be inconsistent (e.g., millimeters instead of meters) or the load magnitudes are too small relative to stiffness.
- **Plots missing forces/supports** – Verify the IDs in `Supports` and `Forces` correspond to valid node IDs.

---

## References
- S. P. Timoshenko, *Strength of Materials*, Part 1.
- J. Gere & S. Timoshenko, *Mechanics of Materials*.
- A. Ghali, A. Neville, T. G. Brown, *Structural Analysis: A Unified Classical and Matrix Approach*.

---

## Contributing
Feel free to fork the project, extend section libraries, or add verification examples. When sharing improvements, include a brief description of the load case or theoretical reference used for validation.

---

## License
This project is distributed for educational and research purposes. Adapt as needed for your academic or professional projects.
