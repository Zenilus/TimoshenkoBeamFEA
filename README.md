## MATLAB Timoshenko Beam FEA Toolkit

This toolkit implements a two-dimensional Timoshenko beam finite element workflow in MATLAB. A single driver script reads structural data from an Excel workbook, assembles and solves the global system, computes section-level response metrics, and renders publication-ready visualizations. The code favors clarity and modularity so that you can adapt it to coursework, research prototypes, or rapid what-if studies without reverse engineering a monolithic script.

## Contents
- [Why Timoshenko?](#why-timoshenko)
- [Feature Highlights](#feature-highlights)
- [Repository Structure](#repository-structure)
- [Software Requirements](#software-requirements)
- [Quick Start](#quick-start)
- [Input Workbook Specification](#input-workbook-specification)
- [Running Analyses](#running-analyses)
- [Understanding the Outputs](#understanding-the-outputs)
- [Customising the Workflow](#customising-the-workflow)
- [Verification Guidance](#verification-guidance)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

---

## Why Timoshenko?
Classical Euler-Bernoulli theory neglects shear deformation, which becomes inaccurate for thick beams, short spans, laminated sections, or materials with low shear stiffness. Timoshenko beam theory augments the formulation with shear flexibility so bending and shear deformations are solved concurrently. The additional fidelity is critical when:
- span-to-depth ratios drop below roughly 10
- laminated timber, composite, or sandwich panels are analysed
- high shear loads drive serviceability checks
- accurate rotations are needed for connection design

This repository packages that theory into a MATLAB code base with a clear control flow, keeping advanced behaviour accessible to students and practitioners.

---

## Feature Highlights
- End-to-end pipeline from Excel input to figures and console summaries.
- Modular source tree with specific directories for I/O, analysis, reporting, and plotting.
- Automatic application of self-weight via consistent nodal loads and end moments.
- Stress post-processing that tracks controlling element, plane, and end metadata for every extreme value.
- Visualization of both undeformed geometry and scaled deformed shapes with displacement magnitude gradients.
- Well-documented functions that support incremental extension (new elements, shapes, or reporting formats).

---

## Repository Structure
| Path | Description |
| ---- | ----------- |
| `scripts/TimoshenkoBeamFEA.m` | Entry point that sets up the MATLAB path, reads input data, runs the solver, and orchestrates reporting/plots. |
| `src/analysis/` | Core finite element routines for stiffness assembly, solving, and stress recovery. |
| `src/io/` | Excel parsing utilities that normalise sheet layouts and property definitions. |
| `src/properties/` | Section property calculators tailored to Timoshenko theory (e.g., shear correction factors). |
| `src/reporting/` | Displacement summaries, extrema extraction helpers, and text formatting for console output. |
| `src/visualization/` | Plotters for the undeformed system, deformed shapes, and colour mapping helpers. |
| `data/` | Drop input workbooks such as `input_data.xlsx`; not tracked by default. |
| `docs/` | Detailed module documentation and implementation notes. |

---

## Software Requirements
- MATLAB R2020a or newer (base MATLAB only). Earlier versions may work but are untested.
- Microsoft Excel (or any editor capable of producing `.xlsx` files) for preparing inputs.
- Optional: familiarity with finite element terminology to interpret outputs and customise loads.

---

## Quick Start
1. Clone or download this repository and open MATLAB in the project root.
2. Copy the sample workbook template found in `data/README.md` or build your own using the [Input Workbook Specification](#input-workbook-specification).
3. From the MATLAB command window run:
   ```matlab
   TimoshenkoBeamFEA
   ```
   Provide an absolute or relative workbook path to analyse other files:
   ```matlab
   TimoshenkoBeamFEA('data/my_load_case.xlsx')
   ```
4. Review the console summary for displacements and stress extrema, then inspect the generated figures for geometry and deformed shape validation.

---

## Input Workbook Specification
The solver expects a single Excel workbook with named sheets. Column order matters, headers must be present on the first row, and values below row 2 define the model. Adopt SI units (metres, Newtons, Pascals, kg/m³) for consistency.

### Nodes
| NodeID | X | Y |
| ------ | - | - |
| 1 | 0.0 | 0.0 |
| 2 | 2.0 | 0.0 |

- `NodeID` must be unique, positive integers.
- Coordinates are 2D Cartesian positions in metres.

### Elements
| ElementID | Node1 | Node2 |
| --------- | ----- | ----- |
| 1 | 1 | 2 |
| 2 | 2 | 3 |

- Each row connects two existing nodes. Orientation is derived from nodal coordinates.
- Elements are straight Timoshenko beams with 3 degrees of freedom per node (`ux`, `uy`, `theta`).

### Supports
| NodeID | Type |
| ------ | ---- |
| 1 | Fixed |
| 3 | Roller |

- Supported types: `Fixed`, `Pinned`, `Roller` (case-insensitive).
- Combine types with repeated rows if a node carries different constraints (e.g., `Pinned` plus `Roller` is equivalent to `Fixed`).

### Forces
| NodeID | Fx | Fy | Mz |
| ------ | -- | -- | -- |
| 2 | 0 | -5000 | 0 |
| 3 | 1000 | 0 | 150 |

- Forces are nodal loads in Newtons; `Mz` is an out-of-plane nodal moment in Newton-metres.
- Positive axes use the right-hand convention: +X to the right, +Y upward, +Z out of the page (counter-clockwise positive moment).
- The column `Mz` must exist even when zero to maintain consistent load vector assembly.

### Properties
At minimum supply:

| YoungsModulus | CrossSectionalArea | Density |
| ------------- | ------------------ | ------- |
| 2.1e11 | 8.0e-3 | 7850 |


| SectionType | Width | Height | Diameter | ShearModulus | PoissonRatio |
| ----------- | ----- | ------ | -------- | ------------ | ------------ |
| Rectangle | 0.1 | 0.2 | | 8.1e10 | 0.3 |

- Supported section types: `Rectangle`, `Square`, `Circle`. Unknown types fall back to conservative rectangular assumptions.



> **Tips**
> - Keep IDs consecutive to simplify debugging.
> - Start each sheet at row 1 with headers; MATLAB treats empty header cells as unnamed variables and the importer will reject them.
> - Set `Density` to zero if you want to omit self-weight from the load case.

---

## Running Analyses
- Execute `TimoshenkoBeamFEA` from MATLAB. The script automatically adds `src/` to the path, reads the workbook, and orchestrates the workflow.
- To reuse the solver in other scripts, call `performTimoshenkoFEA` directly:
  ```matlab
  workbookPath = fullfile('data', 'input_data.xlsx');
  [nodes, elements, supports, forces, E, A, rho, sectionType, w, h, d, G, nu] = readInputData(workbookPath);
  [displacements, stresses] = performTimoshenkoFEA(nodes, elements, supports, forces, E, A, rho, sectionType, w, h, d, G, nu);
  ```
- The solver returns displacement and stress structures that you can feed into custom plotting, optimisation loops, or export functions.

---

## Understanding the Outputs
**Console summary**
- Global counts, material constants, and derived section properties.
- Maximum horizontal and vertical displacements plus extreme rotations with qualitative flags (small/moderate/large).
- Peak bending moments, bending stresses, shear forces, shear stresses, and von Mises stresses with controlling element, plane, and end labels.

**Figures**
- `plotBeamSystem`: undeformed geometry, support types, and applied forces/moments for model validation.
- `plotDeformedShape`: overlaid deformed configuration with automatic scaling, displacement magnitude colour bar, and annotations for maximum displacement.

**Programmatic data**
- `displacements`: global DOF vector ordered `[u1, v1, theta1, u2, ...]`.
- `stresses.elemental`: von-Mises-style element envelope combining axial, bending, and shear contributions.
- `stresses.bending_details`, `stresses.shear_details`, `stresses.von_mises_details`: metadata structs with indices and signed values for each plane/end combination, suitable for automated reporting.

Store the returned struct for post-processing pipelines such as fatigue checks, design code verifications, or report generation.

---

## Customising the Workflow
- **Alternative load cases**: Duplicate the workbook, adjust node forces or densities, and rerun the driver without altering MATLAB files.
- **Element refinement**: Add intermediate nodes to capture curvature and local effects. The solver automatically adapts to the new topology.
- **Section library extensions**: Enhance `calculateTimoshenkoSectionProperties.m` with profiles such as I-beams or hollow rectangles and supply appropriate shear correction factors.
- **Spatially varying properties**: Modify `readInputData` and `performTimoshenkoFEA` to accept per-element vectors for `E`, `A`, or density when modelling tapered or composite members.
- **Reporting/export**: Extend the helpers in `src/reporting/` to write CSV/Excel summaries or integrate with dashboards.
- **Visual styling**: Adjust `plotBeamSystem` or `plotDeformedShape` to match publication branding or to add annotations.

---

## Verification Guidance
- Validate the code with canonical problems: cantilever under tip load, simply supported beam under uniform load, and deep beam scenarios where shear deformation is significant.
- Compare deflections, rotations, shear forces, and bending moments against textbook Timoshenko solutions or high-fidelity FEA packages.
- Document load cases, solver outputs, and plots in `docs/` to maintain a traceable verification record.
- If you extend the section library or solver, add regression tests (e.g., MATLAB scripts comparing analytical values) to ensure future changes remain trustworthy.

---

## Troubleshooting
- **Workbook import error**: Confirm the file path, sheet names, and headers match the specification. Empty header cells cause the importer to fail.
- **Ill-conditioned stiffness warning**: Indicates missing supports or duplicated DOFs. Review boundary conditions and ensure each rigid body motion is restrained.
- **Unexpectedly small displacements**: Check unit consistency (metres vs millimetres) and verify load magnitudes.
- **Missing forces in plots**: Ensure every force references an existing node ID and that `Mz` is present even when zero.
- **Section modulus warning**: Occurs when supplied geometric dimensions do not match the implied area within ±20%. Revisit `Width`, `Height`, or `Diameter` entries.

---

## Contributing
Contributions are welcome. Please open an issue describing the proposed enhancement or bug fix, link relevant references or validation cases, and submit a pull request with clear commit messages. Include verification evidence (plots, numerical comparisons) for solver-level changes.

---

## License
Distributed for educational and research use. Adapt and extend it for academic projects or professional prototypes at your discretion.
