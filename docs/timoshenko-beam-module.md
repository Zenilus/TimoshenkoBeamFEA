# Timoshenko Beam FEA Module Documentation

This document complements the project README by detailing the workflow, data contracts, and responsibilities of each MATLAB component in the Timoshenko beam finite element toolkit. Use it when you are extending the codebase, integrating it into larger studies, or validating solver behaviour.

## Contents
- [1. Overview](#1-overview)
- [2. System Architecture](#2-system-architecture)
- [3. Data Contracts](#3-data-contracts)
- [4. Module Reference](#4-module-reference)
- [5. Algorithmic Notes](#5-algorithmic-notes)
- [6. Diagnostics and Error Handling](#6-diagnostics-and-error-handling)
- [7. Extending the Toolkit](#7-extending-the-toolkit)
- [8. Verification Strategy](#8-verification-strategy)
- [9. Glossary](#9-glossary)

---

## 1. Overview
The toolkit solves planar beam problems using Timoshenko theory, accounting for both bending and shear deformations. A single MATLAB driver handles path setup, workbook parsing, solver execution, reporting, and plotting. Each responsibility is isolated in its own function so that contributors can experiment without destabilising the full workflow.

Key goals of the module architecture:
- keep input preparation lightweight (Excel-based)
- ensure solver internals are transparent and traceable
- surface rich metadata about response maxima for automated reporting
- provide clean extension points for new cross-sections, load cases, or renderers

---

## 2. System Architecture

### 2.1 Execution Pipeline
```
TimoshenkoBeamFEA (script entry)
  ├─ readInputData (Excel ingestion)
  ├─ plotBeamSystem (pre-analysis plot)
  ├─ performTimoshenkoFEA (solver core)
  │    ├─ calculateTimoshenkoSectionProperties
  │    ├─ calculateTimoshenkoStiffness (per element)
  │    └─ calculateTimoshenkoStresses (post-processing)
  ├─ displayResults (displacement summary)
  ├─ plotDeformedShape → getColorFromValue
  └─ prepareStressReportEntries
         ├─ computeBendingExtremaSummary
         ├─ computeShearExtremaSummary
         └─ computeVonMisesExtremaSummary
```

### 2.2 Directory Responsibilities
- `scripts/`: UC-level entry points and orchestration logic.
- `src/io/`: External data ingestion and normalisation.
- `src/properties/`: Cross-section derived properties and shear correction factors.
- `src/analysis/`: Element matrices, global assembly, solution, and stress recovery.
- `src/reporting/`: Summaries, metadata extraction, and formatting helpers.
- `src/visualization/`: Diagnostic and presentation-quality figures.

---

## 3. Data Contracts

### 3.1 Input Workbook
- **Nodes**: columns `NodeID`, `X`, `Y`. IDs must be consecutive integers; coordinates in metres.
- **Elements**: columns `ElementID`, `Node1`, `Node2`. Node references must exist in the Nodes sheet.
- **Supports**: columns `NodeID`, `Type` (`Fixed`, `Pinned`, `Roller`). Types are case-insensitive.
- **Forces**: columns `NodeID`, `Fx`, `Fy`, `Mz`. Forces in Newtons, moment in Newton-metres. The `Mz` column must be present even if zero.
- **Properties**: minimum `YoungsModulus`, `CrossSectionalArea`, `Density`. Optional fields `SectionType`, `Width`, `Height`, `Diameter`, `ShearModulus`, `PoissonRatio` refine section behaviour.

Each sheet begins with headers on row 1. Data spans rows 2..N. Empty trailing rows are ignored by MATLAB's table reader.

### 3.2 Solver Outputs
- `displacements`: global vector `[u1, v1, theta1, u2, ...]` (metres, radians).
- `stresses`: struct containing element-wise response metrics:
  - `elemental`: combined von-Mises-style stress (Pa).
  - `bending_moments`, `shear_forces` (N·m, N) and their derived stresses (Pa).
  - `bending_details`, `shear_details`, `von_mises_details`: metadata describing which element, plane, and end produced the maxima.
  - `element_positions`: mid-span coordinates for plotting or sorting.

---

## 4. Module Reference

### 4.1 `scripts/TimoshenkoBeamFEA.m`
- **Responsibility**: User-facing entry point. Validates inputs, prepares MATLAB paths, invokes solver, triggers reporting and figures.
- **Inputs**: optional workbook path (string). Defaults to `data/input_data.xlsx` relative to the repository root.
- **Outputs**: none (side-effects: console output, figure windows). Returns displacement and stress structs only when invoked through nested calls.
- **Key Steps**:
  1. Add `src/` to MATLAB path if absent.
  2. Validate workbook path, falling back to sample data when omitted.
  3. Load model data via `readInputData`.
  4. Render undeformed geometry (`plotBeamSystem`).
  5. Run solver and collect `displacements`/`stresses`.
  6. Print displacement summary (`displayResults`).
  7. Plot deformed configuration (`plotDeformedShape`).
  8. Format stress extrema (`prepareStressReportEntries`) and emit console lines.
- **Extension Hooks**: Insert custom analysis steps or exports after the solver call but before figures to capture intermediate data.

**MATLAB-style pseudocode**
```matlab
function TimoshenkoBeamFEA(inputWorkbook)
  projectRoot = fileparts(mfilename('fullpath'));
  srcFolder = fullfile(projectRoot, 'src');
  addpath(genpath(srcFolder));

  if nargin < 1 || isempty(inputWorkbook)
    inputWorkbook = fullfile(projectRoot, 'data', 'input_data.xlsx');
  end
  assert(isfile(inputWorkbook), 'Input workbook not found: %s', inputWorkbook);

  [nodes, elements, supports, forces, E, A, density, sectionType, width, height, diameter, G, nu] = ...
    readInputData(inputWorkbook);

  plotBeamSystem(nodes, elements, supports, forces);

  [displacements, stresses] = performTimoshenkoFEA(nodes, elements, supports, forces, ...
    E, A, density, sectionType, width, height, diameter, G, nu);

  displayResults(displacements);
  plotDeformedShape(nodes, elements, displacements);

  entries = prepareStressReportEntries(stresses);
  for iEntry = 1:numel(entries)
    fprintf('%s\n', entries(iEntry).message);
  end
end
```

### 4.2 `src/io/readInputData.m`
- **Responsibility**: Normalise Excel workbook content into MATLAB arrays and scalars.
- **Validations**: checks sheet existence, required columns, and minimum dataset sizes. Throws descriptive errors when structure deviates from the contract.
- **Transformations**:
  - Converts tables to numeric matrices and string arrays.
  - Normalises support type strings and trims whitespace.
  - Derives shear modulus from Poisson's ratio when `ShearModulus` is absent.
  - Emits console summary (counts and material constants) for transparency.
- **Edge Cases**: accepts variant column casing (e.g., `nodeid`) by normalising headers. Warns if area or density is zero, since those values impact self-weight.

**MATLAB-style pseudocode**
```matlab
function [nodes, elements, supports, forces, E, A, density, sectionType, width, height, diameter, G, nu] = ...
  readInputData(inputWorkbook)

  if nargin < 1 || isempty(inputWorkbook)
    inputWorkbook = 'input_data.xlsx';
  end
  assert(isfile(inputWorkbook), 'Workbook missing: %s', inputWorkbook);

  nodesTbl     = readtable(inputWorkbook, 'Sheet', 'Nodes');
  elementsTbl  = readtable(inputWorkbook, 'Sheet', 'Elements');
  supportsTbl  = readtable(inputWorkbook, 'Sheet', 'Supports');
  forcesTbl    = readtable(inputWorkbook, 'Sheet', 'Forces');
  propsTbl     = readtable(inputWorkbook, 'Sheet', 'Properties');

  nodesTbl     = normaliseHeaders(nodesTbl);
  elementsTbl  = normaliseHeaders(elementsTbl);
  supportsTbl  = normaliseHeaders(supportsTbl);
  forcesTbl    = normaliseHeaders(forcesTbl);
  propsTbl     = normaliseHeaders(propsTbl);

  validateNodes(nodesTbl);
  validateElements(elementsTbl, nodesTbl.NodeID);
  validateSupports(supportsTbl, nodesTbl.NodeID);
  validateForces(forcesTbl, nodesTbl.NodeID);
  validateProperties(propsTbl);

  nodes    = table2array(nodesTbl(:, {'NodeID', 'X', 'Y'}));
  elements = table2array(elementsTbl(:, {'ElementID', 'Node1', 'Node2'}));

  supports.NodeID = supportsTbl.NodeID;
  supports.Type   = string(lower(strtrim(supportsTbl.Type)));

  forces = table2array(forcesTbl(:, {'NodeID', 'Fx', 'Fy', 'Mz'}));

  E         = propsTbl.YoungsModulus(1);
  A         = propsTbl.CrossSectionalArea(1);
  density   = propsTbl.Density(1);
  sectionType = string(propsTbl.SectionType(1));
  width     = propsTbl.Width(1);
  height    = propsTbl.Height(1);
  diameter  = propsTbl.Diameter(1);
  G         = propsTbl.ShearModulus(1);
  nu        = propsTbl.PoissonRatio(1);

  if isnan(G) && ~isnan(nu)
    G = E / (2 * (1 + nu));
  end

  logWorkbookSummary(nodes, elements, supports, forces, E, A, density);
end
```

### 4.3 `src/properties/calculateTimoshenkoSectionProperties.m`
- **Responsibility**: Compute second moment of area (`I`), effective shear area (`As`), maximum fibre distance, and shear correction factor (`ky`).
- **Shape Support**: `Rectangle`, `Square`, `Circle`. Fallback path assumes a rectangular profile with area-derived dimensions and warns users about conservative assumptions.
- **Usage**: Called once per analysis. Downstream routines rely on consistent `As` and `maxDistance` for stress recovery.

**MATLAB-style pseudocode**
```matlab
function [I, As, maxDistance, ky] = calculateTimoshenkoSectionProperties(sectionType, A, width, height, diameter)
  sectionType = lower(strip(sectionType));

  switch sectionType
    case 'rectangle'
      I = width * height^3 / 12;
      maxDistance = height / 2;
      ky = 5 / 6;
    case 'square'
      I = width^4 / 12;
      maxDistance = width / 2;
      ky = 5 / 6;
    case 'circle'
      I = (pi * diameter^4) / 64;
      maxDistance = diameter / 2;
      ky = 9 / 10;
    otherwise
      inferredDepth = sqrt(A);
      I = inferredDepth^4 / 12;
      maxDistance = inferredDepth / 2;
      ky = 5 / 6;
      warning('Unknown section type; assuming square with area %g m^2.', A);
  end

  As = ky * A;
end
```

### 4.4 `src/analysis/calculateTimoshenkoStiffness.m`
- **Responsibility**: Return the 6×6 local stiffness matrix for a 2-node Timoshenko beam element.
- **Formulae**:
  - Applies shear flexibility correction `phi = 12 E I / (G As L²)`.
  - Recovers Euler-Bernoulli stiffness when shear effects vanish (`phi → 0`).
- **Numerical Notes**: Terms are symmetric; rounding errors are mitigated by constructing each symmetric pair once and mirroring.

**MATLAB-style pseudocode**
```matlab
function kLocal = calculateTimoshenkoStiffness(E, G, A, As, I, L)
  phi = 12 * E * I / (G * As * L^2);
  kappa = 1 + phi;

  EA = E * A / L;
  EI = E * I / (L^3 * kappa);
  GAs = G * As / L;

  kLocal = zeros(6, 6);

  % Axial terms
  kLocal([1,4],[1,4]) = EA * [ 1, -1; -1, 1];

  % Shear-bending coupling
  kLocal([2,5],[2,5]) = GAs * [ kappa, -kappa; -kappa, kappa];
  kLocal([2,5],[3,6]) = EI * L * [ 6, -6; 6, -6];
  kLocal([3,6],[2,5]) = EI * L * [ 6, 6; -6, -6];

  % Rotational stiffness
  rotBlock = EI * L^2 * [ (4 + phi), (2 - phi); (2 - phi), (4 + phi) ];
  kLocal([3,6],[3,6]) = rotBlock;

  % Enforce symmetry explicitly
  kLocal = (kLocal + kLocal.')./2;
end
```

### 4.5 `src/analysis/performTimoshenkoFEA.m`
- **Responsibility**: End-to-end finite element solve.
- **Subtasks**:
  1. Preallocate global stiffness matrix (`K`) and load vector (`F`).
  2. Derive section properties and print summary.
  3. Loop through elements to assemble `K` using direction cosines and transformed local stiffness matrices.
  4. Map nodal loads and equivalent self-weight contributions into `F` (gravity = 9.80665 m/s² by default).
  5. Apply boundary conditions by extracting free degrees of freedom and solving `K_ff · u_f = F_f`.
  6. Call `calculateTimoshenkoStresses` for post-processing.
- **Diagnostics**: warns when `rcond(K_ff)` is below a tolerance, hinting at singular systems or poor constraint definitions.
- **Outputs**: populates the displacement vector for all DOFs (constrained DOFs remain zero) and assembles a rich `stresses` struct.

**MATLAB-style pseudocode**
```matlab
function [displacements, stresses] = performTimoshenkoFEA(nodes, elements, supports, forces, ...
  E, A, density, sectionType, width, height, diameter, G, nu)

  numNodes = size(nodes, 1);
  numElements = size(elements, 1);
  totalDof = 3 * numNodes;

  displacements = zeros(totalDof, 1);
  globalK = zeros(totalDof);
  globalF = zeros(totalDof, 1);

  [I, As, maxDistance, ky] = calculateTimoshenkoSectionProperties(sectionType, A, width, height, diameter);

  for e = 1:numElements
    nodeIds = elements(e, 2:3);
    n1 = nodeIds(1);
    n2 = nodeIds(2);

    x1 = nodes(nodes(:,1) == n1, 2:3);
    x2 = nodes(nodes(:,1) == n2, 2:3);

    delta = x2 - x1;
    L = norm(delta);
    c = delta(1) / L;
    s = delta(2) / L;

    kLocal = calculateTimoshenkoStiffness(E, G, A, As, I, L);
    T = buildTransformationMatrix(c, s);
    kGlobal = T.' * kLocal * T;

    dofIdx = elementDofIndices(nodeIds);
    globalK(dofIdx, dofIdx) = globalK(dofIdx, dofIdx) + kGlobal;

    elementWeight = density * A * L * 9.80665;
    fBodyLocal = consistentBeamLoad(elementWeight, L);
    fBodyGlobal = T.' * fBodyLocal;
    globalF(dofIdx) = globalF(dofIdx) + fBodyGlobal;
  end

  for f = 1:size(forces, 1)
    nodeId = forces(f, 1);
    dofIdx = nodeDofIndices(nodeId);
    globalF(dofIdx(1:3)) = globalF(dofIdx(1:3)) + forces(f, 2:4).';
  end

  constrained = supportsToConstrainedDof(supports);
  freeDof = setdiff(1:totalDof, constrained);

  Kff = globalK(freeDof, freeDof);
  Ff = globalF(freeDof);
  displacements(freeDof) = Kff \ Ff;

  [stresses.elemental, stresses.bending_moments, stresses.shear_forces, ...
    stresses.shear_stresses, stresses.element_positions, stresses.bending_details, ...
    stresses.shear_details, stresses.von_mises_details, stresses.von_mises] = ...
    calculateTimoshenkoStresses(nodes, elements, displacements, E, G, A, As, I, ...
      maxDistance, sectionType, width, height, diameter, density);
end
```

### 4.6 `src/analysis/calculateTimoshenkoStresses.m`
- **Responsibility**: Translate displacements into internal forces and stress measures.
- **Process**:
  - Builds transformation matrices to rotate global displacements into local element frames.
  - Computes axial strain/force and shear strain/force using Timoshenko relations.
  - Recovers bending moments at each end from the consistent stiffness formulation.
  - Subtracts self-weight contributions so reported internal forces reflect applied loading only.
  - Calculates bending stresses at top/bottom fibres and a von-Mises-style combined stress.
  - Captures metadata (element index, sign, plane, end) for maxima in bending, shear, and combined responses.
- **Sanity Checks**: compares computed section modulus against geometric inputs for supported shapes, issuing warnings when deviations exceed ±20%.

**MATLAB-style pseudocode**
```matlab
function [elemental, bendingMoments, shearForces, shearStresses, elementPositions, ...
  bendingDetails, shearDetails, vonMisesDetails, vonMises] = ...
  calculateTimoshenkoStresses(nodes, elements, displacements, E, G, A, As, I, maxDistance, ...
  sectionType, width, height, diameter, density)

  numElements = size(elements, 1);
  elemental       = zeros(numElements, 1);
  bendingMoments  = zeros(numElements, 1);
  shearForces     = zeros(numElements, 1);
  shearStresses   = zeros(numElements, 1);
  vonMises        = zeros(numElements, 1);
  elementPositions = zeros(numElements, 2);

  sectionModulus = I / max(maxDistance, eps);

  for e = 1:numElements
    nodeIds = elements(e, 2:3);
    dofIdx = elementDofIndices(nodeIds);
    dGlobal = displacements(dofIdx);

    coords = nodes(nodeIndices(nodes, nodeIds), 2:3);
    L = norm(diff(coords));
    c = (coords(2,1) - coords(1,1)) / L;
    s = (coords(2,2) - coords(1,2)) / L;

    T = buildTransformationMatrix(c, s);
    dLocal = T * dGlobal;

    axialStrain = (dLocal(4) - dLocal(1)) / L;
    avgTheta = 0.5 * (dLocal(3) + dLocal(6));
    shearStrain = (dLocal(5) - dLocal(2)) / L - avgTheta;

    axialForce = E * A * axialStrain;
    shearForce = G * As * shearStrain;

    kLocal = calculateTimoshenkoStiffness(E, G, A, As, I, L);
    fLocal = kLocal * dLocal;

    weight = density * A * L * 9.80665;
    fWeightLocal = consistentBeamLoad(weight, L);
    fLocal = fLocal - fWeightLocal;

    momentA = fLocal(3);
    momentB = fLocal(6);

    bendingMoments(e) = max(abs([momentA, momentB]));
    shearForces(e) = max(abs([shearForce, -shearForce]));
    shearStresses(e) = shearForce / As;

    bendingPlaneValues = [ -momentA, -momentB; momentA, momentB ] / sectionModulus;
    bendingEnvelope = max(abs(bendingPlaneValues(:)));

    axialStress = axialForce / A;
    vonMises(e) = sqrt(axialStress^2 + bendingEnvelope^2 + 3 * shearStresses(e)^2);
    elemental(e) = vonMises(e);

    elementPositions(e, :) = mean(coords, 1);

    bendingDetails = updateBendingMetadata(bendingDetails, e, bendingPlaneValues, [momentA, momentB]);
    shearDetails   = updateShearMetadata(shearDetails, e, shearForce, shearStresses(e));
    vonMisesDetails = updateVonMisesMetadata(vonMisesDetails, e, vonMises(e), bendingPlaneValues, shearStresses(e));
  end

  performSectionModulusCheck(sectionType, width, height, diameter, A, sectionModulus);
end
```

### 4.7 `src/reporting/displayResults.m`
- **Responsibility**: Provide a quick displacement overview.
- **Output**: formatted console text summarising maximum horizontal displacement (`ux`), vertical displacement (`uy`), rotation (`theta`), and resultant magnitude per node, with qualitative guidance (small/moderate/large).
- **Implementation Detail**: loops operate on DOF indices using modular arithmetic (`mod(i, 3)` style) to separate translational and rotational components.

**MATLAB-style pseudocode**
```matlab
function displayResults(displacements)
  numNodes = numel(displacements) / 3;

  ux = displacements(1:3:end);
  uy = displacements(2:3:end);
  theta = displacements(3:3:end);

  mag = hypot(ux, uy);

  [maxUx, idxUx] = max(abs(ux));
  [maxUy, idxUy] = max(abs(uy));
  [maxTheta, idxTheta] = max(abs(theta));
  [maxMag, idxMag] = max(mag);

  fprintf('Max |ux| = %.6e m at node %d\n', maxUx, idxUx);
  fprintf('Max |uy| = %.6e m at node %d\n', maxUy, idxUy);
  fprintf('Max |theta| = %.6e rad at node %d\n', maxTheta, idxTheta);
  fprintf('Max |u| = %.6e m at node %d\n', maxMag, idxMag);

  if maxMag > 1e-2
    fprintf('Large displacements detected; review assumptions.\n');
  elseif maxMag > 1e-4
    fprintf('Moderate displacements observed.\n');
  else
    fprintf('Displacements are small relative to beam length.\n');
  end
end
```

### 4.8 `src/visualization/plotBeamSystem.m`
- **Responsibility**: Render the undeformed model for validation before solving.
- **Features**:
  - Element connectivity displayed as polylines.
  - Supports coloured by type (e.g., fixed/green, pinned/red, roller/blue).
  - Forces drawn as quiver arrows scaled heuristically; nodal moments indicated by circular arrows.
  - Data tips expose node IDs and coordinates via a custom callback.
- **Usage**: called automatically by the driver; can be invoked independently for debugging workbook layouts.

**MATLAB-style pseudocode**
```matlab
function plotBeamSystem(nodes, elements, supports, forces)
  figure('Name', 'Undeformed Beam System');
  axis equal; grid on; hold on;

  nodeCoords = containers.Map(nodes(:,1), num2cell(nodes(:,2:3), 2));

  for e = 1:size(elements, 1)
    n1 = elements(e, 2);
    n2 = elements(e, 3);
    p1 = nodeCoords(n1);
    p2 = nodeCoords(n2);
    plot([p1(1), p2(1)], [p1(2), p2(2)], 'k-', 'LineWidth', 2);
  end

  scatter(nodes(:,2), nodes(:,3), 40, 'k', 'filled');

  supportColours = struct('fixed', 'g', 'pinned', 'r', 'roller', 'b');
  for k = 1:numel(supports.NodeID)
    nodeId = supports.NodeID(k);
    pt = nodeCoords(nodeId);
    t = supports.Type(k);
    colour = supportColours.(char(t));
    scatter(pt(1), pt(2), 80, colour, 'filled');
  end

  scale = computeForceScale(forces(:,2:3));
  quiver(nodes(:,2), nodes(:,3), scale*forces(:,2), scale*forces(:,3), 0, 'm', 'LineWidth', 1.5);

  drawMoments(nodes, forces(:,4));
  title('Undeformed geometry, supports, and applied loads');
  xlabel('X [m]'); ylabel('Y [m]');
  hold off;
end
```

### 4.9 `src/visualization/plotDeformedShape.m`
- **Responsibility**: Overlay undeformed and deformed shapes, colouring the deformed geometry by displacement magnitude.
- **Workflow**:
  1. Compute displacement magnitudes per node.
  2. Select a scaling factor (default 20×, capped to avoid nonphysical visuals).
  3. Plot undeformed geometry in grey, then the deformed curve with colour gradations provided by `getColorFromValue`.
  4. Annotate the maximum displacement value and display a colour bar.
- **Extensibility**: adjust scaling heuristics or colour schemes to suit publication standards.

**MATLAB-style pseudocode**
```matlab
function plotDeformedShape(nodes, elements, displacements)
  ux = displacements(1:3:end);
  uy = displacements(2:3:end);
  nodeIds = nodes(:,1);
  nodeUx = containers.Map(nodeIds, num2cell(ux));
  nodeUy = containers.Map(nodeIds, num2cell(uy));

  nodeCoords = containers.Map(nodeIds, num2cell(nodes(:,2:3), 2));

  mags = hypot(ux, uy);
  maxMag = max(mags);
  scale = pickScaleFactor(maxMag);

  figure('Name', 'Deformed Beam Shape');
  hold on; grid on; axis equal;

  for e = 1:size(elements, 1)
    n1 = elements(e,2);
    n2 = elements(e,3);
    p1 = nodeCoords(n1);
    p2 = nodeCoords(n2);
    plot([p1(1), p2(1)], [p1(2), p2(2)], 'Color', 0.8*[1 1 1]);
  end

  cmap = jet(256);

  for e = 1:size(elements, 1)
    n1 = elements(e,2);
    n2 = elements(e,3);
    p1 = nodeCoords(n1);
    p2 = nodeCoords(n2);
    d1 = [nodeUx(n1); nodeUy(n1)];
    d2 = [nodeUx(n2); nodeUy(n2)];

    q1 = p1 + scale * d1.';
    q2 = p2 + scale * d2.';

    color1 = getColorFromValue(norm(d1), maxMag);
    color2 = getColorFromValue(norm(d2), maxMag);

    plot([q1(1), q2(1)], [q1(2), q2(2)], '-', 'Color', mean([color1; color2], 1), 'LineWidth', 2);
  end

  scatter(nodes(:,2) + scale * ux, nodes(:,3) + scale * uy, 36, mags, 'filled');
  colormap(cmap);
  colorbar;
  title(sprintf('Deformed shape (scale %.1f), max displacement %.3e m', scale, maxMag));
  xlabel('X [m]'); ylabel('Y [m]');
  hold off;
end
```

### 4.10 `src/visualization/getColorFromValue.m`
- **Responsibility**: Map a scalar value to an RGB triple using MATLAB's `jet` colormap.
- **Behaviour**: normalises by the supplied maximum, clamps indices to [1, 256], and returns the matching row from the colormap matrix. Handles zero-displacement cases gracefully.

**MATLAB-style pseudocode**
```matlab
function color = getColorFromValue(value, maxValue)
  if maxValue <= 0
    normalized = 0;
  else
    normalized = min(max(value / maxValue, 0), 1);
  end

  cmap = jet(256);
  idx = max(1, min(256, round(normalized * 255) + 1));
  color = cmap(idx, :);
end
```

### 4.11 `src/reporting/prepareStressReportEntries.m`
- **Responsibility**: Convert raw metadata into human-readable report lines.
- **Operation**:
  - Calls specialised summary helpers for bending, shear, and von Mises extrema.
  - Builds structs containing magnitude, signed value, element ID, plane label, end label, and formatted text.
  - Returns a struct that the driver iterates to print consistent console output.
- **Extension**: integrate additional limit states by writing new summary helpers and extending this formatter.

**MATLAB-style pseudocode**
```matlab
function entries = prepareStressReportEntries(stresses)
  bending = computeBendingExtremaSummary(stresses);
  shear = computeShearExtremaSummary(stresses);
  vm = computeVonMisesExtremaSummary(stresses);

  entries = struct([]);

  entries(end+1).key = 'maxBendingMoment';
  entries(end).value = bending.moment.magnitude;
  entries(end).message = sprintf('Max bending moment %.3e N·m at element %d (%s, %s)', ...
    bending.moment.magnitude, bending.moment.elementID, bending.moment.planeLabel, bending.moment.endLabel);

  for planeIdx = 1:numel(bending.stressPlane)
    plane = bending.stressPlane(planeIdx);
    entries(end+1).key = sprintf('maxBendingStressPlane%d', planeIdx);
    entries(end).value = plane.magnitude;
    entries(end).message = sprintf('Max bending stress (plane %s) %.3e Pa at element %d (%s)', ...
      plane.planeLabel, plane.magnitude, plane.elementID, plane.endLabel);
  end

  entries(end+1).key = 'maxBendingStressEnvelope';
  entries(end).value = bending.stressEnvelope.magnitude;
  entries(end).message = sprintf('Max bending stress envelope %.3e Pa at element %d (%s)', ...
    bending.stressEnvelope.magnitude, bending.stressEnvelope.elementID, bending.stressEnvelope.endLabel);

  entries(end+1).key = 'maxShearForce';
  entries(end).value = shear.force.magnitude;
  entries(end).message = sprintf('Max shear force %.3e N at element %d (%s)', ...
    shear.force.magnitude, shear.force.elementID, shear.force.endLabel);

  entries(end+1).key = 'maxShearStress';
  entries(end).value = shear.stress.magnitude;
  entries(end).message = sprintf('Max shear stress %.3e Pa at element %d (%s)', ...
    shear.stress.magnitude, shear.stress.elementID, shear.stress.endLabel);

  entries(end+1).key = 'maxVonMisesStress';
  entries(end).value = vm.value.magnitude;
  entries(end).message = sprintf('Max von Mises stress %.3e Pa at element %d (%s)', ...
    vm.value.magnitude, vm.elementID, vm.endLabel);
end
```

### 4.12 `src/reporting/computeBendingExtremaSummary.m`
- **Responsibility**: Identify the largest bending moment and stresses by plane (typically top/bottom fibres) and envelope.
- **Data Sources**: reads `stresses.bending_details`, which stores plane/end matrices and precomputed absolute maxima.
- **Edge Cases**: returns zeroed defaults when the stresses struct is empty, sustaining robust reporting for degenerate models.

**MATLAB-style pseudocode**
```matlab
function summary = computeBendingExtremaSummary(stresses)
  summary = makeDefaultBendingSummary();
  details = stresses.bending_details;

  if isempty(details) || isempty(details.moment.maxAbs)
    return;
  end

  [summary.moment.magnitude, idxMoment] = max(details.moment.maxAbs);
  summary.moment.signedValue = details.moment.signed(idxMoment);
  summary.moment.elementID = stresses.element_ids(idxMoment);
  summary.moment.planeIdx = details.moment.planeIdx(idxMoment);
  summary.moment.planeLabel = details.planeLabels(summary.moment.planeIdx);
  summary.moment.endIdx = details.moment.endIdx(idxMoment);
  summary.moment.endLabel = details.endLabels(summary.moment.endIdx);

  for planeIdx = 1:numel(details.planeLabels)
    [mag, idxPlane] = max(details.stress.maxAbs(:, planeIdx));
    summary.stressPlane(planeIdx).magnitude = mag;
    summary.stressPlane(planeIdx).signedValue = details.stress.signed(idxPlane, planeIdx);
    summary.stressPlane(planeIdx).elementID = stresses.element_ids(idxPlane);
    summary.stressPlane(planeIdx).planeIdx = planeIdx;
    summary.stressPlane(planeIdx).planeLabel = details.planeLabels(planeIdx);
    summary.stressPlane(planeIdx).endIdx = details.stress.endIdx(idxPlane, planeIdx);
    summary.stressPlane(planeIdx).endLabel = details.endLabels(summary.stressPlane(planeIdx).endIdx);
  end

  [summary.stressEnvelope.magnitude, idxEnv] = max(details.stressEnvelope.maxAbs(:));
  summary.stressEnvelope.signedValue = details.stressEnvelope.signed(idxEnv);
  summary.stressEnvelope.elementID = stresses.element_ids(details.stressEnvelope.elementIdx(idxEnv));
  summary.stressEnvelope.planeIdx = details.stressEnvelope.planeIdx(idxEnv);
  summary.stressEnvelope.planeLabel = details.planeLabels(summary.stressEnvelope.planeIdx);
  summary.stressEnvelope.endIdx = details.stressEnvelope.endIdx(idxEnv);
  summary.stressEnvelope.endLabel = details.endLabels(summary.stressEnvelope.endIdx);
end
```

### 4.13 `src/reporting/computeShearExtremaSummary.m`
- **Responsibility**: Extract controlling shear force and shear stress, preserving plane/end metadata (planes represent positive/negative faces for shear stress bookkeeping).
- **Notable Logic**: clamps plane/end indices to valid ranges to protect against data corruption or coding errors upstream.

**MATLAB-style pseudocode**
```matlab
function summary = computeShearExtremaSummary(stresses)
  summary = makeDefaultShearSummary();
  details = stresses.shear_details;

  if isempty(details) || isempty(details.force.maxAbs)
    return;
  end

  [summary.force.magnitude, idxForce] = max(details.force.maxAbs);
  summary.force.signedValue = details.force.signed(idxForce);
  summary.force.elementID = stresses.element_ids(idxForce);
  summary.force.endIdx = clampIndex(details.force.endIdx(idxForce), numel(details.endLabels));
  summary.force.endLabel = details.endLabels(summary.force.endIdx);

  [summary.stress.magnitude, idxStress] = max(details.stress.maxAbs);
  summary.stress.signedValue = details.stress.signed(idxStress);
  summary.stress.elementID = stresses.element_ids(idxStress);
  summary.stress.endIdx = clampIndex(details.stress.endIdx(idxStress), numel(details.endLabels));
  summary.stress.endLabel = details.endLabels(summary.stress.endIdx);
end
```

### 4.14 `src/reporting/computeVonMisesExtremaSummary.m`
- **Responsibility**: Determine the maximum von Mises stress and its location metadata.
- **Usage**: feeds the console report and provides a quick check when comparing against allowable stresses in design workflows.

**MATLAB-style pseudocode**
```matlab
function summary = computeVonMisesExtremaSummary(stresses)
  summary = makeDefaultVonMisesSummary();
  details = stresses.von_mises_details;

  if isempty(details) || isempty(details.maxAbs)
    return;
  end

  [summary.value.magnitude, idx] = max(details.maxAbs);
  summary.value.signedValue = details.signed(idx);
  summary.elementID = stresses.element_ids(details.elementIdx(idx));
  summary.planeIdx = clampIndex(details.planeIdx(idx), numel(details.planeLabels));
  summary.planeLabel = details.planeLabels(summary.planeIdx);
  summary.endIdx = clampIndex(details.endIdx(idx), numel(details.endLabels));
  summary.endLabel = details.endLabels(summary.endIdx);
end
```

---

## 5. Algorithmic Notes

- **Degree-of-Freedom Ordering**: each node contributes `(ux, uy, theta)` to the global vector. Assembly functions map element DOFs using this ordering, so changes must propagate consistently.
- **Transformation Matrices**: element matrices are constructed in local coordinates then rotated to global via 6×6 transformation matrices built from direction cosines `(c, s)`.
- **Shear Deformation**: the correction factor `phi` modifies the bending portion of the stiffness matrix, distinguishing Timoshenko behaviour from Euler-Bernoulli. Setting `G` very large or `As → ∞` effectively recovers Euler-Bernoulli behaviour.
- **Self-Weight Loading**: each element's distributed weight (`ρ · A · L · g`) is converted into equivalent nodal forces and moments using a consistent formulation. These are superimposed on user-specified nodal loads.
- **Von Mises Approximation**: the combined stress metric uses `σ_vm = sqrt(σ_axial^2 + σ_bending^2 + 3 τ^2)`. This matches common beam design approximations where axial and bending stresses act along the same axis with shear through thickness.

---

## 6. Diagnostics and Error Handling

- **Workbook Errors**: missing sheets or headers raise MATLAB errors with guidance describing expected sheet names and columns.
- **Constraint Warnings**: `performTimoshenkoFEA` checks the condition number of `K_ff`. If the system is singular or poorly conditioned, a warning prompts users to review support definitions.
- **Geometry Sanity Checks**: `calculateTimoshenkoStresses` warns when user-supplied geometry deviates from implied section modulus values beyond ±20%, preventing silent mistakes in section dimensions.
- **Plotting Safeguards**: both plotting functions protect against empty datasets and handle zero displacement cases without throwing errors.

---

## 7. Extending the Toolkit

1. **New Section Types**: add cases to `calculateTimoshenkoSectionProperties`. Provide formulas for `I`, `ky`, and `maxDistance`, then update workbook documentation so users can supply the required dimensions.
2. **Element-Wise Properties**: modify `readInputData` to ingest property columns per element, and update the solver to read the appropriate value inside element loops.
3. **Distributed Loads**: extend `performTimoshenkoFEA` with functions that translate line loads into equivalent nodal forces, similar to the existing self-weight logic.
4. **Exporters**: create utilities under `src/reporting/` to write CSV/Excel summaries. Reuse `prepareStressReportEntries` to avoid duplicating formatting rules.
5. **Automation**: embed the solver within optimisation or parametric studies by calling `performTimoshenkoFEA` directly inside your scripts. Disconnect plotting to run headless analyses when performance matters.

---

## 8. Verification Strategy

- **Analytical Benchmarks**: start with classical cases (cantilever with tip load, simply supported beam with uniform load) and compare against published Timoshenko solutions for deflection and rotation.
- **Mesh Convergence**: refine element counts and observe convergence of displacements and stresses to ensure the formulation behaves as expected.
- **Regression Scripts**: capture benchmark load cases in simple MATLAB scripts (e.g., under a `tests/` folder) that assert tolerance-based comparisons. This guards against regressions when extending functionality.
- **Peer Review**: document load cases, assumptions, and results in `docs/` to maintain traceability, especially when using the toolkit in coursework or research publications.

---

## 9. Glossary
- **DOF**: Degree of Freedom (`ux`, `uy`, `theta` per node).
- **Element Ends**: In reporting, `End A` corresponds to the first node, `End B` to the second node of each element.
- **Planes**: Top/bottom fibres used for bending stress evaluation; positive plane typically represents the positive local y-face.
- **Ky**: Shear correction factor applied to cross-sectional area to account for non-uniform shear stress distribution.
- **Von Mises Stress**: Scalar combination of normal and shear stresses used for yield checks.

---

For additional context, refer back to `README.md` for high-level usage instructions or explore the source files directly to inspect implementation details.
