# Timoshenko Beam FEA Module Documentation

This document explains the purpose, data flow, and implementation details of the MATLAB modules that support the Timoshenko beam finite element analysis workflow. Every section includes a short summary, the expected inputs and outputs, and high-level pseudocode to clarify the computational steps.

---

## System Overview

- The entry point `scripts/TimoshenkoBeamFEA.m` orchestrates the analysis: it loads input data, runs the solver, prints numerical results, and generates visualizations.
- Raw model data is imported from an Excel workbook (`src/io/readInputData.m`).
- Section properties specific to Timoshenko theory are derived in `src/properties/calculateTimoshenkoSectionProperties.m`.
- Element stiffness, global assembly, and response calculations reside in `src/analysis`.
- Reporting utilities in `src/reporting` and `src/visualization` present the results.

Data moves through the system in the following order:

```
Excel workbook
  -> readInputData
    -> TimoshenkoBeamFEA (driver)
      -> plotBeamSystem (pre-analysis visualization)
      -> performTimoshenkoFEA (solver)
          -> calculateTimoshenkoSectionProperties
          -> calculateTimoshenkoStiffness (per element)
          -> calculateTimoshenkoStresses (postprocessing)
      -> displayResults (displacement summary)
      -> plotDeformedShape + getColorFromValue (post-analysis visualization)
```

---

## scripts/TimoshenkoBeamFEA.m

**Role**: Main workflow function users call to execute a complete Timoshenko beam analysis.

**Inputs**
- `inputWorkbook` (string, optional): Path to the Excel workbook containing model data. Defaults to `data/input_data.xlsx` inside the project.

**Outputs**
- None (prints results and generates plots). Returns when analysis completes.

**Key operations**
- Ensures the `src` directory is on the MATLAB path.
- Falls back to the packaged sample workbook when no argument is supplied.
- Delegates to the solver and visualization utilities in sequence.

**Pseudocode**

```text
function TimoshenkoBeamFEA(inputWorkbook)
    add project src directory to MATLAB path
    if inputWorkbook missing
        set inputWorkbook to default sample workbook
    [nodes, elements, supports, forces, material props] = readInputData(inputWorkbook)
    plotBeamSystem(nodes, elements, supports, forces)
    [displacements, stresses] = performTimoshenkoFEA(nodes, elements, supports, forces, props...)
    displayResults(displacements)
    plotDeformedShape(nodes, elements, displacements)
    print stress extrema (elemental, bending, shear)
end
```

**Notable considerations**
- Designed as a function, so it can be invoked from MATLAB scripts or the command window.
- Adds the entire `src` tree to the path to support modular development.

---

## src/io/readInputData.m

**Role**: Load model definition, boundary conditions, and material properties from an Excel workbook.

**Inputs**
- `inputWorkbook` (string, optional): Path to Excel workbook. Defaults to `input_data.xlsx` in the current working directory.

**Outputs**
- `nodes` (`n×3` numeric array): `[NodeID, X, Y]`.
- `elements` (`e×3` numeric array): `[ElementID, Node1, Node2]`.
- `supports` (struct): Fields `NodeID` (vector), `Type` (string array).
- `forces` (`f×4` numeric array): `[NodeID, Fx, Fy, Mz]` with zero-filled `Mz` if absent.
- `E`, `A`, `density`, `sectionType`, `width`, `height`, `diameter`, `G`, `nu`: Material and geometric scalars.

**Pseudocode**

```text
function readInputData(inputWorkbook)
    if workbook path missing -> default name
    if file not found -> error
    try
        read tables from Nodes, Elements, Supports, Forces, Properties sheets
        standardize column names (NodeID, Fx, etc.)
        ensure required force columns exist; add zero Mz if missing
        extract basic material properties (E, A, density)
        initialize default section geometry (rectangle assumptions)
        override defaults with Properties sheet values when present
        compute shear modulus G from sheet or via E and nu
        convert tables to numeric arrays / struct for downstream code
        print summary of loaded data
    catch exception
        print expected workbook structure and rethrow
    end
    return loaded arrays and scalars
end
```

**Notable considerations**
- Accepts minor variations in column naming by normalizing identifiers.
- Gracefully computes shear modulus from Poisson's ratio when necessary.
- Keeps support definitions as a struct because they mix text and numeric fields.

---

## src/properties/calculateTimoshenkoSectionProperties.m

**Role**: Derive section properties (moment of inertia, effective shear area, etc.) tailored to Timoshenko beam theory based on the cross-section shape.

**Inputs**
- `sectionType` (string): `"Rectangle"`, `"Square"`, `"Circle"`, or other.
- `A` (scalar): Cross-sectional area.
- `width`, `height`, `diameter` (scalars): Shape dimensions; unused values can be zero.

**Outputs**
- `I`: Second moment of area about the neutral axis.
- `As`: Effective shear area (`ky * A`).
- `maxDistance`: Distance from neutral axis to the outer fiber (for bending stress).
- `ky`: Shear correction factor.

**Pseudocode**

```text
function [I, As, maxDistance, ky] = calculateTimoshenkoSectionProperties(type, A, width, height, diameter)
    switch lowercase(type)
        case 'rectangle'
            I = width * height^3 / 12
            maxDistance = height / 2
            ky = 5/6
        case 'square'
            I = width^4 / 12
            maxDistance = width / 2
            ky = 5/6
        case 'circle'
            I = pi * diameter^4 / 64
            maxDistance = diameter / 2
            ky = 9/10
        otherwise
            I = A^2 / 12
            maxDistance = sqrt(A) / 2
            ky = 5/6
            warn about unknown section type
    end
    As = ky * A
end
```

**Notable considerations**
- Defaults to conservative rectangular assumptions when the section type is unrecognized.
- Returns `maxDistance` for downstream stress calculations (bending stress = `M * maxDistance / I`).

---

## src/analysis/calculateTimoshenkoStiffness.m

**Role**: Build the 6×6 local stiffness matrix for a single Timoshenko beam element.

**Inputs**
- `E`, `G`: Young's modulus and shear modulus.
- `A`, `As`: Cross-sectional area and effective shear area.
- `I`: Second moment of area.
- `L`: Element length.

**Output**
- `k_local` (`6×6` matrix): Element stiffness in local coordinates with DOF ordering `[u1, v1, theta1, u2, v2, theta2]`.

**Pseudocode**

```text
function k_local = calculateTimoshenkoStiffness(E, G, A, As, I, L)
    phi = 12 * E * I / (G * As * L^2)  // shear flexibility correction
    compute axial stiffness terms (k11, k14, ...)
    compute bending/shear coupling terms scaled by phi
    compute rotational stiffness terms (k33, k36, ...)
    assemble the symmetric 6x6 matrix in standard DOF order
    return k_local
end
```

**Notable considerations**
- When `phi` approaches zero, the terms converge to Euler-Bernoulli stiffness.
- Supports asymmetric force-rotation coupling terms typical for Timoshenko theory.

---

## src/analysis/performTimoshenkoFEA.m

**Role**: Solve the global finite element problem: assemble the stiffness matrix, apply boundary conditions, solve for nodal displacements, and compute derived stresses.

**Inputs**
- Structural definition arrays (`nodes`, `elements`, `supports`, `forces`).
- Material and section properties (`E`, `A`, `density`, `sectionType`, `width`, `height`, `diameter`, `G`, `nu`).

**Outputs**
- `displacements` (`3n×1` vector): Nodal displacements ordered `[u1, v1, theta1, u2, ...]`.
- `stresses` (struct): Fields `elemental`, `bending_moments`, `shear_forces`, `shear_stresses`, `element_positions`.

**Pseudocode**

```text
function [displacements, stresses] = performTimoshenkoFEA(nodes, elements, supports, forces, material props...)
    num_nodes = count rows in nodes
    num_elements = count rows in elements
    total_dof = 3 * num_nodes
    initialize global stiffness K (total_dof x total_dof) and force vector F

    [I, As, maxDistance, ky] = calculateTimoshenkoSectionProperties(sectionType, A, width, height, diameter)
    print summary of model and material properties

    for each element in elements
        locate node rows for connectivity
        compute element length L and orientation cosines (c, s)
        k_local = calculateTimoshenkoStiffness(E, G, A, As, I, L)
        build transformation matrix T from orientation cosines
        k_global = T' * k_local * T
        scatter-add k_global into the appropriate entries of K
    end

    for each applied nodal force
        map Fx, Fy, and optional Mz into F based on node row
    end

    for each element
        compute self-weight = density * A * L * g
        distribute half the weight to each node's vertical DOF in F
    end

    identify constrained DOFs from supports (fixed, pinned, roller)
    free_dof = complement of constrained set

    extract Kff = K(free_dof, free_dof), Ff = F(free_dof)
    if rcond(Kff) very small -> warn about ill-conditioning
    solve displacements on free DOF: displacements(free_dof) = Kff \ Ff

    [elemental_stresses, bending_moments, shear_forces, shear_stresses, element_positions] =
        calculateTimoshenkoStresses(nodes, elements, displacements, E, G, A, As, I, maxDistance)

    pack stress results into struct and return
end
```

**Notable considerations**
- Assumes 2D planar frame elements with DOF order `[u, v, theta]` at each node.
- Self-weight is included as a uniformly distributed load converted to equivalent nodal forces.
- Supports can be mixed; the logic simply appends constrained DOFs before solving.

---

## src/analysis/calculateTimoshenkoStresses.m

**Role**: Post-process nodal displacements to compute element-level axial forces, bending moments, shear forces, and derived stress measures.

**Inputs**
- `nodes`, `elements`, `displacements`: Structural definition and solved displacements.
- Material/section data: `E`, `G`, `A`, `As`, `I`, `maxDistance`.

**Outputs**
- `elemental_stresses`: Scalar von-Mises-like stress per element.
- `bending_moments`: Internal bending moments.
- `shear_forces`: Internal shear forces.
- `shear_stresses`: Shear stress per element (`V / As`).
- `element_positions`: Midpoint x-coordinates (for plotting or reporting).

**Pseudocode**

```text
function calculateTimoshenkoStresses(nodes, elements, displacements, E, G, A, As, I, maxDistance)
    initialize result vectors sized to number of elements
    for each element
        locate node rows for its end nodes
        extract local DOF displacements (u, v, theta at each node)
        compute element length L and direction cos/sin
        build transformation matrix T and form local displacement vector d_local = T * [d1; d2]

        axial_strain = (d_local(4) - d_local(1)) / L
        curvature = (d_local(6) - d_local(3)) / L
        theta_avg = (d_local(3) + d_local(6)) / 2
        shear_strain = (d_local(5) - d_local(2)) / L - theta_avg

        axial_force = E * A * axial_strain
        bending_moment = E * I * curvature
        shear_force = G * As * shear_strain

        axial_stress = axial_force / A
        bending_stress = bending_moment * maxDistance / I
        shear_stress = shear_force / As

        elemental_stress = sqrt(axial_stress^2 + bending_stress^2 + 3 * shear_stress^2)
        store outputs and element midpoint position
    end
    return all result vectors
end
```

**Notable considerations**
- Uses a von-Mises-style combination of axial, bending, and shear stresses to produce a single comparison metric.
- Transformation matrix matches the one used during stiffness assembly to ensure consistent local coordinates.

---

## src/reporting/displayResults.m

**Role**: Print key displacement metrics to the MATLAB console for quick inspection.

**Inputs**
- `displacements` (`3n×1` vector): Global nodal solution from the solver.

**Outputs**
- None (produces console output).

**Pseudocode**

```text
function displayResults(displacements)
    max_x = max(abs(displacements at DOF index 1 mod 3))
    max_y = max(abs(displacements at DOF index 2 mod 3))
    max_theta = max(abs(displacements at DOF index 0 mod 3))

    for each node
        compute magnitude = sqrt(ux^2 + uy^2)
    end

    print formatted summary (max horizontal, vertical, rotation, total magnitude)
    print qualitative note (large/moderate/small) based on magnitude threshold
end
```

**Notable considerations**
- Provides basic interpretive feedback (e.g., "Large displacements detected") to guide users.

---

## src/visualization/plotBeamSystem.m

**Role**: Produce a pre-analysis plot of the undeformed beam system, highlighting supports and applied loads.

**Inputs**
- `nodes`, `elements`: Geometry definition.
- `supports`: Struct produced by `readInputData`.
- `forces`: Numeric array `[NodeID, Fx, Fy, Mz]`.

**Outputs**
- MATLAB figure with interactive data tips exposing node IDs and coordinates.

**Pseudocode**

```text
function plotBeamSystem(nodes, elements, supports, forces)
    create figure and hold on
    for each element
        draw a thick line between its node coordinates
    end
    for each node
        draw a circular marker; if supported, override marker color based on type
    end
    for each force entry
        draw scaled quiver arrows for Fx and Fy
        if Mz present, draw circular arrow indicating moment direction
    end
    set axis equal, enable grid, label axes, add title
    configure custom data cursor callback to show NodeID and coordinates
    hold off
end
```

**Notable considerations**
- Uses color coding (green fixed, red pinned, blue roller) to distinguish support types.
- Plots nodal moments as magenta circular arrows with simple arrowheads.

---

## src/visualization/plotDeformedShape.m

**Role**: Visualize the post-analysis deformed configuration with a color gradient mapping displacement magnitude.

**Inputs**
- `nodes`, `elements`, `displacements`: Structural geometry and solved nodal response.

**Outputs**
- MATLAB figure illustrating original versus deformed geometry and displacement magnitudes.

**Pseudocode**

```text
function plotDeformedShape(nodes, elements, displacements)
    create figure and hold on
    compute displacement magnitude at each node
    determine visualization scale factor (default 20, capped by max displacement)

    for each element
        plot original element in light gray
    end

    configure jet colormap, colorbar, and caxis based on max displacement

    for each element
        extract original coordinates and displacements for both nodes
        compute deformed coordinates = original + scale_factor * displacement
        interpolate intermediate points and displacement magnitudes along element
        draw colored line segments using getColorFromValue for gradient effect
        scatter deformed nodes with color mapped to magnitude
    end

    finalize axes, title, legend, and annotation (max displacement info)
    hold off
end
```

**Notable considerations**
- Scales displacements to remain visually interpretable regardless of magnitude.
- Annotates the plot with maximum displacement and explains that gray lines show the undeformed reference.

---

## src/visualization/getColorFromValue.m

**Role**: Helper that maps a scalar displacement magnitude to an RGB triplet using MATLAB's `jet` colormap.

**Inputs**
- `value`: Scalar (typically displacement magnitude at a point).
- `max_value`: Maximum displacement magnitude (normalization denominator).

**Output**
- `color`: RGB row vector with components in `[0, 1]`.

**Pseudocode**

```text
function color = getColorFromValue(value, max_value)
    if max_value == 0
        normalized = 0
    else
        normalized = value / max_value
    end
    cmap = jet(256)
    index = clamp(round(normalized * 255) + 1, 1, 256)
    color = cmap(index, :)
end
```

**Notable considerations**
- Gracefully handles the degenerate case where the entire model has zero displacement.
- Uses 256 bins to match MATLAB's default size for the `jet` colormap and ensure smooth gradients.

---

## Extending the Module

- **Additional loading types**: The solver currently treats self-weight as a uniformly distributed load. Line or point loads on elements could be added by expanding the force assembly step in `performTimoshenkoFEA`.
- **Material variability**: To support element-wise material properties, modify `readInputData` to ingest per-element data and adjust `performTimoshenkoFEA` to look up properties per element.
- **Result export**: Implement CSV or Excel writers in `src/reporting` to export nodal displacements and element forces for external post-processing.
- **3D extension**: Transitioning to 3D Timoshenko or Timoshenko-Timoshenko frame elements would require expanding the DOF set and updating transformation matrices accordingly.
