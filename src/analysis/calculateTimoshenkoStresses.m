function [elemental_stresses, bending_moments, shear_forces, shear_stresses, element_positions, bending_stress_planes, bending_details, shear_details, von_mises_details] = ...
    calculateTimoshenkoStresses(nodes, elements, displacements, E, G, A, As, I, maxDistance, sectionType, width, height, diameter, density)
    % CALCULATE STRESSES AND INTERNAL FORCES IN TIMOSHENKO BEAM ELEMENTS
    % This function takes the calculated displacements and determines the
    % internal forces and stresses that develop in each beam element.
    % This is crucial for checking if the beam is safe and won't fail.
    %
    % The key difference from Euler-Bernoulli theory is that Timoshenko
    % theory accounts for shear deformation, giving more accurate results
    % especially for thick beams or beams with low shear stiffness.
    
    % Initialize arrays to store results for each element
    num_elements = size(elements, 1);
    elemental_stresses = zeros(num_elements, 1);  % Combined stress measure
    bending_moments = zeros(num_elements, 1);     % Max |bending moment| envelope (Plane/End)
    shear_forces = zeros(num_elements, 1);        % Internal shear forces (signed at controlling location)
    shear_stresses = zeros(num_elements, 1);      % Shear stresses (signed at controlling location)
    element_positions = zeros(num_elements, 1);   % Position for plotting
    bending_stress_planes = zeros(num_elements, 2); % Extreme-fiber stresses per plane (|Plane1|, |Plane2|)

    % Bookkeeping arrays for bending moment/stress plane-end tracking
    bending_moment_plane_end = zeros(num_elements, 2, 2);
    bending_moment_plane_indices = zeros(num_elements, 1);
    bending_moment_end_indices = zeros(num_elements, 1);
    bending_moment_signed_values = zeros(num_elements, 1);
    bending_moment_abs_values = zeros(num_elements, 1);

    bending_stress_plane_end = zeros(num_elements, 2, 2);
    bending_stress_plane_max_end_idx = zeros(num_elements, 2);
    bending_stress_plane_max_signed = zeros(num_elements, 2);
    bending_stress_plane_max_abs = zeros(num_elements, 2);
    bending_stress_max_plane_indices = zeros(num_elements, 1);
    bending_stress_max_end_indices = zeros(num_elements, 1);
    bending_stress_max_signed = zeros(num_elements, 1);
    bending_stress_max_abs = zeros(num_elements, 1);

    shear_plane_end_forces = zeros(num_elements, 2, 2);  % Plane/End shear force matrix per element
    shear_plane_end_stresses = zeros(num_elements, 2, 2);% Plane/End shear stress matrix per element
    shear_plane_indices = zeros(num_elements, 1);
    shear_end_indices = zeros(num_elements, 1);
    shear_abs_values = zeros(num_elements, 1);

    von_mises_plane_end = zeros(num_elements, 2, 2);
    von_mises_max_plane_indices = zeros(num_elements, 1);
    von_mises_max_end_indices = zeros(num_elements, 1);
    von_mises_max_abs = zeros(num_elements, 1);

    % Pre-compute section modulus and validate geometry inputs
    if I <= 0 || maxDistance <= 0
        error('Invalid section properties: I and maxDistance must be positive.');
    end
    section_modulus = I / maxDistance;

    % Runtime geometry sanity checks to catch inconsistent inputs
    s_expected = NaN;
    shape = lower(char(sectionType));
    switch shape
        case 'rectangle'
            if width > 0 && height > 0
                s_expected = width * height^2 / 6;
            end
        case 'square'
            side = width;
            if side <= 0 && height > 0
                side = height;
            end
            if side > 0
                s_expected = side^3 / 6;
            end
        case 'circle'
            if diameter <= 0 && width > 0
                diameter = width;
            end
            if diameter > 0
                s_expected = pi * diameter^3 / 32;
            end
    end

    if ~isnan(s_expected) && s_expected > 0
        if abs(section_modulus / s_expected - 1) > 0.2
            switch shape
                case 'rectangle'
                    warning('Section modulus (rectangle) off by >20%. Check units and geometry.');
                case 'square'
                    warning('Section modulus (square) off by >20%. Check units and geometry.');
                case 'circle'
                    warning('Section modulus (circle) off by >20%. Check units and geometry.');
            end
        end
    end
    
    % Calculate stresses for each element
    for i = 1:num_elements
        % Find the nodes that this element connects
        node1_row = find(nodes(:,1) == elements(i,2), 1);
        node2_row = find(nodes(:,1) == elements(i,3), 1);
        
        % Extract the displacements for both nodes
        % Each node has 3 DOFs: u (axial), v (transverse), θ (rotation)
        d1 = displacements(3*node1_row-2:3*node1_row);
        d2 = displacements(3*node2_row-2:3*node2_row);
        
        % Calculate element geometry
        x1 = nodes(node1_row,2); y1 = nodes(node1_row,3);
        x2 = nodes(node2_row,2); y2 = nodes(node2_row,3);
        L = sqrt((x2-x1)^2 + (y2-y1)^2);  % Element length
        
        % Store element center position for result plotting
        element_positions(i) = (x1 + x2) / 2;
        
        % Transform displacements to local element coordinates
        % Each element has its own local coordinate system
        c = (x2-x1)/L;  % Cosine of element orientation
        s = (y2-y1)/L;  % Sine of element orientation
        T = [c s 0 0 0 0;
            -s c 0 0 0 0;
             0 0 1 0 0 0;
             0 0 0 c s 0;
             0 0 0 -s c 0;
             0 0 0 0 0 1];
        
        % Get displacements in local coordinates
        d_local = T * [d1; d2];
        
        % CALCULATE STRAINS FROM DISPLACEMENTS
        % These represent the deformation of the element
        
        % Calculate axial strain (how much the element stretches or compresses)
        axial_strain = (d_local(4) - d_local(1))/L;
        
        % Calculate shear strain (unique to Timoshenko theory)
        % This represents the "sliding" deformation between node ends
        % γ = (dv/dx - θ) where θ is the rotation of the cross-section
        % CALCULATE INTERNAL FORCES FROM STRAINS
        % These are the forces that the material develops to resist deformation
        
        % Axial force (tension or compression along the beam length)
        axial_force = E * A * axial_strain;
        
        % Extract local end moments to build an A/B envelope
        k_local = calculateTimoshenkoStiffness(E, G, A, As, I, L);
        internal_forces_local = k_local * d_local;

        % Remove equivalent nodal loads from self-weight to recover true internal forces
        f_equiv_local = zeros(6,1);
        grav = 9.81;
        tol = 1e-12;
        if density > 0
            w_global = [0; -density * A * grav];
            qx_local = c * w_global(1) + s * w_global(2);
            qy_local = -s * w_global(1) + c * w_global(2);

            if abs(qx_local) > tol
                axial_eq = qx_local * L / 2;
                f_equiv_local(1) = f_equiv_local(1) + axial_eq;
                f_equiv_local(4) = f_equiv_local(4) + axial_eq;
            end

            if abs(qy_local) > tol
                shear_eq = qy_local * L / 2;
                moment_eq = qy_local * L^2 / 12;
                f_equiv_local(2) = f_equiv_local(2) + shear_eq;
                f_equiv_local(5) = f_equiv_local(5) + shear_eq;
                f_equiv_local(3) = f_equiv_local(3) + moment_eq;
                f_equiv_local(6) = f_equiv_local(6) - moment_eq;
            end
        end

        internal_forces_local = internal_forces_local - f_equiv_local;

        moment_endA = internal_forces_local(3);
        moment_endB = internal_forces_local(6);
        moment_plane_end = [moment_endA, moment_endB; -moment_endA, -moment_endB];
        bending_moment_plane_end(i,:,:) = moment_plane_end;
        bending_moments(i) = computePlaneEndEnvelope(moment_plane_end);
        [~, linearIdxMoment] = max(abs(moment_plane_end), [], 'all', 'linear');
        [plane_idx_m, end_idx_m] = ind2sub(size(moment_plane_end), linearIdxMoment);
        bending_moment_plane_indices(i) = plane_idx_m;
        bending_moment_end_indices(i) = end_idx_m;
        bending_moment_signed_values(i) = moment_plane_end(plane_idx_m, end_idx_m);
        bending_moment_abs_values(i) = abs(bending_moment_signed_values(i));

        % Shear force (causes sliding between layers of the beam)
        shear_endA = internal_forces_local(2);
        shear_endB = internal_forces_local(5);
        shear_plane_end = [shear_endA, shear_endB; -shear_endA, -shear_endB];
        shear_plane_end_forces(i,:,:) = shear_plane_end;

        [shear_selected, plane_idx, end_idx, shear_abs] = extractShearEnvelope(shear_plane_end);
        shear_plane_indices(i) = plane_idx;
        shear_end_indices(i) = end_idx;
        shear_forces(i) = shear_selected;
        shear_plane_end_stresses(i,:,:) = shear_plane_end / As;
        shear_stresses(i) = shear_selected / As;
        shear_abs_values(i) = shear_abs;
        
        % CALCULATE STRESSES FROM FORCES
        % Stresses are forces per unit area - these determine if the material will fail
        
        % Axial stress (normal stress along the beam length)
        axial_stress = axial_force / A;
        
        % Extreme-fiber bending stresses (Plane 1 = top, Plane 2 = bottom)
        sigma_plane1 = -[moment_endA, moment_endB] / section_modulus;
        sigma_plane2 = [moment_endA, moment_endB] / section_modulus;
        plane_sigma = [sigma_plane1; sigma_plane2];
        bending_stress_planes(i,1) = computePlaneEndEnvelope(sigma_plane1);
        bending_stress_planes(i,2) = computePlaneEndEnvelope(sigma_plane2);
        bending_stress_plane_end(i,:,:) = plane_sigma;

        [plane1_abs, plane1_endIdx] = max(abs(sigma_plane1));
        bending_stress_plane_max_abs(i,1) = plane1_abs;
        bending_stress_plane_max_end_idx(i,1) = plane1_endIdx;
        bending_stress_plane_max_signed(i,1) = sigma_plane1(plane1_endIdx);

        [plane2_abs, plane2_endIdx] = max(abs(sigma_plane2));
        bending_stress_plane_max_abs(i,2) = plane2_abs;
        bending_stress_plane_max_end_idx(i,2) = plane2_endIdx;
        bending_stress_plane_max_signed(i,2) = sigma_plane2(plane2_endIdx);

        [~, linearIdxStress] = max(abs(plane_sigma), [], 'all', 'linear');
        [plane_idx_s, end_idx_s] = ind2sub(size(plane_sigma), linearIdxStress);
        bending_stress_max_plane_indices(i) = plane_idx_s;
        bending_stress_max_end_indices(i) = end_idx_s;
        bending_stress_max_signed(i) = plane_sigma(plane_idx_s, end_idx_s);
        bending_stress_max_abs(i) = abs(bending_stress_max_signed(i));

        bending_stress_extreme = computePlaneEndEnvelope(plane_sigma);

        shear_plane_end_stress = squeeze(shear_plane_end_stresses(i,:,:));
        von_mises_local = sqrt(plane_sigma.^2 + 3 * (shear_plane_end_stress.^2));
        von_mises_plane_end(i,:,:) = von_mises_local;
        [von_mises_max, linearIdxVM] = max(von_mises_local, [], 'all', 'linear');
        [plane_idx_vm, end_idx_vm] = ind2sub(size(von_mises_local), linearIdxVM);
        von_mises_max_plane_indices(i) = plane_idx_vm;
        von_mises_max_end_indices(i) = end_idx_vm;
        von_mises_max_abs(i) = von_mises_max;
        
        % Shear stress (varies across the cross-section, maximum at neutral axis)
        % Already stored as shear_stresses(i) above.

        % COMBINED STRESS MEASURE
        % Use von Mises equivalent stress as a measure of overall stress state
        % This combines all stress components into a single value for comparison
        % with material strength limits
        shear_component = abs(shear_stresses(i));
        elemental_stresses(i) = sqrt(axial_stress^2 + bending_stress_extreme^2 + 3 * shear_component^2);
    end

    bending_moment_details = struct('planeEndMoments', bending_moment_plane_end, ...
                                    'maxPlaneIdx', bending_moment_plane_indices, ...
                                    'maxEndIdx', bending_moment_end_indices, ...
                                    'maxSigned', bending_moment_signed_values, ...
                                    'maxAbs', bending_moment_abs_values);

    bending_stress_details = struct('planeEndStresses', bending_stress_plane_end, ...
                                    'planeMaxEndIdx', bending_stress_plane_max_end_idx, ...
                                    'planeMaxSigned', bending_stress_plane_max_signed, ...
                                    'planeMaxAbs', bending_stress_plane_max_abs, ...
                                    'maxPlaneIdx', bending_stress_max_plane_indices, ...
                                    'maxEndIdx', bending_stress_max_end_indices, ...
                                    'maxSigned', bending_stress_max_signed, ...
                                    'maxAbs', bending_stress_max_abs);

    bending_details = struct('moment', bending_moment_details, ...
                             'stress', bending_stress_details);

    shear_details = struct('planeEndForces', shear_plane_end_forces, ...
                           'planeEndStresses', shear_plane_end_stresses, ...
                           'maxPlaneIdx', shear_plane_indices, ...
                           'maxEndIdx', shear_end_indices, ...
                           'maxAbs', shear_abs_values(:));

    von_mises_details = struct('planeEndValues', von_mises_plane_end, ...
                               'maxPlaneIdx', von_mises_max_plane_indices, ...
                               'maxEndIdx', von_mises_max_end_indices, ...
                               'maxAbs', von_mises_max_abs(:));
end

function envelope = computePlaneEndEnvelope(values)
%COMPUTEPLANEENDENVELOPE Return the maximum absolute entry across planes and ends.
%   This helper is used to report extreme responses by considering both
%   section planes (e.g., top/bottom fibers) and the A/B ends of an element.
%   It provides a single scalar suitable for max-value reporting.
%
%   envelope = computePlaneEndEnvelope(values)
%       values : matrix of size (nPlanes x nEnds)
%       envelope : scalar max(abs(values(:))) or 0 if input is empty

    if isempty(values)
        envelope = 0;
        return;
    end

    envelope = max(abs(values), [], 'all');
end

function [selectedValue, planeIdx, endIdx, maxAbs] = extractShearEnvelope(values)
%EXTRACTSHEARENVELOPE Return the controlling shear entry and its metadata.
%   [value, planeIdx, endIdx, maxAbs] = extractShearEnvelope(values) finds the
%   entry with the largest absolute magnitude, returning the signed value, its
%   plane index (1 = Plane 1, 2 = Plane 2), end index (1 = End A, 2 = End B),
%   and the absolute magnitude. The function accepts any non-empty matrix and
%   treats NaNs as invalid entries.

    if isempty(values)
        selectedValue = 0;
        planeIdx = 1;
        endIdx = 1;
        maxAbs = 0;
        return;
    end

    if any(isnan(values), 'all')
        error('extractShearEnvelope:NaNInput', 'Input values contain NaNs.');
    end

    [maxAbs, linearIdx] = max(abs(values), [], 'all', 'linear');
    [planeIdx, endIdx] = ind2sub(size(values), linearIdx);
    selectedValue = values(planeIdx, endIdx);
end
