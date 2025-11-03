function [displacements, stresses] = performTimoshenkoFEA(nodes, elements, supports, forces, E, A, density, sectionType, width, height, diameter, G, nu)
    % PERFORM TIMOSHENKO BEAM FINITE ELEMENT ANALYSIS
    % This is the main computational function that solves the structural problem.
    %
    % The Timoshenko beam theory is more accurate than simpler beam theories because
    % it accounts for both bending deformation AND shear deformation. This is
    % especially important for thick beams or beams with low shear stiffness.
    
    % DETERMINE PROBLEM SIZE
    % Count how many nodes and elements we have to work with
    num_nodes = size(nodes, 1);
    num_elements = size(elements, 1);
    
    % DEGREES OF FREEDOM (DOF) SETUP
    % In Timoshenko beam theory, each node has 3 possible movements:
    % 1. Horizontal displacement (u)
    % 2. Vertical displacement (v)
    % 3. Rotation (θ) - how much the cross-section rotates
    num_dof_per_node = 3;
    total_dof = num_dof_per_node * num_nodes;
    
    % INITIALIZE STRUCTURAL MATRICES
    % K = stiffness matrix (relates forces to displacements)
    % F = force vector (all external forces applied to the structure)
    K = zeros(total_dof, total_dof);
    F = zeros(total_dof, 1);
    
    % CALCULATE CROSS-SECTION PROPERTIES
    % These geometric properties determine how the beam responds to loads
    [I, As, maxDistance, ky] = calculateTimoshenkoSectionProperties(sectionType, A, width, height, diameter);
    
    % DISPLAY PROBLEM INFORMATION
    % Let the user know what we're analyzing
    fprintf('\nBeam Structure Information:\n');
    fprintf('Number of nodes: %d\n', num_nodes);
    fprintf('Number of elements: %d\n', num_elements);
    fprintf('Number of supports: %d\n', length(supports.NodeID));
    fprintf('Number of forces: %d\n', size(forces, 1));
    num_moments = 0;
    if size(forces, 2) >= 4
        num_moments = sum(forces(:,4) ~= 0);
    end
    fprintf('Number of moments: %d\n', num_moments);
    
    % DISPLAY MATERIAL AND GEOMETRIC PROPERTIES
    % Show the beam characteristics that affect its behavior
    fprintf('\nBeam Properties:\n');
    fprintf('Cross-section Type: %s\n', sectionType);
    fprintf('Cross-sectional Area: %.3e m²\n', A);
    fprintf('Effective Shear Area: %.3e m²\n', As);
    fprintf('Moment of Inertia: %.3e m⁴\n', I);
    fprintf('Shear Correction Factor: %.3f\n', ky);
    fprintf('Young''s Modulus: %.3e Pa\n', E);
    fprintf('Shear Modulus: %.3e Pa\n', G);
    fprintf('Poisson''s Ratio: %.3f\n', nu);
    fprintf('Material Density: %.3e kg/m³\n', density);
    
    % DISPLAY DETAILED GEOMETRY INFORMATION
    % Show specific dimensions based on the cross-section shape
    switch lower(char(sectionType))
        case 'rectangle'
            fprintf('Rectangle dimensions: Width = %.3e m, Height = %.3e m\n', width, height);
        case 'square'
            fprintf('Square dimension: Width = %.3e m\n', width);
        case 'circle'
            fprintf('Circle dimension: Diameter = %.3e m\n', diameter);
        otherwise
            fprintf('Using default approximation for section properties\n');
    end
    
    fprintf('Maximum distance from neutral axis: %.3e m\n', maxDistance);
    fprintf('Section modulus: %.3e m³\n', I/maxDistance);
    
    % BUILD THE GLOBAL STIFFNESS MATRIX
    % This is done by calculating each element's stiffness matrix and
    % assembling them into one large matrix representing the entire structure
    for i = 1:num_elements
        % Get the two nodes that this element connects
        node1_id = elements(i,2);
        node2_id = elements(i,3);
        
        % Find the coordinates of these nodes
        node1_row = find(nodes(:,1) == node1_id, 1);
        node2_row = find(nodes(:,1) == node2_id, 1);
        
        % CALCULATE ELEMENT GEOMETRY
        % Each element has its own length and orientation in 2D space
        x1 = nodes(node1_row,2); y1 = nodes(node1_row,3);
        x2 = nodes(node2_row,2); y2 = nodes(node2_row,3);
        L = sqrt((x2-x1)^2 + (y2-y1)^2);  % Element length
        c = (x2-x1)/L;  % Cosine of angle (for coordinate transformation)
        s = (y2-y1)/L;  % Sine of angle (for coordinate transformation)
        
        % CREATE ELEMENT STIFFNESS MATRIX
        % This represents how this individual element resists deformation
        k_local = calculateTimoshenkoStiffness(E, G, A, As, I, L);
        
        % TRANSFORM FROM LOCAL TO GLOBAL COORDINATES
        % Each element has its own "local" coordinate system, but we need
        % to work in the global coordinate system for the entire structure
        T = [c s 0 0 0 0;
             -s c 0 0 0 0;
             0 0 1 0 0 0;
             0 0 0 c s 0;
             0 0 0 -s c 0;
             0 0 0 0 0 1];
        
        % Apply the coordinate transformation
        k_global = T' * k_local * T;
        
        % ASSEMBLE INTO GLOBAL STIFFNESS MATRIX
        % Add this element's contribution to the overall structural stiffness
        % Get the degrees of freedom for this element
        dof = [
            3*node1_row-2; 3*node1_row-1; 3*node1_row;    % u1, v1, θ1
            3*node2_row-2; 3*node2_row-1; 3*node2_row     % u2, v2, θ2
        ];
        
        % Add the element stiffness to the global matrix
        for ii = 1:6
            for jj = 1:6
                K(dof(ii), dof(jj)) = K(dof(ii), dof(jj)) + k_global(ii,jj);
            end
        end
    end
    
    % ADD EXTERNAL FORCES TO THE FORCE VECTOR
    % These are the loads that the user specified
    for i = 1:size(forces,1)
        % Find which node this force is applied to
        node_id = forces(i,1);
        node_row = find(nodes(:,1) == node_id, 1);
        
        % Add the force components to the global force vector
        F(3*node_row-2) = F(3*node_row-2) + forces(i,2);  % X-direction force
        F(3*node_row-1) = F(3*node_row-1) + forces(i,3);  % Y-direction force
        
        % If there's a moment (rotational force), add that too
        if size(forces,2) >= 4
            F(3*node_row) = F(3*node_row) + forces(i,4);  % Applied moment
        end
    end
    
    % ADD SELF-WEIGHT FORCES AND CONSISTENT MOMENTS
    % The beam's own weight is treated as a uniform distributed load that
    % induces both shear forces and end moments in each element.
    grav = 9.81;  % Acceleration due to gravity (m/s^2)
    tol = 1e-12;  % Numerical tolerance to avoid adding tiny values
    for i = 1:num_elements
        % Get the nodes for this element
        node1_id = elements(i,2);
        node2_id = elements(i,3);
        node1_row = find(nodes(:,1) == node1_id, 1);
        node2_row = find(nodes(:,1) == node2_id, 1);

        % Calculate element geometry and orientation
        x1 = nodes(node1_row,2); y1 = nodes(node1_row,3);
        x2 = nodes(node2_row,2); y2 = nodes(node2_row,3);
        L = sqrt((x2-x1)^2 + (y2-y1)^2);
        c = (x2-x1)/L;
        s = (y2-y1)/L;

        % Transformation matrix for this element (local -> global)
        T = [c s 0 0 0 0;
             -s c 0 0 0 0;
             0 0 1 0 0 0;
             0 0 0 c s 0;
             0 0 0 -s c 0;
             0 0 0 0 0 1];

        % Uniform load per unit length in global coordinates
        w_global = [0; -density * A * grav];  % [wx; wy]

        % Express distributed load in element local axes
        qx_local = c * w_global(1) + s * w_global(2);
        qy_local = -s * w_global(1) + c * w_global(2);

        % Build consistent nodal load vector in local coordinates
        f_local = zeros(6,1);

        if abs(qx_local) > tol
            axial_eq = qx_local * L / 2;
            f_local(1) = f_local(1) + axial_eq;
            f_local(4) = f_local(4) + axial_eq;
        end

        if abs(qy_local) > tol
            shear_eq = qy_local * L / 2;
            moment_eq = qy_local * L^2 / 12;
            f_local(2) = f_local(2) + shear_eq;
            f_local(5) = f_local(5) + shear_eq;
            f_local(3) = f_local(3) + moment_eq;
            f_local(6) = f_local(6) - moment_eq;
        end

        % Convert local nodal loads back to global coordinates
        f_global = T' * f_local;

        % Scatter into the global force vector
        dof = [
            3*node1_row-2; 3*node1_row-1; 3*node1_row;
            3*node2_row-2; 3*node2_row-1; 3*node2_row
        ];
        F(dof) = F(dof) + f_global;
    end
    
    % APPLY BOUNDARY CONDITIONS (SUPPORTS)
    % Supports prevent certain movements at specific nodes
    % We need to remove the corresponding degrees of freedom from our equations
    free_dof = 1:total_dof;  % Start with all DOFs free to move
    constrained_dof = [];    % List of DOFs that are constrained by supports
    
    % Go through each support and constrain the appropriate DOFs
    for i = 1:length(supports.NodeID)
        node_id = supports.NodeID(i);
        node_row = find(nodes(:,1) == node_id, 1);
        support_type = lower(char(supports.Type(i)));
        
        % Different support types constrain different movements:
        switch support_type
            case 'fixed'
                % Fixed support: no translation AND no rotation
                constrained_dof = [constrained_dof, 3*node_row-2:3*node_row];
            case 'pinned'
                % Pinned support: no translation but can rotate
                constrained_dof = [constrained_dof, 3*node_row-2:3*node_row-1];
            case 'roller'
                % Roller support: can't move vertically but can move horizontally
                constrained_dof = [constrained_dof, 3*node_row-1];
        end
    end
    
    % Update the list of free (unconstrained) degrees of freedom
    free_dof = setdiff(free_dof, constrained_dof);
    
    % SOLVE THE SYSTEM OF EQUATIONS
    % This is where we actually calculate the displacements
    displacements = zeros(total_dof, 1);  % Initialize displacement vector
    
    % Extract only the parts of the stiffness matrix and force vector
    % that correspond to free (unconstrained) degrees of freedom
    Kff = K(free_dof, free_dof);
    Ff = F(free_dof);
    
    % CHECK FOR NUMERICAL PROBLEMS
    % If the stiffness matrix is ill-conditioned, the solution may be unreliable
    if rcond(Kff) < 1e-12
        warning('Stiffness matrix is ill-conditioned. Check boundary conditions.');
    end
    
    % SOLVE FOR DISPLACEMENTS
    % This is the fundamental equation: K × u = F, so u = K^(-1) × F
    displacements(free_dof) = Kff \ Ff;
    
    % CALCULATE STRESSES AND INTERNAL FORCES
    % Now that we know the displacements, we can calculate the internal
    % forces and stresses that develop in each element
    [elemental_stresses, bending_moments, shear_forces, shear_stresses, element_positions, bending_stress_planes, bending_details, shear_details, von_mises_details] = ...
        calculateTimoshenkoStresses(nodes, elements, displacements, E, G, A, As, I, maxDistance, sectionType, width, height, diameter, density);
    
    % ORGANIZE RESULTS
    % Package all the stress information into a convenient structure
    stresses = struct('elemental', elemental_stresses, ...
                     'bending_moments', bending_moments, ...
                     'bending_stresses', bending_stress_planes, ...
                     'bending_details', bending_details, ...
                     'shear_forces', shear_forces, ...
                     'shear_stresses', shear_stresses, ...
                     'von_mises', von_mises_details.maxAbs, ...
                     'von_mises_details', von_mises_details, ...
                     'element_positions', element_positions, ...
                     'element_ids', elements(:,1), ...
                     'shear_details', shear_details);
end
