function [elemental_stresses, bending_moments, shear_forces, shear_stresses, element_positions] = ...
    calculateTimoshenkoStresses(nodes, elements, displacements, E, G, A, As, I, maxDistance)
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
    bending_moments = zeros(num_elements, 1);     % Internal bending moments
    shear_forces = zeros(num_elements, 1);        % Internal shear forces
    shear_stresses = zeros(num_elements, 1);      % Shear stresses
    element_positions = zeros(num_elements, 1);   % Position for plotting
    
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
        
        % Calculate bending curvature using Timoshenko theory
        % This represents how much the beam curves under bending loads
        curvature = (d_local(6) - d_local(3))/L;
        
        % Calculate shear strain (unique to Timoshenko theory)
        % This represents the "sliding" deformation between node ends
        % γ = (dv/dx - θ) where θ is the rotation of the cross-section
        theta_avg = (d_local(3) + d_local(6))/2;
        shear_strain = (d_local(5) - d_local(2))/L - theta_avg;
        
        % CALCULATE INTERNAL FORCES FROM STRAINS
        % These are the forces that the material develops to resist deformation
        
        % Axial force (tension or compression along the beam length)
        axial_force = E * A * axial_strain;
        
        % Bending moment (causes the beam to curve)
        bending_moments(i) = E * I * curvature;
        
        % Shear force (causes sliding between layers of the beam)
        shear_forces(i) = G * As * shear_strain;
        
        % CALCULATE STRESSES FROM FORCES
        % Stresses are forces per unit area - these determine if the material will fail
        
        % Axial stress (normal stress along the beam length)
        axial_stress = axial_force / A;
        
        % Bending stress (varies across the cross-section, maximum at outer fibers)
        bending_stress = bending_moments(i) * maxDistance / I;
        
        % Shear stress (varies across the cross-section, maximum at neutral axis)
        shear_stresses(i) = shear_forces(i) / As;
        
        % COMBINED STRESS MEASURE
        % Use von Mises equivalent stress as a measure of overall stress state
        % This combines all stress components into a single value for comparison
        % with material strength limits
        elemental_stresses(i) = sqrt(axial_stress^2 + bending_stress^2 + 3 * shear_stresses(i)^2);
    end
end
