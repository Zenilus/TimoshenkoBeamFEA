function k_local = calculateTimoshenkoStiffness(E, G, A, As, I, L)
    % CALCULATE TIMOSHENKO BEAM ELEMENT STIFFNESS MATRIX
    % This function creates the stiffness matrix for a single beam element
    % using Timoshenko beam theory.
    %
    % The stiffness matrix relates forces to displacements: F = K × u
    % where F is force vector, K is stiffness matrix, u is displacement vector
    
    % CALCULATE SHEAR PARAMETER
    % This parameter φ (phi) accounts for the additional flexibility
    % due to shear deformation. When φ = 0, we get Euler-Bernoulli theory.
    phi = 12 * E * I / (G * As * L^2);
    
    % AXIAL STIFFNESS TERMS
    % These terms relate axial forces to axial displacements
    % (stretching and compression along the beam length)
    k11 = E * A / L;  % Axial stiffness at node 1
    k44 = k11;        % Axial stiffness at node 2
    k14 = -k11;       % Coupling between nodes (opposite sign)
    k41 = k14;        % Symmetric coupling
    
    % BENDING AND SHEAR STIFFNESS TERMS
    % These terms are modified from Euler-Bernoulli theory to include
    % the effects of shear deformation (the φ parameter)
    
    % Transverse (perpendicular) stiffness - reduced by shear effects
    k22 = 12 * E * I / (L^3 * (1 + phi));
    k55 = k22;        % Same at both nodes
    k25 = -k22;       % Coupling between nodes
    k52 = k25;        % Symmetric coupling
    
    % Force-rotation coupling terms
    k23 = 6 * E * I / (L^2 * (1 + phi));
    k26 = k23;        % Force at node 1 coupled to rotation at node 2
    k32 = k23;        % Rotation at node 1 coupled to force at node 2
    k62 = k23;        % Similar coupling at node 2
    
    k56 = -6 * E * I / (L^2 * (1 + phi));
    k53 = k56;        % Negative coupling terms
    k65 = k56;
    k35 = k56;
    
    % ROTATIONAL STIFFNESS TERMS
    % These relate moments to rotations
    k33 = (4 + phi) * E * I / (L * (1 + phi));  % Rotational stiffness at node 1
    k66 = k33;                                   % Same at node 2
    
    k36 = (2 - phi) * E * I / (L * (1 + phi));  % Coupling between rotations
    k63 = k36;                                   % Symmetric coupling
    
    % ASSEMBLE THE 6×6 ELEMENT STIFFNESS MATRIX
    % The matrix is organized as follows:
    % Rows/Columns 1-3: Node 1 (u1, v1, θ1)
    % Rows/Columns 4-6: Node 2 (u2, v2, θ2)
    % where u = axial displacement, v = transverse displacement, θ = rotation
    k_local = [
        k11  0    0    k14  0    0   ;  % u1 equation
        0    k22  k23  0    k25  k26 ;  % v1 equation
        0    k32  k33  0    k35  k36 ;  % θ1 equation
        k41  0    0    k44  0    0   ;  % u2 equation
        0    k52  k53  0    k55  k56 ;  % v2 equation
        0    k62  k63  0    k65  k66 ;  % θ2 equation
    ];
end
