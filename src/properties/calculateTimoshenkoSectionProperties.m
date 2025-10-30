function [I, As, maxDistance, ky] = calculateTimoshenkoSectionProperties(sectionType, A, width, height, diameter)
    % CALCULATE CROSS-SECTION PROPERTIES FOR TIMOSHENKO BEAM THEORY
    % This function determines the geometric properties that affect how
    % the beam responds to bending and shear forces. Different cross-section
    % shapes have different properties and shear correction factors.
    %
    % Returns:
    % I = Second moment of area (resistance to bending)
    % As = Effective shear area (resistance to shear deformation)
    % maxDistance = Maximum distance from neutral axis (for stress calculation)
    % ky = Shear correction factor (accounts for non-uniform shear distribution)
    
    switch lower(char(sectionType))
        case 'rectangle'
            % RECTANGULAR CROSS-SECTION
            % Most common in construction (like wooden beams, concrete beams)
            I = (width * height^3) / 12;  % Moment of inertia formula for rectangle
            maxDistance = height / 2;     % Maximum stress occurs at top/bottom
            
            % Shear correction factor: accounts for the fact that shear stress
            % is not uniform across the cross-section (parabolic distribution)
            ky = 5/6;  % Theoretical value for rectangular sections
            As = ky * A;  % Effective shear area
            
        case 'square'
            % SQUARE CROSS-SECTION
            % Special case of rectangle where width = height
            I = (width^4) / 12;        % Same formula as rectangle
            maxDistance = width / 2;   % Distance to corner
            ky = 5/6;                  % Same correction factor as rectangle
            As = ky * A;
            
        case 'circle'
            % CIRCULAR CROSS-SECTION
            % Common for pipes, rods, and some structural members
            I = (pi * diameter^4) / 64;  % Moment of inertia for circle
            maxDistance = diameter / 2;  % Maximum stress at outer fiber
            
            % Circular sections have a different shear stress distribution
            % leading to a different correction factor
            ky = 9/10;  % Higher value than rectangular sections
            As = ky * A;
            
        otherwise
            % DEFAULT/UNKNOWN SECTION TYPE
            % When the section type is not recognized, make conservative estimates
            I = (A^2) / 12;           % Rough approximation based on area
            maxDistance = sqrt(A) / 2; % Assume square-like proportions
            ky = 5/6;                 % Use rectangular assumption (conservative)
            As = ky * A;
            fprintf('Warning: Unknown section type. Using rectangular assumptions.\n');
    end
end
