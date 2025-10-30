function displayResults(displacements)
    % DISPLAY ANALYSIS RESULTS IN AN EASY-TO-UNDERSTAND FORMAT
    % This function takes the calculated displacements and presents them
    % in a clear way that helps users understand how much the beam deforms
    % under the applied loads.
    
    % EXTRACT DISPLACEMENT COMPONENTS
    % Separate the different types of displacements for analysis
    max_x_disp = max(abs(displacements(1:3:end)));  % Maximum horizontal movement
    max_y_disp = max(abs(displacements(2:3:end)));  % Maximum vertical movement (deflection)
    max_rot = max(abs(displacements(3:3:end)));     % Maximum rotation
    
    % CALCULATE TOTAL DISPLACEMENT MAGNITUDE
    % This combines horizontal and vertical movements into a single measure
    % of how far each point moves from its original position
    displacement_magnitude = zeros(length(displacements)/3, 1);
    for i = 1:length(displacement_magnitude)
        ux = displacements(3*i-2);  % Horizontal displacement at node i
        uy = displacements(3*i-1);  % Vertical displacement at node i
        displacement_magnitude(i) = sqrt(ux^2 + uy^2);  % Total movement magnitude
    end
    
    % DISPLAY RESULTS IN A USER-FRIENDLY FORMAT
    % Present the key results with clear explanations
    fprintf('\nTimoshenko Beam Analysis Results:\n');
    fprintf('=====================================\n');
    fprintf('Maximum Horizontal Displacement: %.6e m\n', max_x_disp);
    fprintf('Maximum Vertical Displacement:   %.6e m\n', max_y_disp);
    fprintf('Maximum Rotation:                %.6e rad (%.3f degrees)\n', max_rot, max_rot*180/pi);
    fprintf('Maximum Total Displacement:      %.6e m\n', max(displacement_magnitude));
    
    % PROVIDE CONTEXT FOR THE RESULTS
    % Help users understand whether these displacements are significant
    if max(displacement_magnitude) > 0.001  % More than 1 mm
        fprintf('Note: Large displacements detected. Check if these are acceptable for your application.\n');
    elseif max(displacement_magnitude) < 1e-6  % Less than 1 micrometer
        fprintf('Note: Very small displacements. The structure is very stiff relative to the applied loads.\n');
    else
        fprintf('Note: Moderate displacements. Structure appears to be responding normally to loads.\n');
    end
end
