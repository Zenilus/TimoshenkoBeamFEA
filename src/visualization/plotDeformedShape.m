function plotDeformedShape(nodes, elements, displacements)
    % PLOT THE DEFORMED BEAM SHAPE WITH DISPLACEMENT VISUALIZATION
    % This function creates a visual representation of how the beam looks
    % after the loads are applied. It shows both the original shape and
    % the deformed shape, with colors indicating displacement magnitudes.
    %
    % This visualization helps engineers and students understand:
    % 1. Where the beam deforms the most
    % 2. The overall deformation pattern
    % 3. How the structure responds to the applied loads
    
    % Create a new figure for the deformed shape plot
    figure('Name', 'Timoshenko Beam Deformed Shape with Displacement Gradient');
    hold on;
    
    % CALCULATE DISPLACEMENT MAGNITUDES FOR VISUALIZATION
    % We need to know how much each node moves to create the color gradient
    node_count = size(nodes, 1);
    displacement_magnitudes = zeros(node_count, 1);
    
    for i = 1:node_count
        % Get the displacement components for this node
        ux = displacements(3*i-2);  % Horizontal displacement
        uy = displacements(3*i-1);  % Vertical displacement
        % Calculate the total movement magnitude (Pythagorean theorem)
        displacement_magnitudes(i) = sqrt(ux^2 + uy^2);
    end
    
    % DETERMINE SCALING FOR VISUALIZATION
    % We need to scale the displacements to make them visible on the plot
    max_displacement = max(displacement_magnitudes);
    
    % Choose an appropriate scale factor
    scale_factor = 20;  % Default scaling
    if max_displacement > 0
        % Adjust scale factor to make displacements clearly visible
        % but not so large that they dominate the plot
        scale_factor = min(20, 0.2 / max_displacement);
    end
    
    % PLOT THE ORIGINAL STRUCTURE FOR REFERENCE
    % Show the undeformed beam in light gray so users can compare
    % the original and deformed shapes
    for i = 1:size(elements, 1)
        % Get the nodes that this element connects
        node1_id = elements(i,2);
        node2_id = elements(i,3);
        
        % Find the coordinates of these nodes
        node1_row = find(nodes(:,1) == node1_id, 1);
        node2_row = find(nodes(:,1) == node2_id, 1);
        
        % Draw the original element as a thin gray line
        plot([nodes(node1_row,2), nodes(node2_row,2)], ...
             [nodes(node1_row,3), nodes(node2_row,3)], 'color', [0.8, 0.8, 0.8], 'LineWidth', 1);
    end
    
    % SET UP COLOR MAPPING FOR DISPLACEMENT VISUALIZATION
    % Use colors to show displacement magnitudes - hot colors for large displacements
    colormap(jet);  % Red = high displacement, Blue = low displacement
    c_bar = colorbar;
    ylabel(c_bar, 'Displacement Magnitude (m)');
    caxis([0, max_displacement]);  % Set color scale limits
    
    % PLOT THE DEFORMED BEAM ELEMENTS WITH COLOR GRADIENT
    % Show how the beam looks after deformation, with colors indicating displacement
    for i = 1:size(elements, 1)
        % Get the nodes that this element connects
        node1_id = elements(i,2);
        node2_id = elements(i,3);
        
        % Find the coordinates and displacements of these nodes
        node1_row = find(nodes(:,1) == node1_id, 1);
        node2_row = find(nodes(:,1) == node2_id, 1);
        
        % Original node coordinates
        x1 = nodes(node1_row,2);
        y1 = nodes(node1_row,3);
        x2 = nodes(node2_row,2);
        y2 = nodes(node2_row,3);
        
        % Get the displacement components for both nodes
        u1x = displacements(3*node1_row-2);  % Node 1 horizontal displacement
        u1y = displacements(3*node1_row-1);  % Node 1 vertical displacement
        u2x = displacements(3*node2_row-2);  % Node 2 horizontal displacement
        u2y = displacements(3*node2_row-1);  % Node 2 vertical displacement
        
        % Calculate the new (deformed) positions
        x1_def = x1 + scale_factor * u1x;
        y1_def = y1 + scale_factor * u1y;
        x2_def = x2 + scale_factor * u2x;
        y2_def = y2 + scale_factor * u2y;
        
        % Get displacement magnitudes for color mapping
        disp_mag1 = displacement_magnitudes(node1_row);
        disp_mag2 = displacement_magnitudes(node2_row);
        
        % CREATE GRADIENT COLORING ALONG THE ELEMENT
        % Divide each element into segments with gradually changing colors
        num_points = 10;  % Number of segments per element
        x_points = linspace(x1_def, x2_def, num_points);
        y_points = linspace(y1_def, y2_def, num_points);
        c_points = linspace(disp_mag1, disp_mag2, num_points); % Color values
        
        % Plot each segment with its own color
        for j = 1:num_points-1
            line([x_points(j), x_points(j+1)], ...
                 [y_points(j), y_points(j+1)], ...
                 'Color', getColorFromValue(c_points(j), max_displacement), ...
                 'LineWidth', 3);
        end
        
        % Mark the deformed node positions
        scatter(x1_def, y1_def, 25, disp_mag1, 'filled', 'MarkerEdgeColor', 'k');
        scatter(x2_def, y2_def, 25, disp_mag2, 'filled', 'MarkerEdgeColor', 'k');
    end
    
    % SET UP THE PLOT APPEARANCE
    grid on;                              % Add grid for easier reading
    axis equal;                           % Equal scaling in both directions
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title('Timoshenko Beam Deformed Shape with Displacement Gradient');
    
    % ADD INFORMATIVE LEGEND
    % Help users understand what they're seeing
    legend_str = sprintf('Deformation Scale: %.0fx', scale_factor);
    legend(legend_str, 'Location', 'best');
    
    % ADD HELPFUL TEXT FOR INTERPRETATION
    if max_displacement > 0
        text_str = sprintf('Max Displacement: %.2e m\nGray lines show original position', max_displacement);
        text(0.02, 0.98, text_str, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
             'BackgroundColor', 'white', 'EdgeColor', 'black');
    end
    
    hold off;
end
