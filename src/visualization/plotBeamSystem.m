function plotBeamSystem(nodes, elements, supports, forces)
    % PLOT THE ORIGINAL BEAM STRUCTURE
    % This function creates a visual representation of the beam system before
    % any analysis is performed. It's like drawing the engineering blueprint
    % that shows the beam layout, support locations, and applied forces.
    
    fig = figure('Name', 'Beam System Configuration');
    hold on;  % Keep adding elements to the same plot
    
    % DRAW THE BEAM ELEMENTS (the main structural members)
    % Each element is a straight line segment connecting two nodes
    for i = 1:size(elements, 1)
        % Get the node numbers that this element connects
        node1_id = elements(i,2);
        node2_id = elements(i,3);
        
        % Find the actual coordinates of these nodes
        node1_row = find(nodes(:,1) == node1_id, 1);
        node2_row = find(nodes(:,1) == node2_id, 1);
        
        % Draw a thick black line representing the beam element
        plot([nodes(node1_row,2), nodes(node2_row,2)], ...
             [nodes(node1_row,3), nodes(node2_row,3)], 'k-', 'LineWidth', 2);
    end
    
    % DRAW THE NODES (connection points)
    % These are the key points where beam elements connect together
    nodeHandles = gobjects(size(nodes, 1), 1);
    for i = 1:size(nodes, 1)
        x = nodes(i,2);      % X coordinate of this node
        y = nodes(i,3);      % Y coordinate of this node
        nodeID = nodes(i,1); % ID number of this node
        
        % Draw a basic node as a white circle with black outline
        node_plot = plot(x, y, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 10);
        node_plot.UserData = struct('NodeID', nodeID);
        nodeHandles(i) = node_plot;
        
        % CHECK IF THIS NODE HAS A SUPPORT
        support_idx = find([supports.NodeID] == nodeID);
        if ~isempty(support_idx)
            support_type = lower(char(supports.Type(support_idx)));
            % Use different colors to show different types of supports:
            switch support_type
                case 'fixed'
                    % Fixed support: Can't move or rotate (like welded to a wall)
                    support_plot = plot(x, y, 'ko', 'MarkerEdgeColor', 'g', 'MarkerSize', 12);
                case 'pinned'
                    % Pinned support: Can't move but can rotate (like a hinge)
                    support_plot = plot(x, y, 'ko', 'MarkerEdgeColor', 'r', 'MarkerSize', 12);
                case 'roller'
                    % Roller support: Can move in one direction but not the other
                    support_plot = plot(x, y, 'ko', 'MarkerEdgeColor', 'b', 'MarkerSize', 12);
            end
            support_plot.UserData = struct('NodeID', nodeID);
            nodeHandles(i) = support_plot;
        end
    end
    
    % DRAW THE APPLIED FORCES
    % These are external loads applied to the beam
    for i = 1:size(forces, 1)
        % Find where this force is applied
        node_id = forces(i,1);
        node_row = find(nodes(:,1) == node_id, 1);
        x = nodes(node_row,2);
        y = nodes(node_row,3);
        
        % Get the force components
        Fx = forces(i,2);  % Force in X direction
        Fy = forces(i,3);  % Force in Y direction
        
        % Draw arrows to show force direction and magnitude
        scale = 0.25; % Scale factor to make arrows visible
        if Fx ~= 0
            % Draw horizontal force arrow
            quiver(x, y, sign(Fx)*scale, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.4);
        end
        if Fy ~= 0
            % Draw vertical force arrow
            quiver(x, y, 0, sign(Fy)*scale, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.4);
        end

        % Draw a small circular indicator for applied nodal moments.
        if size(forces, 2) >= 4
            Mz = forces(i,4);
            if Mz ~= 0
                radius = 0.2;
                direction = sign(Mz);
                if direction == 0
                    direction = 1;
                end
                if direction > 0
                    theta = linspace(pi/2, -3*pi/2, 60);  % Counter-clockwise arc
                else
                    theta = linspace(-pi/2, 3*pi/2, 60);  % Clockwise arc
                end
                arc_x = x + radius * cos(theta);
                arc_y = y + radius * sin(theta);
                plot(arc_x, arc_y, 'm', 'LineWidth', 1.5);

                % Add a simple arrow head to show rotation sense.
                head_x = arc_x(end);
                head_y = arc_y(end);
                arrow_angle = theta(end);
                head_scale = 0.07;
                arm1_x = head_x - head_scale * cos(arrow_angle + direction * pi/6);
                arm1_y = head_y - head_scale * sin(arrow_angle + direction * pi/6);
                arm2_x = head_x - head_scale * cos(arrow_angle - direction * pi/6);
                arm2_y = head_y - head_scale * sin(arrow_angle - direction * pi/6);
                line([head_x, arm1_x], [head_y, arm1_y], 'Color', 'm', 'LineWidth', 1.5);
                line([head_x, arm2_x], [head_y, arm2_y], 'Color', 'm', 'LineWidth', 1.5);

            end
        end
    end
    
    % SET UP THE PLOT APPEARANCE
    axis equal;  % Make sure the scale is the same in X and Y directions
    grid on;     % Add a grid to help read coordinates
    xlabel('X Position');
    ylabel('Y Position');
    title('2D Timoshenko Beam System Configuration');
    hold off;    % Finished adding elements to this plot

    dcm = datacursormode(fig);
    set(dcm, 'Enable', 'on', 'UpdateFcn', @customNodeDataTip);

    function outputTxt = customNodeDataTip(~, eventObj)
        pos = eventObj.Position;
        target = eventObj.Target;
        outputTxt = {sprintf('X %.4g', pos(1)), sprintf('Y %.4g', pos(2))};
        if any(target == nodeHandles) && isfield(target.UserData, 'NodeID')
            outputTxt = [{sprintf('NodeID %d', target.UserData.NodeID)}, outputTxt];
        end
    end
end
