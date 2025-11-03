function [nodes, elements, supports, forces, E, A, density, sectionType, width, height, diameter, G, nu] = readInputData(inputWorkbook)
    % READ INPUT DATA FROM EXCEL FILE
    % This function reads all the structural information needed for the analysis
    % from an Excel file. This acts as the "blueprint" of our beam structure.

    if nargin < 1 || isempty(inputWorkbook)
        inputWorkbook = 'input_data.xlsx';
    end

    if ~isfile(inputWorkbook)
        error('Input workbook not found at %s', inputWorkbook);
    end

    try
        % LOAD DATA FROM DIFFERENT SHEETS IN THE EXCEL FILE
        % Each sheet contains different types of information:
        % - Nodes: The key points along the beam
        % - Elements: The beam segments that connect nodes together
        % - Supports: Where the beam is supported/constrained using fixed, pinned, or roller supports
        % - Forces: External loads applied to the beam
        % - Properties: Material characteristics

        nodes = readtable(inputWorkbook, 'Sheet', 'Nodes');
        elements = readtable(inputWorkbook, 'Sheet', 'Elements');
        supports = readtable(inputWorkbook, 'Sheet', 'Supports');
    forces = readtable(inputWorkbook, 'Sheet', 'Forces');
        properties = readtable(inputWorkbook, 'Sheet', 'Properties');
        
        % STANDARDIZE COLUMN NAMES
        % Make sure the data columns have the expected names, regardless of
        % how the user named them in Excel. This prevents errors later.
        
        if ~isequal(nodes.Properties.VariableNames, {'NodeID', 'X', 'Y'})
            % Standard format: ID number, X position, Y position
            nodes.Properties.VariableNames = {'NodeID', 'X', 'Y'};
        end
        
        if ~isequal(elements.Properties.VariableNames, {'ElementID', 'Node1', 'Node2'})
            % Standard format: ID number, first node, second node
            elements.Properties.VariableNames = {'ElementID', 'Node1', 'Node2'};
        end
        
        if ~isequal(supports.Properties.VariableNames, {'NodeID', 'Type'})
            % Standard format: which node, what type of support
            supports.Properties.VariableNames = {'NodeID', 'Type'};
        end
        
        % Standardize forces column names, allowing an optional nodal moment column.
        forceVarNames = forces.Properties.VariableNames;
        standardizedForceNames = cell(size(forceVarNames));
        for idx = 1:numel(forceVarNames)
            nameLower = lower(forceVarNames{idx});
            switch nameLower
                case {'nodeid', 'node', 'id'}
                    standardizedForceNames{idx} = 'NodeID';
                case {'fx', 'forcex'}
                    standardizedForceNames{idx} = 'Fx';
                case {'fy', 'forcey'}
                    standardizedForceNames{idx} = 'Fy';
                case {'mz', 'moment', 'momentz', 'torsion', 'torque'}
                    standardizedForceNames{idx} = 'Mz';
                otherwise
                    standardizedForceNames{idx} = forceVarNames{idx};
            end
        end
        forces.Properties.VariableNames = standardizedForceNames;

        requiredForceColumns = {'NodeID', 'Fx', 'Fy', 'Mz'};
        if ~all(ismember(requiredForceColumns, forces.Properties.VariableNames))
            error('Forces sheet must contain NodeID, Fx, Fy, and Mz columns.');
        end

        forces = forces(:, requiredForceColumns);

        for idx = 1:numel(requiredForceColumns)
            colName = requiredForceColumns{idx};
            if ~isnumeric(forces.(colName))
                error('Forces sheet column %s must be numeric.', colName);
            end
            if any(~isfinite(forces.(colName)))
                error('Forces sheet column %s contains missing or invalid values.', colName);
            end
        end
        
        % CHECK FOR ADVANCED MATERIAL PROPERTIES
        expected_props = {'YoungsModulus', 'CrossSectionalArea', 'Density', 'SectionType', 'Width', 'Height', 'Diameter', 'ShearModulus', 'PoissonRatio'};
        if ~all(ismember(expected_props, properties.Properties.VariableNames))
            % If advanced properties aren't provided, use basic format
            if size(properties, 2) < 9
                properties.Properties.VariableNames = {'YoungsModulus', 'CrossSectionalArea', 'Density'};
            end
        end
        
        % EXTRACT BASIC MATERIAL PROPERTIES
        % These are the fundamental properties every structural analysis needs:
        E = properties.YoungsModulus(1);     % How stiff the material is (like steel vs rubber)
        A = properties.CrossSectionalArea(1); % How much area the beam cross-section has
        density = properties.Density(1);      % How heavy the material is per unit volume
        
        % SET DEFAULT VALUES FOR CROSS-SECTION GEOMETRY
        % If the user didn't specify the exact beam shape, we'll make reasonable assumptions
        sectionType = 'Rectangle';  % Assume rectangular cross-section
        width = sqrt(A);           % Make it square with the given area
        height = sqrt(A);
        diameter = sqrt(4*A/pi);   % Equivalent circle diameter
        
        % EXTRACT DETAILED GEOMETRY IF PROVIDED
        % If the user specified the exact beam shape, use those values instead
        if ismember('SectionType', properties.Properties.VariableNames)
            sectionType = string(properties.SectionType(1));
        end
        
        if ismember('Width', properties.Properties.VariableNames)
            width = properties.Width(1);
        end
        
        if ismember('Height', properties.Properties.VariableNames)
            height = properties.Height(1);
        end
        
        if ismember('Diameter', properties.Properties.VariableNames)
            diameter = properties.Diameter(1);
        end
        
        % DETERMINE SHEAR PROPERTIES
        % Timoshenko theory needs to know how the material responds to shear forces
        if ismember('ShearModulus', properties.Properties.VariableNames)
            G = properties.ShearModulus(1);  % Direct shear modulus value
        else
            % If not provided, calculate it from Young's modulus and Poisson's ratio
            % This is a standard relationship in materials science
            nu = 0.3;  % Typical value for steel
            G = E / (2 * (1 + nu));
            fprintf('Warning: Shear modulus not provided. Calculated G = %.3e Pa assuming nu = %.2f\n', G, nu);
        end
        
        if ismember('PoissonRatio', properties.Properties.VariableNames)
            nu = properties.PoissonRatio(1);  % How much the material "squeezes" when stretched
            if ~exist('G', 'var')
                G = E / (2 * (1 + nu));  % Recalculate shear modulus with actual Poisson's ratio
            end
        else
            if ~exist('nu', 'var')
                nu = 0.3; % Default assumption for steel
                fprintf('Warning: Poisson ratio not provided. Using default nu = %.2f\n', nu);
            end
        end
        
        % CONVERT DATA TO NUMERICAL ARRAYS
        % MATLAB works better with numerical arrays than with table formats
        nodes = table2array(nodes);
        elements = table2array(elements);
        % Keep supports as a special structure because it has both numbers and text
        supports = struct('NodeID', table2array(supports(:,1)), ...
                         'Type', string(supports.Type));
        forces = table2array(forces);
        if size(forces, 2) >= 4
            num_moments = sum(forces(:,4) ~= 0);
        else
            num_moments = 0;
        end
        
        % DISPLAY SUMMARY OF LOADED DATA
        % Let the user know what we successfully loaded
        fprintf('Beam Data loaded successfully!\n');
        %fprintf('Nodes: %d\n', size(nodes, 1));
        %fprintf('Elements: %d\n', size(elements, 1));
        %fprintf('Supports: %d\n', length(supports.NodeID));
        %fprintf('Forces: %d\n', size(forces, 1));
        %fprintf('Moments: %d\n', num_moments);
        
    catch ME
        fprintf('Error reading input file: %s\n', ME.message);
        fprintf('Attempted workbook path: %s\n', inputWorkbook);
        fprintf('\nExpected Excel file structure for Timoshenko analysis:\n');
        fprintf('1. Nodes sheet:\n   NodeID, X, Y\n');
        fprintf('2. Elements sheet:\n   ElementID, Node1, Node2\n');
        fprintf('3. Supports sheet:\n   NodeID, Type\n');
        fprintf('4. Forces sheet:\n   NodeID, Fx, Fy\n');
        fprintf('5. Properties sheet:\n   YoungsModulus, CrossSectionalArea, Density, SectionType, Width, Height, Diameter, ShearModulus, PoissonRatio\n');
        fprintf('   SectionType should be "Rectangle", "Square", or "Circle"\n');
        fprintf('   ShearModulus and PoissonRatio are required for accurate Timoshenko analysis\n');
        rethrow(ME);
    end
end
