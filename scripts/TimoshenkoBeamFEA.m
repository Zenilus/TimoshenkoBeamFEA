function TimoshenkoBeamFEA(inputWorkbook)
    % TIMOSHENKO BEAM FINITE ELEMENT ANALYSIS
    % This is the main function that performs a complete structural analysis
    % of a beam using Timoshenko beam theory, which accounts for both bending
    % and shear deformation (unlike simpler Euler-Bernoulli theory).

    % Ensure the src tree is available on the MATLAB path no matter where the
    % user launches the run script from.
    currentFile = mfilename('fullpath');
    projectRoot = fileparts(fileparts(currentFile));
    addpath(genpath(fullfile(projectRoot, 'src')));

    % STEP 1: Load the problem setup from an Excel file
    % This reads all the structural information: where nodes are located,
    % how elements connect nodes, where supports are placed, what forces
    % are applied, and what material properties the beam has.
    if nargin < 1 || isempty(inputWorkbook)
        dataDir = fullfile(projectRoot, 'data');
        inputWorkbook = fullfile(dataDir, 'input_data.xlsx');
    end
    [nodes, elements, supports, forces, E, A, density, sectionType, width, height, diameter, G, nu] = readInputData(inputWorkbook);

    % STEP 2: Show the original beam structure
    % This creates a visual representation of the beam before any loads are applied,
    % showing the beam elements, support locations, and applied forces.
    plotBeamSystem(nodes, elements, supports, forces);

    % STEP 3: Solve the structural problem using Timoshenko beam theory
    % This is the heart of the analysis - it calculates how much each point
    % on the beam will move (displacements) and what internal forces and
    % stresses develop when the loads are applied.
    [displacements, stresses] = performTimoshenkoFEA(nodes, elements, supports, forces, E, A, density, sectionType, width, height, diameter, G, nu);

    % STEP 4: Show the numerical results
    % This displays the maximum displacements helping us
    % understand how much the beam deflects under load.
    displayResults(displacements);

    % STEP 5: Show the deformed beam shape
    % This creates a visual representation of how the beam looks after
    % the loads are applied, with colors showing displacement magnitudes.
    plotDeformedShape(nodes, elements, displacements);

    % STEP 6: Display stress analysis results
    % These are the internal forces and stresses that develop in the beam.
    stressEntries = prepareStressReportEntries(stresses);
    reportOrder = {
        'maxBendingMoment', ...
        'maxBendingStressPlane1', ...
        'maxBendingStressPlane2', ...
        'maxBendingStressEnvelope', ...
        'maxShearForce', ...
        'maxShearStress', ...
        'maxVonMisesStress'
    };

    for idx = 1:numel(reportOrder)
        entry = stressEntries.(reportOrder{idx});
        fprintf('%s\n', entry.text);
    end
end
