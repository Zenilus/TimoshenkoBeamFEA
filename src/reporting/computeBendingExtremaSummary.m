function summary = computeBendingExtremaSummary(stresses)
%COMPUTEBENDINGEXTREMASUMMARY Gather max bending responses with metadata.
%
%   summary = computeBendingExtremaSummary(stresses) returns a structure
%   containing the extreme bending moment and bending stress values along
%   with their controlling element/plane/end locations. It assumes the
%   stresses struct produced by performTimoshenkoFEA includes the
%   bending_details bookkeeping data populated by calculateTimoshenkoStresses.

    plane_labels = {'Plane 1', 'Plane 2'};
    end_labels = {'End A', 'End B'};

    if ~isfield(stresses, 'bending_details') || isempty(stresses.bending_details)
        error('computeBendingExtremaSummary:MissingDetails', ...
            'Bending details are not available in the stresses struct.');
    end

    details = stresses.bending_details;
    num_elements = numel(stresses.element_ids);

    % Default placeholders handle empty models gracefully.
    defaultMoment = struct('value', 0, 'signedValue', 0, 'elementID', NaN, ...
        'planeIdx', 1, 'planeLabel', plane_labels{1}, ...
        'endIdx', 1, 'endLabel', end_labels{1});
    defaultStress = defaultMoment;

    summary = struct('moment', defaultMoment, ...
                     'stressPlane', repmat(defaultStress, 1, 2), ...
                     'stressEnvelope', defaultMoment);

    if num_elements == 0
        return;
    end

    % ---- Bending moment envelope across planes/ends ----
    [moment_value, idx_elem] = max(details.moment.maxAbs);
    if ~isempty(moment_value) && idx_elem > 0
        plane_idx = clampIndex(details.moment.maxPlaneIdx(idx_elem), 2);
        end_idx = clampIndex(details.moment.maxEndIdx(idx_elem), 2);
        summary.moment = struct( ...
            'value', moment_value, ...
            'signedValue', details.moment.maxSigned(idx_elem), ...
            'elementID', stresses.element_ids(idx_elem), ...
            'planeIdx', plane_idx, ...
            'planeLabel', plane_labels{plane_idx}, ...
            'endIdx', end_idx, ...
            'endLabel', end_labels{end_idx});
    end

    % ---- Plane-by-plane bending stress extrema ----
    plane_values = details.stress.planeMaxAbs;
    plane_signed = details.stress.planeMaxSigned;
    plane_end_idx = details.stress.planeMaxEndIdx;
    for plane = 1:min(size(plane_values, 2), 2)
        [plane_val, idx_elem_plane] = max(plane_values(:, plane));
        if isempty(plane_val) || idx_elem_plane == 0
            continue;
        end
        end_idx = clampIndex(plane_end_idx(idx_elem_plane, plane), 2);
        summary.stressPlane(plane) = struct( ...
            'value', plane_val, ...
            'signedValue', plane_signed(idx_elem_plane, plane), ...
            'elementID', stresses.element_ids(idx_elem_plane), ...
            'planeIdx', plane, ...
            'planeLabel', plane_labels{plane}, ...
            'endIdx', end_idx, ...
            'endLabel', end_labels{end_idx});
    end

    % ---- Global bending stress envelope across planes/ends ----
    [envelope_val, idx_env] = max(details.stress.maxAbs);
    if ~isempty(envelope_val) && idx_env > 0
        plane_idx = clampIndex(details.stress.maxPlaneIdx(idx_env), 2);
        end_idx = clampIndex(details.stress.maxEndIdx(idx_env), 2);
        summary.stressEnvelope = struct( ...
            'value', envelope_val, ...
            'signedValue', details.stress.maxSigned(idx_env), ...
            'elementID', stresses.element_ids(idx_env), ...
            'planeIdx', plane_idx, ...
            'planeLabel', plane_labels{plane_idx}, ...
            'endIdx', end_idx, ...
            'endLabel', end_labels{end_idx});
    end

    % Ensure every entry has a valid location label when elements exist.
    if num_elements > 0
        for plane = 1:2
            if isnan(summary.stressPlane(plane).elementID)
                summary.stressPlane(plane).elementID = stresses.element_ids(1);
                summary.stressPlane(plane).planeIdx = plane;
                summary.stressPlane(plane).planeLabel = plane_labels{plane};
                summary.stressPlane(plane).endIdx = 1;
                summary.stressPlane(plane).endLabel = end_labels{1};
            end
        end

        if isnan(summary.stressEnvelope.elementID)
            summary.stressEnvelope.elementID = stresses.element_ids(1);
            summary.stressEnvelope.planeIdx = 1;
            summary.stressEnvelope.planeLabel = plane_labels{1};
            summary.stressEnvelope.endIdx = 1;
            summary.stressEnvelope.endLabel = end_labels{1};
        end
    end
end

function idx = clampIndex(idx, upper)
%CLAMPINDEX Ensure indices fall within the valid labeling range.

    if isnan(idx) || idx < 1
        idx = 1;
    elseif idx > upper
        idx = upper;
    end
end
