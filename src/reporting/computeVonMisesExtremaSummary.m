function summary = computeVonMisesExtremaSummary(stresses)
%COMPUTEVONMISESEXTREMASUMMARY Return max von Mises stress metadata.

    plane_labels = {'Plane 1', 'Plane 2'};
    end_labels = {'End A', 'End B'};

    defaultEntry = struct('value', 0, 'elementID', NaN, ...
        'planeIdx', 1, 'planeLabel', plane_labels{1}, ...
        'endIdx', 1, 'endLabel', end_labels{1});

    summary = defaultEntry;

    if ~isfield(stresses, 'von_mises_details') || isempty(stresses.von_mises_details)
        return;
    end

    details = stresses.von_mises_details;
    if isempty(details.maxAbs)
        return;
    end

    [~, elem_idx] = max(details.maxAbs);
    if isempty(elem_idx) || elem_idx < 1 || elem_idx > numel(details.maxAbs)
        return;
    end

    plane_idx = clampIndex(details.maxPlaneIdx(elem_idx), numel(plane_labels));
    end_idx = clampIndex(details.maxEndIdx(elem_idx), numel(end_labels));

    element_ids = stresses.element_ids;
    if isempty(element_ids)
        return;
    end
    element_id = element_ids(elem_idx);

    plane_end_values = squeeze(details.planeEndValues(elem_idx, :, :));
    vm_value = plane_end_values(plane_idx, end_idx);

    summary = struct( ...
        'value', vm_value, ...
        'elementID', element_id, ...
        'planeIdx', plane_idx, ...
        'planeLabel', plane_labels{plane_idx}, ...
        'endIdx', end_idx, ...
        'endLabel', end_labels{end_idx});
end

function idx = clampIndex(idx, upper)
%CLAMPINDEX Keep index values within expected reporting limits.

    if isnan(idx) || idx < 1
        idx = 1;
    elseif idx > upper
        idx = upper;
    end
end
