function color = getColorFromValue(value, max_value)
    % HELPER FUNCTION: MAP DISPLACEMENT VALUE TO COLOR
    % This function converts a displacement magnitude into a color for visualization
    % Normalize the displacement value to a 0-1 range for color mapping
    if max_value == 0
        normalized = 0;  % Handle case where there's no displacement
    else
        normalized = value / max_value;  % Scale to 0-1 range
    end
    
    % Use MATLAB's jet colormap (blue to red spectrum)
    % Blue = low displacement, Red = high displacement
    cmap = jet(256);  % Create 256-color jet colormap
    
    % Map the normalized value to a colormap index
    index = max(1, min(256, round(normalized * 255) + 1));
    
    % Return the corresponding RGB color
    color = cmap(index, :);
end
