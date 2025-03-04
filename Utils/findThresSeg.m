function segments = findThresSeg(array, threshold)
    % Initialize variables
    segments = [];
    startIdx = -1;
    endIdx = -1;
    
    % Iterate over the array
    for i = 1:length(array)
        % Check if current element is above the threshold
        if array(i) > threshold
            % If this is the start of a new segment
            if startIdx == -1
                startIdx = i;
                endIdx = i;
            else
                % If this is part of an existing segment
                endIdx = i;
            end
        else
            % If the current element is below the threshold
            if startIdx ~= -1
                % Add the segment to the result
                segments = [segments; [startIdx endIdx]];
                startIdx = -1;
                endIdx = -1;
            end
        end
    end
    
    % Check if the last segment is open-ended
    if startIdx ~= -1
        segments = [segments; [startIdx endIdx]];
    end
end
