function segments = segmentDirArray(data, threshold)
    segments = [];
    segmentStart = 1;
    segmentSum = 0;
    
    for i = 2:length(data)
        diffValue = abs(data(i) - data(i-1));
        segmentSum = segmentSum + diffValue;
        
        if segmentSum > threshold
            segment = [segmentStart, i-1];
            segments = [segments; segment];
            segmentStart = i;
            segmentSum = 0;
        end
    end
    
    % Check if there is a segment at the end
    if segmentStart <= length(data)
        segment = [segmentStart, length(data)];
        segments = [segments; segment];
    end
end
