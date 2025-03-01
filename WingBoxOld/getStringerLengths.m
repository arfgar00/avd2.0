function [L1, L2, L3] = getStringerLengths(wbN, stripy, D)
    % getStringerLengths calculates the spanwise (y-direction) length
    % of three stringer increment regions.
    %
    % Inputs:
    %   wbN    - Array of stringer counts along the wing (e.g., wb.N).
    %   stripy - Array of y-coordinates (from wing tip to wing tip).
    %   D      - Diameter.
    %
    % Only the positive stripy values larger than D/2 are considered.
    
    % Find indices corresponding to positive y-values larger than D/2.
    idx = find(stripy > D/2);
    
    % Extract corresponding sections of the stringer count and y coordinates.
    Npos = wbN(idx);
    stripyPos = stripy(idx);
    
    % Find where the stringer count changes (i.e. the boundaries between regions).
    changeIdx = find(diff(Npos) ~= 0);
    
    % Check that we have exactly two changes (thus three regions).
    if numel(changeIdx) < 2
        error('Expected three regions (two changes) but found fewer.');
    end
    
    % Compute stringer length as the difference between the last and first y in each block.
    % Region 1: from the first positive value to the change point.
    L1 = stripyPos(changeIdx(1)) - stripyPos(1);
    
    % Region 2: from the first index after the first change to the second change.
    L2 = stripyPos(changeIdx(2)) - stripyPos(changeIdx(1) + 1);
    
    % Region 3: from the first index after the second change to the end.
    L3 = stripyPos(end) - stripyPos(changeIdx(2) + 1);
end