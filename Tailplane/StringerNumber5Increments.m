function Nc = StringerNumber5Increments(N,y,D)
    % Convert continuous stringer distribution to 5 discrete increments
    % N: continuous number of stringers
    % y: spatial coordinates
    % D: characteristic dimension (defines lower bound at D/2)
    
    % Example indices marking increment boundaries
    idx = [136, 120, 90, 60, 30];
    
    % Get boundary values and corresponding N values
    yd = abs(y(idx));  % Get y-values at indices
    Nd = N(idx);       % Get N-values at indices
    
    % Sort boundaries in descending order (highest y first)
    [yd_sorted, sortIdx] = sort(yd, 'descend');
    Nd_sorted = Nd(sortIdx);
    
    % Initialize output array
    Nc = zeros(size(y));
    
    for i = 1:length(y)
        yabs = abs(y(i));
        
        % Check if below minimum operational value
        if yabs < D/2
            Nc(i) = 0;
            continue;
        end
        
        % Find which increment the point belongs to
        if yabs >= yd_sorted(1)
            Nc(i) = Nd_sorted(1);
        elseif yabs >= yd_sorted(2)
            Nc(i) = Nd_sorted(2);
        elseif yabs >= yd_sorted(3)
            Nc(i) = Nd_sorted(3);
        elseif yabs >= yd_sorted(4)
            Nc(i) = Nd_sorted(4);
        else % Between D/2 and lowest boundary
            Nc(i) = Nd_sorted(5);
        end
    end
    
    Nc = round(Nc);  % Final rounding to nearest integer
end