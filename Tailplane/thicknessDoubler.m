function tc = thicknessDoubler(t,y,D)
    % input an ideal continuous distribution, convert to discrete with skin doubler
    % d : length of the doubler
    
    % Example indices
    idx1 = 136;
    idx2 = 100;
    idx3 = 50;
    
    td = [t(idx1), t(idx2), t(idx3)];
    yd = abs([y(idx1), y(idx2),y(idx3)]);
    
    % Initialize tc to be the same size as y
    tc = zeros(size(y));
    
    for i = 1:length(y)
        % Use only the i-th value of y
        yabs = abs(y(i));
        
        % thickness 1
        if (yabs > D/2) && (yabs < yd(2))
            tc(i) = td(1);
        end
        
        % thickness 2
        if (yabs >= yd(2)) && (yabs < yd(3))
            tc(i) = td(2);
        end

        % thickness 3
        if (yabs >= yd(3))
            tc(i) = td(3);
        end
    end
    tc = ceil(tc * 1e3) / 1e3;
end
