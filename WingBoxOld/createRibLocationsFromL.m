function [yrib, total_ribs] = createRibLocationsFromL(L, wing, D)
% createRibLocationsFromL_auto estimates rib locations from a local rib
% spacing distribution L, automatically determining the "ideal" number
% of ribs from the integral of 1/L over each wing half.
%
% Inputs:
%   L     - local rib spacing vector (same size as wing.stripy)
%           (NaN for points in fuselage region |y| <= D/2)
%   wing  - structure with field 'stripy' containing spanwise sample points
%   D     - fuselage width
%
% Outputs:
%   yrib       - vector of estimated rib locations (spanwise positions)
%   total_ribs - total number of ribs used (sum of left + right)
%
% Strategy:
%   1) On each wing half, define S(y) = ∫(1./L(y)) dy from the tip to the
%      fuselage boundary. The final value of S(y) suggests how many ribs
%      that half "wants." We use round() to convert to an integer.
%   2) We then place that many ribs by splitting [0, S_end] into N equal
%      increments in S-space and inverting to find y-locations.
%   3) We combine left and right results, sort them, and output.

    % Extract the spanwise grid
    y = wing.stripy(:);  % force column

    % Identify which points are on left/right wings (outside fuselage)
    left_mask  = (y < -D/2);
    right_mask = (y >  D/2);

    %-------------------------------------------------------------
    % Handle left wing
    %-------------------------------------------------------------
    y_left = y(left_mask);
    L_left = L(left_mask);

    % Remove NaN or zero/negative L
    valid_left = ~isnan(L_left) & (L_left > 0);
    y_left = y_left(valid_left);
    L_left = L_left(valid_left);

    % Sort from tip (most negative) to fuselage edge
    [y_left, sortIdx] = sort(y_left);
    L_left = L_left(sortIdx);

    if isempty(y_left)
        % No valid left wing points
        rib_left = [];
        N_left = 0;
    else
        % Compute cumulative integral of 1/L
        % This effectively counts "rib units" from the left tip to the fuselage
        S_left = cumtrapz(y_left, 1 ./ L_left);

        % The final value of S_left indicates total "rib units" on left side
        S_left_end = S_left(end);

        % "Ideal" integer # of ribs for left wing (rounding to nearest integer)
        N_left = round(S_left_end);

        if N_left < 1
            % If S_left_end < 0.5, rounding might give 0
            % Force at least 1 rib if you want at least a single bracket
            N_left = 1;
        end

        % Define target S-values from 0 to S_left_end with N_left+1 points
        % (So we get N_left intervals)
        S_target_left = linspace(0, S_left_end, N_left+1);

        % Invert S_left to get y-locations
        rib_left = interp1(S_left, y_left, S_target_left, 'linear');
    end

    %-------------------------------------------------------------
    % Handle right wing
    %-------------------------------------------------------------
    y_right = y(right_mask);
    L_right = L(right_mask);

    valid_right = ~isnan(L_right) & (L_right > 0);
    y_right = y_right(valid_right);
    L_right = L_right(valid_right);

    [y_right, sortIdx] = sort(y_right);
    L_right = L_right(sortIdx);

    if isempty(y_right)
        % No valid right wing points
        rib_right = [];
        N_right = 0;
    else
        S_right = cumtrapz(y_right, 1 ./ L_right);
        S_right_end = S_right(end);

        N_right = round(S_right_end);
        if N_right < 1
            N_right = 1;
        end

        S_target_right = linspace(0, S_right_end, N_right+1);
        rib_right = interp1(S_right, y_right, S_target_right, 'linear');
    end

    %-------------------------------------------------------------
    % Combine
    %-------------------------------------------------------------
    yrib_all = [rib_left(:); rib_right(:)];

    % Remove duplicates if you want (e.g. if the fuselage edge is repeated)
    yrib = unique(yrib_all);

    % Sort from left tip to right tip
    yrib = sort(yrib);

    % The total number of "segments" is (N_left + N_right)
    % The total # of distinct rib STATIONS is (N_left + 1) + (N_right + 1)
    % minus any duplicates. But let's define "total_ribs" as the sum:
    total_ribs = N_left + N_right;

end
