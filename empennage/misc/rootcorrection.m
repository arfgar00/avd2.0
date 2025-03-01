% Your data (x is the second column and y is the first column)
x = [
    0.00
    0.50
    1.00
    1.50
    2.00
    2.50
    3.00
    3.50
    4.00
    4.50
    5.00
    5.50
    6.00
    6.50
    7.00
    7.50
    8.00
    8.50
    9.00
    9.50
    10.00
    10.50
    11.00
    11.50
    12.00
    12.50
    13.00
    13.50
    14.00
    14.50
    15.00
    15.50
    16.00
    16.50
    17.00
    17.50
    18.00
    18.50
    19.00
    19.50
    20.00
    20.50
    21.00
    21.50
    22.00
    22.50
    23.00
    23.50
    24.00
    24.50
    25.00
    25.50
    26.00
    26.50
    27.00
    27.50
    28.00
    28.50
    29.00
    29.30
];

y = [
    1.00
    0.979093099
    0.958186199
    0.937279298
    0.916372397
    0.895465496
    0.874558596
    0.853651695
    0.832744794
    0.8132409
    0.801499214
    0.790487652
    0.77947609
    0.768464528
    0.757452965
    0.746441403
    0.735429841
    0.724418279
    0.713406716
    0.702395154
    0.691383592
    0.68037203
    0.669360467
    0.658348905
    0.647337343
    0.636325781
    0.625314219
    0.614302656
    0.603291094
    0.592279532
    0.58126797
    0.570256407
    0.559244845
    0.548233283
    0.537221721
    0.526210158
    0.515198596
    0.504187034
    0.493175472
    0.482163909
    0.471152347
    0.460140785
    0.449129223
    0.438117661
    0.427106098
    0.416094536
    0.405082974
    0.394071412
    0.383059849
    0.372048287
    0.361036725
    0.350025163
    0.3390136
    0.328002038
    0.316990476
    0.305978914
    0.294967351
    0.283955789
    0.272944227
    0.00
];

% Normalize the x-values such that the maximum value is 1
x_normalized = x / max(x);

% Interpolate the data
% Create a finer x grid for interpolation (for example, 100 points between 0 and 1)
x_fine = linspace(min(x_normalized), max(x_normalized), 100);

% Interpolate y-values at the new normalized x positions
y_interpolated = interp1(x_normalized, y, x_fine, 'pchip');  % pchip for smooth interpolation

% To generate 'n' points, for example n = 10:
n = 10;  % You can change this value
x_generated = linspace(min(x_normalized), max(x_normalized), n);
y_generated = interp1(x_normalized, y, x_generated, 'pchip');

% Plot the generated points
figure;
plot(x_generated, y_generated, 'rx', 'MarkerSize', 10); % Generated points
title('Generated Interpolated Points');
xlabel('Normalized x');
ylabel('y');
grid on;
