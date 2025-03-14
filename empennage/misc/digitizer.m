clear
pngname = 'shear.png';
scr_s0 = struct(...
        'ts_t_ratio', {}, ...
        'x_data', {}, ...     % Store x-values (As/bt)
        'y_data', {});        % Store y-values (σ_cr/σ₀)

%input data points from excel, and use spline interpolation.
img = imread(pngname);
% Display the image


figure(1);
clf;
imshow(img);
hold on; % Keep the image displayed while adding points
display("click origin, max x, max y in sequence")
[x,y] = ginput(3)

x_pixel_min = x(1); % Pixel x-coordinate of the origin
x_pixel_max = x(2); % Pixel x-coordinate of the maximum x-value
y_pixel_min = y(3); % Pixel y-coordinate of the maximum y-value
y_pixel_max = y(1); % Pixel y-coordinate of the origin

x_data_min = 0;   % Data x-coordinate of the origin
x_data_max = input("input x max");  % Data x-coordinate of the maximum x-value
y_data_min = 0;   % Data y-coordinate of the origin
y_data_max = input("input y max");   % Data y-coordinate of the maximum y-value

% Create mapping functions
map_x = @(x_pixel) x_data_min + (x_data_max - x_data_min) * ...
                   (x_pixel - x_pixel_min) / (x_pixel_max - x_pixel_min);

map_y = @(y_pixel) y_data_min + (y_data_max - y_data_min) * ...
                   (y_pixel_max - y_pixel) / (y_pixel_max - y_pixel_min);

% Prompt the user to select points
while true
    disp('Click on the plot to digitize points. Press ENTER when done.');
    [x_pixel, y_pixel] = ginput;
    
    % Map pixel coordinates to data coordinates
    x_data = map_x(x_pixel);
    y_data = map_y(y_pixel);
    
    % Display the digitized points
    plot(x_pixel, y_pixel, 'ro'); % Plot on the image for verification
    disp('Digitized points (x_data, y_data):');
    disp([x_data, y_data]);
    
    
    n = length(scr_s0);
    ts_t= input("enter ts/t ratio");
    scr_s0(n+1).ts_t_ratio = ts_t;
    scr_s0(n+1).x_data = x_data;
    scr_s0(n+1).y_data = y_data;
    figure(1)
    save("scr_s0.mat","scr_s0")
end