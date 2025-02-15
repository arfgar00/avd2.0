classdef material
    properties
        E %Youngs modulus
        sigma_y %yield stress
        nu     %Poissons ratio
        G      %shear modulus
        K      %BulkModulus
        rho    %Density
    end

    methods
        function obj = untitled2(inputArg1,inputArg2)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

% materialProperties = struct(...
%     'Al2024', struct(...
%         'YoungsModulus', 70, ...  % E (GPa)
%         'PoissonRatio', 0.3, ...  % v
%         'ShearModulus', 27, ...   % G (GPa)
%         'BulkModulus', 58, ...    % K (GPa)
%         'Density', 2780), ...     % Density (kg/m^3)
%     'Al7075', struct(...
%         'YoungsModulus', 73, ...  % E (GPa)
%         'PoissonRatio', 0.36, ... % v
%         'ShearModulus', 27, ...   % G (GPa)
%         'BulkModulus', 87, ...    % K (GPa)
%         'Density', 2810));       % Density (kg/m^3)