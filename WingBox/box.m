clear
clc

%spar section design

%bending moment parameters

D = 200;
t1 = 20;
t2 = 1.6;
b = 30;
c = 100;
Ixx = 1/6 * b * t1^3 + 1/2 * b * t1 * (D - t1)^2 + (1/12) * t2 * (D - 2*t1)^3;


M_root = 23286.8;
M_root_sigmax = M_root * 1000 * c / Ixx;
M_1m = 11668.9;
M_1m_sigmax = M_1m * 1000 * c / Ixx;
M_2m = 4332.6;
M_2m_sigmax = M_2m * 1000 * c / Ixx;
M_3m = 782.2;
M_3m_sigmax = M_3m * 1000 * c / Ixx;

% material properties

materialProperties = struct(...
    'Al2024', struct(...
        'YoungsModulus', 70, ...  % E (GPa)
        'PoissonRatio', 0.3, ...  % v
        'ShearModulus', 27, ...   % G (GPa)
        'BulkModulus', 58, ...    % K (GPa)
        'Density', 2780), ...     % Density (kg/m^3)
    'Al7075', struct(...
        'YoungsModulus', 73, ...  % E (GPa)
        'PoissonRatio', 0.36, ... % v
        'ShearModulus', 27, ...   % G (GPa)
        'BulkModulus', 87, ...    % K (GPa)
        'Density', 2810));       % Density (kg/m^3)

tensile_yield = 324;
shear_yield = 187.1;
compressive_yield = 350

%torsion parameters

torsion_parameters = struct(...
    'An', 11000, ...         % Cell 1 area (mm^2)
    'Ar', 16500, ...         % Cell 2 area (mm^2)
    't_w', 1.6, ...          % Spar thickness (mm)
    'tn', 1.2, ...           % Skin thickness (mm)
    'tr', 1.2, ...           % Tr (mm)
    'd', 120, ...            % Spar distance to SC (mm)
    'P', -17274.76262, ...   % Load (N)
    'Sn', 670, ...           % Circumference 1 (Sn) (mm)
    'Sr', 310, ...           % Circumference 1 (Sr) (mm)
    'h_w', 200, ...          % Spar height (mm)
    'G', 25000);             % Shear modulus (N/mm^2)

torque = torsion_parameters.P / torsion_parameters.d






matrix1 = [2*torsion_parameters.An 2*torsion_parameters.Ar 0; torsion_parameters.h_w -torsion_parameters.h_w -torsion_parameters.h_w; 1/2/torsion_parameters.An/torsion_parameters.G * torsion_parameters.Sn / torsion_parameters.tn -1/2/torsion_parameters.Ar/torsion_parameters.G * torsion_parameters.Sr / torsion_parameters.tr torsion_parameters.h_w/torsion_parameters.t_w*1/2/torsion_parameters.An/torsion_parameters.G +  torsion_parameters.h_w/torsion_parameters.t_w*1/2/torsion_parameters.Ar/torsion_parameters.G ]
vector1 = [torque; torsion_parameters.P; 0];

solution = matrix1 \ vector1;

tau_n = solution(1) / torsion_parameters.tn;
tau_r = solution(2) / torsion_parameters.tr;
tau_w = solution(3) / torsion_parameters.t_w;


%add checks here
P = (solution(1) - solution(2) - solution(3))*torsion_parameters.h_w;
T = 2*torsion_parameters.An*solution(1) + 2*torsion_parameters.Ar*solution(2);
theta_xn =

