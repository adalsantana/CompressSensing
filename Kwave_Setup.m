%k-wave main file for compress sensing integration with ultrasound

%input_args = [1 2];

% ----------- create the copmutational grid --------------------- 
Nx = 64;        %number of grid points in the x direction
Ny = 64;        %number of grid points in the y direction      
Nz = 64;        %number of grid points in the z direction
dx = 0.1e-3;    %grid point spacing in the x direction [m]
dy = 0.1e-3;    %grid point spacing in the y direction [m]
dz = 0.1e-3;    %grid point spacing in the z direction [m]
cx = dx / 2;
cy = dy / 2;
cz = dz / 2;

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
number_scan_lines = 96; 
% ----------- computational grid created -----------------------


%  ---- define a series of Cartesian points to collect data ---- 
% sensor mask used to define where our pressure field will be recorded
x = (-22:2:22) * dx; 
y = 22*dy*ones(size(x));
z = (-22:2:22) * dz; 
sensor.mask = [x;y;z]; 

% ------------ binary sensor mask defined ----------------------

% ------- define the properties of the homogeneous propagation medium ----
% ------- parameters should match our experimental setup -----------------
medium.sound_speed = 1500;  % [m/s]                 
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;   
medium.density = 1000;     % [kg/m^3] density for water being used for now

%parameter for nonlinearity. Why should it be set is covered in the
%ultrasound transducer kwave documentation but am still unsure what to do
%with it
%medium.BonA = 1;

% -------- define medium properties for non homogenous medium
%medium.sound_speed = 1500 * ones(Nx, Ny);   % [m/s]
%medium.sound_speed(1:Nx/2, :) = 1800;       % [m/s]
%medium.density = 1000 * ones(Nx, Ny);
%medium.density(:, Ny/4:Ny) = 1200; 

% ---------- defining time array -------------------
kgrid.makeTime(medium.sound_speed);

%----------- defining properties of input signal -------------------
source_strength = 1e6;     % [MPa]
tone_burst_freq = 0.5e6;   % [Hz]
tone_burst_cycles = 5; 

%---------- create the input signal using toneBurst ----------------
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source strength divided by the 
% impedance (the source is assigned to the particle velocity)
input_signal = input_signal(1:Nx); %Nx matches length of multiplication below
input_signal = (source_strength ./(medium.sound_speed * medium.density)) .* input_signal; 

% ---------- input signal defined and created -------------------------

%------------ define physical properties of transducer -------------------
transducer.number_elements = 1;      % total number of transducer elements
transducer.element_width = 1;       % width of each element (in grid points)
transducer.element_length = 12;     % length of each element (in grid points)
transducer.element_spacing = 0;     % spacing (kerf width) between the elements (in grid points)
% inf in this context means the transudcer is not curved (e.g. flat,
% linear)
transducer.radius = inf;            % radius of curvature of transducer

%calculate the width of the transudcer in grid points 
transducer_width = transducer.number_elements * transducer.element_width + (transducer.number_elements-1) * transducer.element_spacing; 

%use this to position the transducer in the middle of the computational
%grid
transducer.position = round([1,Ny/2-transducer_width/2,Nz/2-transducer.element_length/2]);

%properties used to derive beamforming delays
transducer.sound_speed = 1540;              % sound speed [m/s]
transducer.focus_distance = 20e-3;          %focus distance [m]
transducer.elevation_focus_distance = 19e-3;%focus distance in the elevation plane [m]
transducer.steering_angle = 0;              %steering angle 

%apodization
transducer.transmit_apodization = 'Rectangular'; 
transducer.receive_apodization = 'Rectangular'; 

%define the transudcer elements that are currently active; 
transducer.active_elements = zeros(transducer.number_elements, 1); 
transducer.active_elements = 1;

%append the input signal used to derive the transducer
transducer.input_signal = input_signal; 

%create the transducer using the defined settings 
transducer = kWaveTransducer(kgrid, transducer);

%create ball to be simulated and imaged
ball_radius = 5; %in grid points 
ball = makeBall(Nx, Ny, Nz , cx, cy, cz, ball_radius);
input_args={'PlotLayout', true, 'PlotPML', false, 'DataCast', 'gpuArray-single', 'CartInterp', 'nearest'};
scan_lines = size(number_scan_lines);
medium_position = 0; 

for scan_line_index = 1:number_scan_lines
    
    %load current section of the medium 
    medium.sound_speed = sound_speed_map(:, medium_position:medium_position + Ny - 1, :);
    medium.density = density_map(:, medium_position:medium_position + Ny - 1, :);
    
    %run the simulation
    sensor_data = kSpaceFirstOrder3D(kgrid, medium, transducer, transducer, input+args{:});
    
    %extract the scan line from the sensor data
    scan_lines(scan_line_index, :) = transducer.scan_line(sensor_data);
    
    medium_position = medium_position + transducer.element_width; 
    
    t0 = length(input_signal)*kgrid.dt/2;
    r = c0 * ((1:length(kgrid.t_array))*kgrid.dt/2-t0);
    
    %create time gain compensation function  based on attenuation value,
    % transmit freqeuncy, and round trip distance
    tgc_alpha = 0.4;
    tgc = exp(2*tgc_alpha*tone_burse_freq*1e-6*r*100);
    
    %apply the time gain compensation to each of the scan lines
    scan_lines = bsxfun(@times, tgc, scan_lines);
    
    %filter the scan lines using both the transmit frequency and the second
    %harmonic
    scan_lines_fund = gaussianFilter(scan_lines, 1/kgrid.dt, tone_burst_freq, 100, true);
    scan_lines_harm = gaussianFilter(scan_lines, 1/kgrid.dt, 2*tone_burst_freq, 30, true);
    
    %envelope detection 
    scan_lines_fund = envelopeDetection(scan_lines_fund);
    scan_lines_harm = envelopeDetection(scan_lines_harm);

    %normalized log compression
    compression_ratio = 3;
    scan_lines_fund = logCompression(scan_lines_fund, compress_ratio, true);
    scan_lines_harm = logCompression(scan_lines_harm, compress_ratio, true);
    
    %upsample the image using linear interpolation 
    scale_factor = 2; 
    scan_lines_fund = interp2(1:kgrid.Nt, (1:number_scan_lines).', scan_lines_fund, 1:kgrid.Nt, (1:1/scale_factor:number_scan_lines).');
    scan_lines_harm = interp2(1:kgrid.Nt, (1:number_scan_lines).', scan_lines_harm, 1:kgrid.Nt, (1:1/scale_factor:number_scan_lines).');
end

%run the simulation
% input_args={'PlotLayout', true, 'PlotPML', false, 'DataCast', 'gpuArray-single', 'CartInterp', 'nearest'};
% sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, sensor, input_args{:});