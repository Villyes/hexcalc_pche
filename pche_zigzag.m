HTXzigzag = pche();
HTXzigzag.name = "HTXzigzag";

% charge=1
% 1: Side A (hot)
% 2: Steel 
% 3: Side B (cold)

HTXzigzag.lambda_st=17.5; %SA 240 Grade 316L Stainless Steel, 16.3 @100°C, 21.5 @500°C, W/mK

HTXzigzag.flow_type = "custom"; %inlet and outlet area linear transition cross-flow to countercurrent, otherwise countercurrent
HTXzigzag.process_medium1 = "CO2";
HTXzigzag.process_medium3 = "CO2";

HTXzigzag.d1=2e-3;          % diameter of channel area (zigzag)
HTXzigzag.w1=0;             % distance between quarter circles for "other_zigzag", see paper "streched semicircle" for a sketch; for semicircular: w=0  
HTXzigzag.pitch1=0.01;      % zigzag channel pitch
HTXzigzag.alpha1=110;       % zigzag channel pitch, make sure to take obtuse angle
HTXzigzag.nChannel1=19;     % number of channels per plate
HTXzigzag.nPlates1=30;      % number of plates
HTXzigzag.ePlate1=3e-3;     % total thickness of plate

HTXzigzag.d3=2e-3;          % diameter of channel area
HTXzigzag.w3=0;             % distance between quarter circles for "other_zigzag"; for semicircular: w=0    
HTXzigzag.pitch3=0.01;      % zigzag channel pitch
HTXzigzag.alpha3=110;       % zigzag channel pitch, make sure to take obtuse angle
HTXzigzag.nChannel3=25;     % number of channels per plate
HTXzigzag.nPlates3=30;      % number of plates
HTXzigzag.ePlate3=3e-3;     % total thickness of plate

HTXzigzag.l=0.35;           % effective hex length (m)                       [1,1]

HTXzigzag.l_in=0.2;         % percentage meaning: length of inlet and outlet area where linear transition between cross-flow and countercurrent