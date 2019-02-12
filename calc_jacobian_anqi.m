
function [] = calc_jacobian
%function [] = calc_jacobian
%A modification of Jörg Gumbels code, the code calulcates K for MATS
%and saves the folowing variables to a file:
%K = Jacobian (pathlength through each gridpoint for each measurement)
%O = Observation vector (y)
%beta = True atmosphere
%x = Angle along orbit (degrees)
%y = distance across orbit (km)
%z = altitude above a spherical earth (km)

% Ole Martin Christensen 2017.01.11


%------------------------------------------------------------------------

% General definitions

% Earth radius:
  Re = 6370;
  deg2rad = pi/180;
  
%-------------------------------------------------------------------------

% Choose what signals to include in the simulation:
  add_NLC = 1;
  add_Rayleigh = 0;
  
%-------------------------------------------------------------------------
% Define instrument paraneters (field of view etc.)
% Note: If DY = 0 is chosen, 2D rather than 3D tomography is simulated,
% i.e. without across track resolution (LY = 1). 

% central tangent altitude (center of FOV):
  Zt = 90;
% central tangent altitude expressed as radius from Earth center:
  Rt = Re + Zt;
  
% height and with of field of view:
  DZ = 40; % (i.e 60-100 km)
  DY = 100; % (i.e -+125 km)
  
% resolution of field of view (expressed as vertical and horisontal
% distance at tangent point):
  dZ = 0.2;
  dY = 5;

% number of vertical and horizontal image pixels:
  LZ = round (DZ / dZ) + 1;
  LY = round (DY / dY) + 1;
  
% position of vertical and horizontal image pixels:
  Z = -DZ/2 + [0:LZ-1] * dZ;
  Y = -DY/2 + [0:LY-1] * dY;
  
% Choose measurement positions (as index of x):
  imeas = [200:3:450];
% For numerical efficiency, the measurement positions (AAO of tangent 
% points) must be synchronized with the the atmospheric grid (AAO of cloud
% grid). 


% altitude of satellite above Earth surface:
  zsat = 600;
% distance between satellite and central tangent point:
  Dsat = sqrt ((Re + zsat)^2 - Rt^2);
  
% %-------------------------------------------------------------------------
% % Prepare cloud field
% 
% % Load atmospheric cloud field:
%   cloud = load('anqi_xi_90.mat');
% % variables: AAO AXO alt cloud
% % resolution 6.7x6.7x0.1 km 
% %sub indexes to use from the cloud grid (to save memory)
%  xI = 1:1:length(cloud.y); %length(cloud.AAO); (i.e ca 20 km resolution)
%  yI = 1:1:length(cloud.x);% (i.e ca 6.7 km resolution)
%  zI = 1:1:length(cloud.z); % (i.e ca 0.2 km resolution)
% 
% % rename variables:
%   x = (cloud.y(xI)+6000)/ (pi/180 * Rt);
%   y = cloud.x(yI);
%   z = cloud.z(zI);
%   beta = permute(cloud.Vat_Pertur(yI,xI,zI),[2,1,3]); %permute size anqis code has different x-y
% 
  
%-------------------------------------------------------------------------
% Prepare cloud field

% Load atmospheric cloud field:
  cloud = load('3D_cloud.mat');
% variables: AAO AXO alt cloud
% resolution 6.7x6.7x0.1 km 

%sub indexes to use from the cloud grid (to save memory)
 xI = 1:3:790; %length(cloud.AAO); (i.e ca 20 km resolution)
 yI = 1:1:length(cloud.AXO);% (i.e ca 6.7 km resolution)
 zI = 1:2:length(cloud.alt); % (i.e ca 0.2 km resolution)

% rename variables:
  AAO = cloud.AAO(xI);
  AXO = cloud.AXO(yI);
  z_km = cloud.alt(zI);
  beta = cloud.beta(xI,yI,zI);

  
% convert the y-coordinate from an angle into a distance:
  y = y * pi/180 * Rt;
 
% convert z-coordinate to height relative to central tangent altitude:
  z = z - Zt;

% resolution of cloud field:  
  dx = x(2)-x(1);
  dy = y(2)-y(1);
  dz = z(2)-z(1);

% dimensions of cloud field:  
  Lx = length (x);  
  Ly = length (y);  
  Lz = length (z);  

% start and end values of cloud field:
  xmin = x(1);
  xmax = x(Lx);
  ymin = y(1);
  ymax = y(Ly);
  zmin = z(1);
  zmax = z(Lz);
  
  
% Choose whether NLC is to be included in simulation:
  if ~add_NLC
    beta = zeros (Lx, Ly, Lz);
  end

%-------------------------------------------------------------------------
% Add molecular Rayleigh scattering

% solar flux
% from Saskatoon's STIS/Kurucz solar spectrum:
  ESun = 3.59e13;  % 270 +/- 1 nm
%  ESun = 6.60e13;  % 300 +/- 1 nm

% Rayleigh cross section (differential cross section for dcattering angle
% 90 degrees, units cm2 str-1)
  sigRayleigh = 5.36e-27;  % 270 nm
%  sigRayleigh = 3.39e-27;  % 300 nm

% Determine latitude as a function of orbit angle AAO.
% maximum latitude based on assuned orbit inclination (from Odin):
  maxlat = 82.3; 
% conversion AAO to latitude:
  lat = asin (sin(x*deg2rad) * sin(maxlat*deg2rad)) / deg2rad;  
  
% get latitude-dependent density profile from MSIS (assuming fixed date)  
% load MSIS look-up table 
% (latMsis, monthMsis, zMsis [km], NMsis [cm-3], TMsis):
  load msisatm
% choose fixed date July 1:
  NMsis = squeeze (NMsis (:,6,:));
% interpolate to atmospheric altitude grid:
  NMsis = interp1 (zMsis, NMsis, Zt+z);
% interpolate to measurement latitudes:
  NMsis = interp1 (latMsis, NMsis', lat)';    % size: Lz Lx
  
% Rayleigh scattering coefficient (as function of AAO and z, but
% independent of AXO):
  betaRayleigh  = 1e2 * NMsis * sigRayleigh;  % [m-1 str-1]

  
% add Rayleigh scattering to beta field:
  if add_Rayleigh
    for iy = 1:Ly
      beta(:,iy,:) = squeeze (beta(:,iy,:)) + betaRayleigh';
    end
  end


%-------------------------------------------------------------------------
% Calculate geometry of LOS integrations

  disp ('Prepare geometry.')

% simulation resolution along central LOS: (ppath step)
  dX = 0.1;  %1; %0.01;  % [km]
  
% limit for LOS integration:
  Xmax = sqrt ((Rt+zmax)^2 - (Rt-DZ/2)^2);
  Xmin = -Xmax;
  
% number of points along central LOS:
  LX = round (2*Xmax / dX) + 1;

% position of points along central LOS:
  X = Xmin + [0:LX-1] * dX; 
  
% for each image pixel, calculate points along LOS in atmospheric 
% coordinates xP, yP, zP. Calculate also corresponding indeces in 
% atmospheric coordinates ix, iy, iz.:  
  xP = atan (X/Rt) / deg2rad;   % size(xP) = 1 LX
  ix = round (xP/dx);
  
  yP = Y' * (1 + X/Dsat);       % size(yP) = LY LX
  iy = round ((yP - ymin)/dy + 1);
  
  for iY = 1:LY
    for iZ = 1:LZ
      zP = sqrt ( ((1+Z(iZ)/Rt)./(1-X/Dsat*Z(iZ)/Rt)).^2  .* ...
                  (Rt^2 + X.^2) + ...
                  (1+X/Dsat).^2 * Y(iY)^2 ) - Rt;
      iz = round ((zP - zmin)/dz + 1);
%     total linear index in 3D atmospheric matrix (x-index of current
%     tangent point to be added later):
      iP = ix + (iy(iY,:)-1)*Lx + (iz-1)*Lx*Ly;
      
      
%     Go along LOS and check what atmospheric cells contribute to LOS 
%     integral.
%     First, remove parts of LOS that are above the atmospheric field:
      i = find (zP > zmax);
      iP(i) = [];
      
      clear iAtm wAtm
      cnt = 1;
      % iAtm is the linear index of all the block a LOS passes through
      iAtm(cnt) = iP(1);
      % wAtm is the number of times the path is within a block (i.e. roughly pathlength through grid)
      wAtm(cnt) = 1; 
      for i = 2:length(iP)
        if iP(i) == iP(i-1)
          wAtm(cnt) = wAtm(cnt)+1;
        else
          cnt = cnt+1;
          iAtm(cnt) = iP(i);
          wAtm(cnt) = 1;
        end
      end
      
      pix(iY,iZ).iAtm = iAtm;
      pix(iY,iZ).wAtm = wAtm;

    end
  end

%-------------------------------------------------------------------------
% Perform LOS integrations

  disp ('Perform LOS integrations.')


% number of measurement positions (tangent points)
  Lm = length(imeas);

% matrix for resulting LOS radiances:
  L = zeros(Lm, LY, LZ);
  O = zeros(Lm*LY*LZ,1);
  %assign memory to the linear indecies that will create K later
  Intersection = zeros(Lm*LY*LZ*400,1);
  Measurements = zeros(Lm*LY*LZ*400,1);
  Pathlength = zeros(Lm*LY*LZ*400,1);

  
  n=1;
  m=1;
  
  for im = 1:Lm
    disp ([im Lm])
    for iY = 1:LY
      for iZ = 1:LZ
%         disp ([iY iZ])
        %Radiance of pixel iY,iZ of meansurement im = sum(scattering from
        %idexes of the atmosphere given by pix(iY,iZ).iAtm + some constant
        %from measurements
        L(im,iY,iZ) = sum (beta(pix(iY,iZ).iAtm+imeas(im)) .* pix(iY,iZ).wAtm);
       
        
        %calulate weighting funtion
%         K(n,pix(iY,iZ).iAtm+imeas(im)) = pix(iY,iZ).wAtm;
        
        %Stuff for new MART
        O(n,1) = L(im,iY,iZ); %measurement with linear index
        nelements = size(pix(iY,iZ).iAtm,2)-1; %number of non zero elements for meausurement n
        Intersection(m:m+nelements) = pix(iY,iZ).iAtm+imeas(im); %linear index of which atmospheric gridpoints measurement n crosses (column of jacobian)
        Measurements(m:m+nelements) = n;
        Pathlength(m:m+nelements) = pix(iY,iZ).wAtm; %pathlength through each gridpoint for measurement n
 
        m = m+nelements+1; %change linear index
        n=n+1; %change measurement number
      end
    end  
  end
  
  %truncate the linear vector to only contain values assigned. 
  Intersection = Intersection(1:m-1);
  Measurements = Measurements(1:m-1);
  Pathlength = Pathlength(1:m-1);
  
  
  clearvars -EXCEPT Measurements Measurements_coarse Intersection Intersection_coarse Pathlength Pathlength_coarse O L beta x y z xI yI zI
  
  K = sparse([Measurements],[Intersection],[Pathlength],length(O),length(beta(:)));
  
  clear Measurements Measurements_coarse Intersection Intersection_coarse Pathlength Pathlength_coarse 

save('MATS_fm_jacobian_anqi','K','L','O','beta','x','y','z','xI','yI','zI','-v7.3') %save results



