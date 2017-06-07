clear all; %#ok<CLALL>
close all;
clc;

%-------- Definitions --------

Ts = 293.15;    % Surface temperature [K]
Ps = 0.84;      % Surface pressure [atm] at 4793 ft (Bozeman)
gamma = -7;     % Adiabatic lapse rate [K/km]

bg = 0.6e-6;    % Backscatter at ground
bl = 3;         % Boundary layer height [km]
blw = 0.5;      % Boundary layer falloff width [km]
sa = 50;        % Aerosol lidar ratio
sm = 8.3776;    % Molecular lidar ratio (8/3)*pi

c = 3.00e8;          % Speed of light [m/s]
n0 = 2.687e25;       % Loschmidt constant [1/m^3] at 0 C and 1 atm
kb = 1.38065e-23;    % Boltzmann constant [J/K]
h = 6.626e-34;       % Planck constant [J s]

tpulse = 0.67e-6;       % Pulse duration [s]
rangebin = c*tpulse/2;  % Range bin size [m]
r = 300:rangebin:6000;  % Heights [m]
rkm = r./1000;          % Heights [km]

m = 5.314e-26;                  % Mass of O2 molecule [kg]
wl = 769.233;                   % Wavelength of laser [nm]
wn0 = 1.299996e6;               % Wavenumber of laser [m^-1]
n = 301;                        % Number of steps in wavenumber range
wn1 = 1.299982e6;               % Beginning of wavenumber range
wn2 = 1.30001e6;                % End wavenumber range
wn = wn1:((wn2-wn1)/(n-1)):wn2; % Wavenumbers [m^-1]

%% --------- Profiles ---------
T = Ts+gamma.*rkm;                      % Temperature [K]
P = Ps*(Ts./T).^(-5.2199);              % Pressure [atm]
Pkpa = 101.325*P;                       % Pressure [kPa]
nO2 = 0.2096*n0*T.*P/(273.15*1);        % O2 number density [m^-3]

betam = 374280*Pkpa./(T*769.233^4);     % Molecular Backscatter
extm = sm*betam;                        % Molecular Extinction

% Aerosol Backscatter
k = size(r);
k = k(2);
for j=1:1:k
    range = rkm(j);
    betaa(j) = bg/100;
    if range<bl+blw/2
        betaa(j) = bg*(1-(bl-blw/2)/10)*(1-(range-(bl-blw/2))/blw);
    end
    if range<bl-blw/2
        betaa(j) = bg*(1-range/10);
    end
end
exta = sa*betaa;                      % Aerosol Extinction

% Doppler Broadened Backscatter Lineshape
D = (((m*c^2)./(8*pi*wn0^2*kb.*T)).*exp(-((m*c^2)./(8*pi*wn0^2*kb.*T)).*(wn'-wn0).^2)).^0.5;

%% -------- Absorption Cross Section Calculation ---------
fwidth = 20;                        % Range to scan [GHz]
fstepnum = 301;                     % Number of steps
fstart = c*wn0-fwidth*1e9/2;        % Beginning of freq range [s^-1]
fstep = fwidth*1e9/(fstepnum-1);    % Step width [s^-1]
fstop = c*wn0+fwidth*1e9/2;         % End of freq range [s^-1]

nu0 = 1/(769.233e-9*100);           % Laser wavenumber [cm^-1]
S0 = 1.08e-25;                      % (HITRAN) line strength [cm/molecule]
t0 = 293.15;                        % Reference temperature [K]
p0 = 1;                             % Reference pressure [atm]
ccm = 3e10;                         % Speed of light [cm/s]
E = 1248.204;                       % Energy above the g.s. of the absorption line's lower level [cm^-1]
gl0 = 0.0362;                       % Halfwidth [cm^-1] (gl0*c = 1.085 GHz)


j=1;
for freq = fstart:fstep:fstop
    f(j) = (freq-wn0*c)/1e9;        % Distance from laser wavenumber [GHz]
    lambda(j) = c/freq;             % Wavelengths within range [m]
    j=j+1;
end

for j=1:1:k
    t = T(j);
    pa = P(j);
    a = S0*(t0/t);
    b = (1-exp(-h*ccm*nu0/(kb*t)))/(1-exp(-h*ccm*nu0/(kb*t0)));
    d = exp(h*ccm*E*(1/t0-1/t)/kb);
    st = a*b*d;                                         % Nehrir Eq. 2.15
    
    gammaL = gl0*(pa/p0)*(t0/t)^0.71;                   % Nehrir Eq. 2.16
    gammaD = (nu0/ccm)*(2*kb*t*log(2)/m)^0.5*100;       % Nehrir Eq. 2.21
    i=1;
    
    for freq = fstart:fstep:fstop
        nu = freq/c/100;
        
        % Calculate the convolution integral
        x = (nu-nu0)*0.8325546/gammaD;
        y = gammaL*0.8325546/gammaD;
        ttt = -3:0.0001:3;
        ft = exp(-ttt.^2)./(y.*y+(x-ttt).^2);
        convint = trapz(ttt,ft);
        sigmaV(j,i) = st*0.12448*gammaL/(gammaD^2)*convint;
        
        i=i+1;
    end
end

cs = (sigmaV(:,151))';                          % cross section of absorbers at online wavelength
alpha = cs.*nO2/(100^2);                        % absorption coefficient for online wavelength
alphab = sigmaV.*repmat(nO2',1,301)/(100^2);    % repeat array nO2 to create matrix same size as sigmaV
                                                % alphab is absorption coefficient. 
                                                % Contains info for any height and wavelength in our ranges.

%% -------- Calculate Returns ---------
% First find optical depth from absorption coefficient (alphab).
numrow = size(alphab,1);                    % Number of rows in alpha matrix
optdepth = [];

for i=1:numrow
    if i==1
        optdepth(i,:) = alphab(1,:);
    else 
        alphai = alphab(1:i,:);     % Selects range of rows from i to last row
        optdepth(i,:) = trapz(alphai);
    end
end

% Next find molecular absorption: Tm = exp[-int(alpha)dr from 0 to r]
Tm = exp(-optdepth);       % Fraction of light reaching a range r

% Find atmospheric transmission
extm = sm*betam;           % Molecular extinction (Casey ~Eq. 45)
exta = sa*betaa;           % Aerosol extinction (Casey ~Eq. 45)
ncol = size(betam,2);      % Number of columns in vector betam (and betaa)
for i=1:ncol
    exti = extm(:,1:i)+exta(:,1:i); % Total extinction
    int_ext(:,i) = trapz(exti);     % Integrate total extinction from 0 to r (keep range info)
end

Ta = exp(-int_ext);             % Atmospheric Transmission

% This was for a non-delta laser lineshape.
% FWHM of diode laser is about 1.5 nm. Translates to 0.0076 GHz. 
% FWHM = 2.355*stddev
%   stddev = 0.00322718;            % standard deviation [GHz]
%   hi = exp((-f.^2)/(2*stddev^2)); % Laser spectral line shape

%   int_Tu = hi'.*Tm';              % Integrad of Tu integral. Transposed so 'trapz' integrates over frequency
%   Tu = trapz(int_Tu);             % Transmittance going up thru atmosphere


% For a delta laser lineshape hi = delta(nu-nu0),
% Tu is integral(hi Tm) over frequency, but delta function picks out nu0:
Tu = Tm(:,151)';                % Fraction of nu0 light reaching range r

N0 = 3.8697e13;                 % Number of photons initially emitted if pulse is 10 uJ
Nu = N0*Ta.*Tu;                 % Number of photons that reach the scatterer

betat = betaa+betam;            % Total backscatter

gi = (betaa./betat)*hi(151)+(betam./betat).*D;  % Return lineshape
Ns = Nu.*betat.*gi;             % Number of photons that backscatter

A = 0.25;                       % Area of telescope entrance pupil [m^2]
Nt = Nu.*betat.*gi.*Ta.*Tm'.*(A./r.^2)*rangebin;    %Number of photons seen by the telescope

overlap = 1;                    % Overlap function
epsilon = 0.25;                 % Reciever transmission*detector efficiency*max etalon tramsmission

% Etalon Transmission calculations
wl_start = wl-0.01972;          % Beginning of wavelength range [nm] (equivalent to 10 GHz away from online)
wl_end = wl+0.01972;            % End of wavelength range [nm]
lambda = wl_start:(wl_end-wl_start)/(n-1):wl_end;   % Wavelength range (made to have 301 points)

R = 0.92956;                    % Etalon reflectivity
FSR = 0.0592;                   % Free Spectral Range (FSR=c/(2*nL)) difference between online and offline wl [nm]
theta = (2*pi*(lambda-wl))/FSR; % Round trip phase accumulation
etalonT = (1-R)^2./(pi*(1+R^2-2*R*cos(theta)));   % Normalized etalon transmission
%plot(lambda,etalonT)           % Should be centered at 769.233 nm

Nd = Nt'.*overlap*epsilon.*etalonT;  % Number of photons seen by the detector

int_Td = gi.*Tm'.*etalonT';
Td = trapz(int_Td);             % Transmittance down thru atmosphere

Nr = Nu.*betat.*Ta.*overlap*epsilon.*(A./r.^2)*rangebin.*Td;  % Number of photons seen by the receiver
%semilogx(Nr,rkm)               % Log plot of Nr (# of photons) vs height (km)

%% The DIAL Equation
rangebinkm = rangebin/1000;     % Rangebin in km
dgi_dr = gradient(gi);          % partial derivative gi wrt r
int_numG = trapz(dgi_dr.*etalonT'.*Tm');
int_denG = trapz(gi.*etalonT'.*Tm');
G = int_numG./int_denG;         % Correction factor

int_num_alphaD = trapz(gi.*etalonT'.*alphab'.*Tm');
int_den_alphaD = trapz(gi.*etalonT'.*Tm');
alphaD = int_num_alphaD./int_den_alphaD;

int_num_alphaU = alphab(:,151).*Tm(:,151);
int_den_alphaU = Tm(:,151);
alphaU = int_num_alphaU./int_den_alphaU;


%% --------- Plots ------------
figure
surf(f,rkm,optdepth,'edgecolor','none')
title('\textbf{Optical Depth}','Interpreter','latex')
xlabel('Frequency ($f-f_0$) [GHz]','Interpreter','latex')
ylabel('Height AGL [km]','Interpreter','latex')
view(2)
colorbar

figure
surf(f,rkm,alphab,'edgecolor','none')
colormap('parula')
title('\textbf{Absorption Coefficient}','Interpreter','latex')
xlabel('Frequency ($f-f_0$) [GHz]','Interpreter','latex')
ylabel('Height AGL [km]','Interpreter','latex')
view(2)
colorbar

alpha1 = alphab(1,:);           % Absorption coefficient at 300 m (min)
alpha2 = alphab(57,:);          % Absorption coefficient at 6 km (max)
plot(f,alpha1,f,alpha2)
title('\textbf{Absorption Coefficient}','Interpreter','latex')
xlabel('Frequency ($f-f_0$) [GHz]','Interpreter','latex')
ylabel('Lineshape','Interpreter','latex')

surf(rkm,wn,D,'edgecolor','none')
title('Doppler Broadened Backscatter Lineshape')
xlabel('r [km]')
ylabel('Wavenumber [m^{-1}]')
colorbar
view(2)

figure
cs1 = D(:,1);       % First column, all rows
cs2 = D(:,57);      % Last column, all rows
plot(wn,cs1,wn,cs2) % Cross sections of initial and final lineshape
title('Doppler Broadened Backscatter Lineshape')
xlabel('Wavenumber [m^{-1}]')
ylabel('Intensity')

figure
subplot(2,1,1)
plot(T,rkm)
title('Atmospheric Temperature')
xlabel('T [K]')
ylabel('r [km]')
subplot(2,1,2)
plot(P,rkm)
title('Atmospheric Pressure')
xlabel('P [atm]')
ylabel('r [km]')

figure
subplot(2,2,1)
plot(betam,rkm)
title('Molecular Backscatter')
xlabel('\beta_m')
ylabel('r [km]')

subplot(2,2,2)
plot(betaa,rkm)
title('Aerosol Backscatter')
xlabel('\beta_a')
ylabel('r [km]')

subplot(2,2,3)
plot(extm,rkm)
title('Molecular Extinction')
xlabel('\sigma_m')
ylabel('r [km]')

subplot(2,2,4)
plot(exta,rkm)
title('Aerosol Extinction')
xlabel('\sigma_a')
ylabel('r [km]')