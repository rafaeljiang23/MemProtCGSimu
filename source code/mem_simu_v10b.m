%%%%%%%%
%%% modified from DOI: 10.1126/sciadv.abo5295
%%% each molecule has x, y, a coordinate
%%% force fields:
%%%     1. LJ potential for VDW interactions
%%%     2. A energy landscape EE (theta, phi, rs) for site-site interactions
%%% distance normalized to protein size
%%% general geometry function (discretized)
%%% mobility and Kb set to unity

%%% output to text file
%%% version v10b: address errors raised by crowding
%%%%%%%%
%% simulation parameters
%%% system dimension
Lx = Lxinit;
Ly = Lyinit;

%%% total number of molecules
N = Ninit;    % default: 64*4, number of protomer

%%% defining positions and angle of the particles
x = zeros(N, 1);
y = zeros(N, 1);
a = zeros(N, 1);

xv = zeros(N, 1);
yv = zeros(N, 1);
av = zeros(N, 1);
%%% defining useful simulation parameters
rdelta = 0.1;     % skin depth
adelta = 0.01;    % rotation allowance

% h = 0.001;      % delta time, defined in run_MemCGSimu_Assembly.m
% Nsteps =  1500000;       % total number of time-steps, defined in run_MemCGSimu_Assembly.m

gamma = 5;       % friction
gamma2 = gamma;     % angular friction
T = 1;        % temperature
kb = 1;        % boltzman
beta = 1/(kb*T);

%%% bond between protomers
% rmin = 0.1;      % minimal distance between sites, defined in run_MemCGSimu_Assembly.m
% rbond = 0.9;     % reaction distance, defined in run_MemCGSimu_Assembly.m

bond = zeros(N, N);

%%% bond making/breaking
kbond = 100000000;       % prefacdtor rate which controls the time scale of bond
rv = 2*rbond + 1 + rdelta;         % neighborlist distance
rc = zeros(N, N);       % potential bond counter, 1 means can form a bond
nlist = zeros(N, 1);        % tracking number of neighbor
vlist = zeros(N, N);        % tracking neighbor particle indice
nattni = zeros(N, 1);        % tracking neighbor particle (not bound) indice
nattbi = zeros(N, 1);        % tracking neighbor particle (bound) indice

%%% force calculation
rcut = 2^(1/6);      % cutoff for the interaction in the WCA(LJ) potential, default: 2^(1/6)
ecut = 1;       % epislon in the WCA(LJ) potential

Klj = zeros(N, 3);        % force field LJ potential
Keb = zeros(N, 3);        % force field harmonic bond potential

%%% geometry funciton 
Grr = geometry_trpv3_v2(1000);
%%% setup random number generator
% seed = 1;  % simulation seed, defined in run_MemCGSimu_Assembly.m
rng1 = rng(seed);

%% initialization positions
x = xinit;
y = yinit;
a = ainit;

for j = 1:N
    % periodic boundaries
    if x(j) > Lx
        x(j) = x(j) - Lx*floor(x(j)/Lx);
    end
    if x(j) < 0
        x(j) = x(j) + Lx*(1+floor(-x(j)/Lx));
    end
    if y(j) > Ly
        y(j) = y(j) - Ly*floor(y(j)/Ly);
    end
    if y(j) < 0
        y(j) = y(j) + Ly*(1+floor(-y(j)/Ly));
    end
    if a(j) > 2*pi
        a(j) = a(j) - 2*pi*floor(a(j)/(2*pi));
    end
    if a(j) < 0
        a(j) = a(j) + 2*pi*(1+floor(-a(j)/(2*pi)));
    end
end

bond = bondinit;
EE = EEinit;
%% output file name
icr = 100;      % time-step increment for writing output files

currenttime = datetime('now','TimeZone','local','Format','d-MMM-y_HH-mm-ss_Z');
filename = string(DataHash(currenttime));
dirname = string(currenttime) + "_ID-" + filename;

mkdir(dirname);
cd(dirname);

filename_info = filename + '_info.txt';
filename_xya = filename + '_xya.txt';
filename_hb = filename + '_harmonic.txt';

%% write initial files
fileID_info = fopen(filename_info,'w');

fprintf(fileID_info,'%10s\n', string(currenttime));
fprintf(fileID_info,'%1s\n', '');
fprintf(fileID_info,'%1s\n','general parameters');
fprintf(fileID_info,'%1s: %10d\n','dimension X', Lx);
fprintf(fileID_info,'%1s: %10d\n','dimension Y', Ly);
fprintf(fileID_info,'%1s: %10d\n','particle number', N);
fprintf(fileID_info,'%1s: %10.2f\n','skin depth', rdelta);
fprintf(fileID_info,'%1s: %10.2f\n','angle allowance', adelta);
fprintf(fileID_info,'%1s: %10.6f\n','time-step size', h);
fprintf(fileID_info,'%1s: %10d\n','total time-step', Nsteps);
fprintf(fileID_info,'%1s: %10.2f\n','friction', gamma);
fprintf(fileID_info,'%1s: %10.2f\n','friction (angular)', gamma2);
fprintf(fileID_info,'%1s: %10.2f\n','temperature', T);
fprintf(fileID_info,'%1s: %10.2f\n','boltzman', kb);
fprintf(fileID_info,'%1s\n', '');
fprintf(fileID_info,'%1s\n','protomer bond (site to site)');
fprintf(fileID_info,'%1s: %1s\n','energy landscape', "saved as other file");
fprintf(fileID_info,'%1s: %10.2f\n','minimal distance', rmin);
fprintf(fileID_info,'%1s: %10.2f\n','reaction distance', rbond);
fprintf(fileID_info,'%1s: %10.2f\n','reaction rate prefactor', kbond);
fprintf(fileID_info,'%1s\n', '');
fprintf(fileID_info,'%1s\n','VDW force calculation (WCA/LJ potential)');
fprintf(fileID_info,'%1s: %10.2f\n','reaction distance', rcut);
fprintf(fileID_info,'%1s: %10.2f\n','energy well depth', ecut);
fprintf(fileID_info,'%1s: %10d\n','time-step increment for writing output files', icr);
fprintf(fileID_info,'%1s: %10.2f\n','seed', seed);

fileID_xya = fopen(filename_xya,'w');
A = [ones(N,1)*0, (1:N)', x, y, a];
fprintf(fileID_xya,'%10s %10s %10s %10s %10s\n','step', 'ID','x','y','a');
fprintf(fileID_xya,'%10d %10d %10.6f %10.6f %10.6f\n', A');

fileID_hb = fopen(filename_hb,'w');
fprintf(fileID_hb,'%10s %10s %10s\n','step', 'particle 1', 'particle 2');

for k = 1:N-1
    for j = k+1:N
        if bond(k, j) == 1
            fprintf(fileID_hb,'%10d %10d %10d\n', [0 k j]);
        end
    end
end
%% run the dynamics
write = 0;       % flag for writing the output files
flag = 1;       % flag for making neighborlist, 1 in the first step
fprintf(fileID_info,'%1s\n', '');
fprintf(fileID_info,'%1s\n','log');
fprintf(fileID_info,'%10s %1s\n', string(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')), 'start');

count_asso = 0;
count_disso = 0;
count_abnormal_LJ = 0;

sz = size(EE);
dd1 = 2*pi/sz(1);   %theta
dd2 = 2*pi/sz(2);   %phi
dd3 = 2*rbond/sz(3);   %rs
for f = 1:Nsteps
    %%% step 1: making neighbor list
    f*h
    if flag == 1
        xv = x;
        yv = y;
        av = a;

        nlist = nlist .* 0;
        vlist = vlist .* 0;
        for k = 1:N-1
            for j = k+1:N
                 if bond(k, j) == 1
                    vlist(k, nlist(k) + 1) = j;       % neighbor particle indice
                    nlist(k) = nlist(k) + 1;
                    vlist(j, nlist(j) + 1) = k;       % neighbor particle indice
                    nlist(j) = nlist(j) + 1;
                    rc(k, j) = 1;          % capable of forming a bond;
                    rc(j, k) = 1;
                else
                    rc(k, j) = 0;
                    rc(j, k) = 0;
                    xr = x(k) - x(j);
                    if xr > Lx*0.5
                        xr = xr - Lx;
                    elseif xr < (-Lx*0.5)
                        xr = xr + Lx;
                    end

                    yr = y(k) - y(j);
                    if yr > Ly*0.5
                        yr = yr - Ly;
                    elseif yr < (-Ly *0.5)
                        yr = yr + Ly;
                    end

                    r = sqrt(xr^2 + yr^2);

                    if r <= rv        % cutoff for neighbor
                        vlist(k, nlist(k) + 1) = j;       % neighbor particle indice
                        nlist(k) = nlist(k) + 1;
                        vlist(j, nlist(j) + 1) = k;       % neighbor particle indice
                        nlist(j) = nlist(j) + 1;
                        if r <= 2*rbond + 1
                            rc(k, j) = 1;          % capable of forming a bond;
                            rc(j, k) = 1;
                        end

                    end
                end
            end
        end
    end
    flag = 0;       % zero the flag for the neigborlist hereW
    % When particles are found outside of the cutoff
    % on the neighborlist calculation the flag will go
    % to 1 signalling a recalculation of the neighborlist
    % at the beginning of the next time step.

    %%% step 2: calculate the force fields (ROTATION: CONSTRUCTION)
    Klj = Klj .* 0;        % initialize the force field
    Keb = Keb .* 0;

    for i = 1:N
        for j = 1:nlist(i)
            jj = vlist(i, j);        % neighbor indice
            

            % LJ potential
            xr = x(i) - x(jj);
            if xr > Lx*0.5
                xr = xr - Lx;
            elseif xr < -Lx*0.5
                xr = xr + Lx;
            end
            yr = y(i) - y(jj);
            if yr > Ly*0.5
                yr = yr - Ly;
            elseif yr < -Ly*0.5
                yr = yr + Ly;
            end
            r = sqrt(xr^2 + yr^2);      % particle distance

            theta = atan2(-yr, -xr) - a(i);          % theta angle in EE
            phi = atan2(yr, xr) - a(jj);         % phi angle in EE

            if theta > 2*pi
                theta = theta - 2*pi*floor(theta/(2*pi));
            elseif theta < 0
                theta = theta + 2*pi*(1+floor(-theta/(2*pi)));
            end
            if theta > pi
                theta = theta - 2*pi;
            end
            if phi > 2*pi
                phi = phi - 2*pi*floor(phi/(2*pi));
            elseif phi < 0
                phi = phi + 2*pi*(1+floor(-phi/(2*pi)));
            end
            if phi > pi
                phi = phi - 2*pi;
            end

            rri = round((theta + pi)/(2*pi/numel(Grr)));
            if rri <= 0
                rri = rri + numel(Grr);
            elseif rri > numel(Grr)
                rri = rri - numel(Grr);
            end
            rrjj = round((phi + pi)/(2*pi/numel(Grr)));
            if rrjj <= 0
                rrjj = rrjj + numel(Grr);
            elseif rrjj > numel(Grr)
                rrjj = rrjj - numel(Grr);
            end

            Grri = Grr(rri);
            Grrjj = Grr(rrjj);

            sigma = Grri + Grrjj;    % rest distance between the particle, CONSTRUCTION
            if r < rcut - 1 + sigma
                ri = 1/r;
%                 prelj = 48*ecut*ri^8*(ri^6 - 0.5);
                prelj = 48*ecut*ri^8*sigma^6*(ri^6*sigma^6 - 0.5); 
                if prelj > 1000
                    prelj = 1000;
                    warning("abnormal LJ potential");
                    count_abnormal_LJ = count_abnormal_LJ + 1;
                end
                Klj(i, 1) = Klj(i, 1) + prelj*xr;      % update pair force
                Klj(i, 2) = Klj(i, 2) + prelj*yr;
            end

            if bond(i, jj)
                xi = x(i) + Grri*cos(atan2(-yr, -xr));
                xjj = x(jj) + Grrjj*cos(atan2(yr, xr));
                xsr = xi - xjj;  
                if xsr == 0
                    xsr = 0.0001;
                elseif xsr > Lx*0.5
                    xsr = xsr - Lx;
                elseif xsr < -Lx*0.5
                    xsr = xsr + Lx;
                end
                yi = y(i) + Grri*sin(atan2(-yr, -xr));
                yjj = y(jj) + Grrjj*sin(atan2(yr, xr));
                ysr = yi - yjj; 
                if ysr == 0
                    ysr = 0.0001;
                elseif ysr > Ly*0.5
                    ysr = ysr - Ly;
                elseif ysr < -Ly*0.5
                    ysr = ysr + Ly;
                end
                rs = sqrt(xsr^2 + ysr^2);      % site distance, rs in EE

                i1a = round((theta+pi)/dd1);
                i1b = i1a - 1;
                i1c = i1a + 1;
                if i1a <= 0
                    i1a = i1a + sz(1);
                elseif i1a > sz(1)
                    i1a = i1a - sz(1);
                end
                if i1b <= 0
                    i1b = i1b + sz(1);
                elseif i1b > sz(1)
                    i1b = i1b - sz(1);
                end
                if i1c <= 0
                    i1c = i1c + sz(1);
                elseif i1c > sz(1)
                    i1c = i1c - sz(1);
                end

                i2a = round((phi+pi)/dd2);
                i2b = i2a - 1;
                i2c = i2a + 1;
                if i2a <= 0
                    i2a = i2a + sz(2);
                elseif i2a > sz(2)
                    i2a = i2a - sz(2);
                end
                if i2b <= 0
                    i2b = i2b + sz(2);
                elseif i2b > sz(2)
                    i2b = i2b - sz(2);
                end
                if i2c <= 0
                    i2c = i2c + sz(2);
                elseif i2c > sz(2)
                    i2c = i2c - sz(2);
                end

                i3a = round(rs/dd3);
                i3b = i3a - 1;
                i3c = i3a + 1;
                if i3a < 1
                    i3a = 1;
                elseif i3a > sz(3) - 1
                    i3a = sz(3) - 1;
                end
                if i3b < 1
                    i3b = 1;
                elseif i3b > sz(3) - 1
                    i3b = sz(3) - 1;
                end
                if i3c < 2
                    i3c = 2;
                elseif i3c > sz(3)
                    i3c = sz(3);
                end

                % calculate force
                f1 = -(EE(i1c, i2a, i3a) - EE(i1b, i2a, i3a))/(2*dd1);   %force in theta
                f2 = -(EE(i1a, i2c, i3a) - EE(i1a, i2b, i3a))/(2*dd2);   %force in phi
                f3 = -(EE(i1a, i2a, i3c) - EE(i1a, i2a, i3b))/(2*dd3);   %force in rs


                % update force
                Keb(i, 1) = Keb(i, 1) + f3*xsr/rs;
                Keb(i, 2) = Keb(i, 2) + f3*ysr/rs;
                Keb(i, 3) = Keb(i, 3) - f1;
                Keb(jj, 3) = Keb(jj, 3) - f2;
            end
        end
    end

    %%% step 3: update positions (ROTATION: CONSTRUCTION)
    K = Klj+Keb;
    for j = 1:N
        % Langevin
        noise_x = sqrt(2*h/(beta*gamma)) * randn(1);
        noise_y = sqrt(2*h/(beta*gamma)) * randn(1);
        noise_a = sqrt(2*h/(beta*gamma2)) * randn(1);
        x(j) = x(j) + (h/gamma)*K(j, 1) + noise_x;
        y(j) = y(j) + (h/gamma)*K(j, 2) + noise_y;
        a(j) = a(j) + (h/gamma2)*K(j, 3) + noise_a;

        % periodic boundaries
        if x(j) > Lx
            x(j) = x(j) - Lx*floor(x(j)/Lx);
        end
        if x(j) < 0
            x(j) = x(j) + Lx*(1+floor(-x(j)/Lx));
        end
        if y(j) > Ly
            y(j) = y(j) - Ly*floor(y(j)/Ly);
        end
        if y(j) < 0
            y(j) = y(j) + Ly*(1+floor(-y(j)/Ly));
        end

        if a(j) > 2*pi
            a(j) = a(j) - 2*pi*(floor(a(j)/(2*pi)));
        end
        if a(j) < 0
            a(j) = a(j) + 2*pi*(1+floor(-a(j)/(2*pi)));
        end

        xx = xv(j) - x(j);
        yy = yv(j) - y(j);
        aa = av(j) - a(j);
        if xx > Lx*0.5
            xx = xx - Lx;
        elseif xx < -Lx*0.5
            xx = xx + Lx;
        end
        if yy > Ly*0.5
            yy = yy - Ly;
        elseif yy < -Ly*0.5
            yy = yy + Ly;
        end
        rr = sqrt(xx^2 + yy^2);
        aa = sqrt(aa^2);
        if rr > rdelta/2 || aa > adelta/2
            flag = 1;
        end
    end

    %%% step 4: make/break bonds
    for m = 1:N
        nattn = 0;      % number of neighbor particles without bonds to m
        nattb = 0;      % number of neighbor bound particles to m

        for k = 1:N
            if rc(m, k) == 1 && m~=k     % neighbor particle
                if bond(m, k) == 0   % not bound
                    nattni(nattn + 1) = k;
                    nattn = nattn + 1;
                else         % bound
                    nattbi(nattb + 1) = k;
                    nattb = nattb + 1;
                end
            end
        end

        for npick = 1:nattb+nattn

            if (npick <= nattn)      % make a bond if possible
                p = nattni(npick);
                % propose a bond making
                m1 = rand(1);
                if m1 < (1 - exp(-kbond*h))        % ensure time scale
                    xr = x(m) - x(p);
                    if xr > Lx*0.5
                        xr = xr - Lx;
                    elseif xr < -Lx*0.5
                        xr = xr + Lx;
                    end
                    yr = y(m) - y(p);
                    if yr > Ly*0.5
                        yr = yr - Ly;
                    elseif yr < -Ly*0.5
                        yr = yr + Ly;
                    end

                    theta = atan2(-yr, -xr) - a(m);          % theta angle in EE
                    phi = atan2(yr, xr) - a(p);         % phi angle in EE

                    if theta > 2*pi
                        theta = theta - 2*pi*floor(theta/(2*pi));
                    elseif theta < 0
                        theta = theta + 2*pi*(1+floor(-theta/(2*pi)));
                    end
                    if theta > pi
                        theta = theta - 2*pi;
                    end
                    if phi > 2*pi
                        phi = phi - 2*pi*floor(phi/(2*pi));
                    elseif phi < 0
                        phi = phi + 2*pi*(1+floor(-phi/(2*pi)));
                    end
                    if phi > pi
                        phi = phi - 2*pi;
                    end

                    rrm = round((theta + pi)/(2*pi/numel(Grr)));
                    if rrm <= 0
                        rrm = rrm + numel(Grr);
                    elseif rrm > numel(Grr)
                        rrm = rrm - numel(Grr);
                    end
                    rrp = round((phi + pi)/(2*pi/numel(Grr)));
                    if rrp <= 0
                        rrp = rrp + numel(Grr);
                    elseif rrp > numel(Grr)
                        rrp = rrp - numel(Grr);
                    end

                    Grrm = Grr(rrm);
                    Grrp = Grr(rrp);

                    xm = x(m) + Grrm*cos(atan2(-yr, -xr));
                    xp = x(p) + Grrp*cos(atan2(yr, xr));
                    xsr = xm - xp;
                    if xsr == 0
                        xsr = 0.0001;
                    elseif xsr > Lx*0.5
                        xsr = xsr - Lx;
                    elseif xsr < -Lx*0.5
                        xsr = xsr + Lx;
                    end
                    ym = y(m) + Grrm*sin(atan2(-yr, -xr));
                    yp = y(p) + Grrp*sin(atan2(yr, xr));
                    ysr = ym - yp; 
                    if ysr == 0
                        ysr = 0.0001;
                    elseif ysr > Ly*0.5
                        ysr = ysr - Ly;
                    elseif ysr < -Ly*0.5
                        ysr = ysr + Ly;
                    end
                    rs = sqrt(xsr^2 + ysr^2);         % rs distance in EE    

                    i1 = round((theta+pi)/dd1);
                    if i1 <= 0
                        i1 = i1 + sz(1);
                    elseif i1 > sz(1)
                        i1 = i1 - sz(1);
                    end
                    i2 = round((phi+pi)/dd2);
                    if i2 <= 0
                        i2 = i2 + sz(2);
                    elseif i2 > sz(2)
                        i2 = i2 - sz(2);
                    end
                    i3 = round(rs/dd3);
                    if i3 < 1
                        i3 = 1;
                    elseif i3 > sz(3)
                        i3 = sz(3);
                    end

                    pe = EE(i1, i2, i3);
                    % accept a bond
                    b = rand(1);
%                      b = 100;
                    if b < exp(-beta*pe)/(1+exp(-beta*pe))
                        count_asso = count_asso + 1;
                        bond(m, p) = 1;
                        bond(p, m) = 1;
                    end
                end

            elseif npick <= nattb+nattn         % break a bond if possible
                p = nattbi(npick - nattn);
                % propose a bond breaking
                m1 = rand(1);
                if m1 < (1 - exp(-kbond*h))        % ensure time scale

                    xr = x(m) - x(p);
                    if xr > Lx*0.5
                        xr = xr - Lx;
                    elseif xr < -Lx*0.5
                        xr = xr + Lx;
                    end
                    yr = y(m) - y(p);
                    if yr > Ly*0.5
                        yr = yr - Ly;
                    elseif yr < -Ly*0.5
                        yr = yr + Ly;
                    end

                    theta = atan2(-yr, -xr) - a(m);          % theta angle in EE
                    phi = atan2(yr, xr) - a(p);         % phi angle in EE

                    if theta > 2*pi
                        theta = theta - 2*pi*floor(theta/(2*pi));
                    elseif theta < 0
                        theta = theta + 2*pi*(1+floor(-theta/(2*pi)));
                    end
                    if theta > pi
                        theta = theta - 2*pi;
                    end
                    if phi > 2*pi
                        phi = phi - 2*pi*floor(phi/(2*pi));
                    elseif phi < 0
                        phi = phi + 2*pi*(1+floor(-phi/(2*pi)));
                    end
                    if phi > pi
                        phi = phi - 2*pi;
                    end

                    rrm = round((theta + pi)/(2*pi/numel(Grr)));
                    if rrm <= 0
                        rrm = rrm + numel(Grr);
                    elseif rrm > numel(Grr)
                        rrm = rrm - numel(Grr);
                    end
                    rrp = round((phi + pi)/(2*pi/numel(Grr)));
                    if rrp <= 0
                        rrp = rrp + numel(Grr);
                    elseif rrp > numel(Grr)
                        rrp = rrp - numel(Grr);
                    end

                    Grrm = Grr(rrm);
                    Grrp = Grr(rrp);

                    xm = x(m) + Grrm*cos(atan2(-yr, -xr));
                    xp = x(p) + Grrp*cos(atan2(yr, xr));
                    xsr = xm - xp;
                    if xsr == 0
                        xsr = 0.0001;
                    elseif xsr > Lx*0.5
                        xsr = xsr - Lx;
                    elseif xsr < -Lx*0.5
                        xsr = xsr + Lx;
                    end
                    ym = y(m) + Grrm*sin(atan2(-yr, -xr));
                    yp = y(p) + Grrp*sin(atan2(yr, xr));
                    ysr = ym - yp;
                    if ysr == 0
                        ysr = 0.0001;
                    elseif ysr > Ly*0.5
                        ysr = ysr - Ly;
                    elseif ysr < -Ly*0.5
                        ysr = ysr + Ly;
                    end
                    rs = sqrt(xsr^2 + ysr^2);        % rs distance in EE

                    i1 = round((theta+pi)/dd1);
                    if i1 <= 0
                        i1 = i1 + sz(1);
                    elseif i1 > sz(1)
                        i1 = i1 - sz(1);
                    end
                    i2 = round((phi+pi)/dd2);
                    if i2 <= 0
                        i2 = i2 + sz(2);
                    elseif i2 > sz(2)
                        i2 = i2 - sz(2);
                    end
                    i3 = round(rs/dd3);
                    if i3 < 1
                        i3 = 1;
                    elseif i3 > sz(3)
                        i3 = sz(3);
                    end

                    pe = EE(i1, i2, i3);
                    % break a bond
                    b = rand(1);
%                      b = 100;
                    if b < exp(beta*pe)/(1+exp(beta*pe))
                        count_disso = count_disso + 1;
                        bond(m, p) = 0;
                        bond(p, m) = 0;
                    end
                end
            end
        end

     end

    %%% step 5: write files
    write = rem(f,icr) == 0;
    if write
        A = [ones(N,1)*f, (1:N)', x, y, a];
        fprintf(fileID_xya,'%10d %10d %10.6f %10.6f %10.6f\n', A');
        for k = 1:N-1
            for j = k+1:N
                if bond(k, j) == 1
                    fprintf(fileID_hb,'%10d %10d %10d\n', [f k j]);
                end
            end
        end
        write = 0;
    end
end

fclose(fileID_xya);
fclose(fileID_hb);

fprintf(fileID_info,'%10s %1s\n', string(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')), 'finish');
if count_abnormal_LJ > 0
    fprintf(fileID_info,'%1s: %10d\n','abnormal LJ potential', count_abnormal_LJ);
end

fclose(fileID_info);







