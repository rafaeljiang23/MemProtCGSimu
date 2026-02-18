%% simulation setup
%%% Note that these parameters MUST match defined values in mem_simu mat
%%% files (i.e. mem_simu_v10b.mat)

% rmin = 0.1;      % minimal distance between sites, defined in run_MemCGSimu_Assembly.m
% rbond = 0.9;     % reaction distance, defined in run_MemCGSimu_Assembly.m

%% tetramer geometry_trpv3_v2

dd1 = 0.05;
dd2 = 0.05;
dd3 = 0.01;

theta = [0 ...
    pi/4 ...
    -pi/4];
phi = -theta;
rs = [4.4*rmin ...
    rmin ...
   rmin];


ee = [-5 ...
    -25 ...
    -25];


kappa = [50 50 100; ...
    100 100 200;...
    100 100 200];

EEinit = zeros(ceil(2*pi/dd1), ceil(2*pi/dd2),  ceil(2*rbond/dd3));
sz = size(EEinit);
num_well = numel(ee);
values = zeros(num_well, 1);

for ii = 1:sz(1)
    for jj = 1:sz(2)
        for kk = 1:sz(3)
            values = values .* 0;
            for aa = 1:num_well
                dii_min = sz(1);
                djj_min = sz(2);
                for rr = 0
                    dii = ii + rr*sz(1) - (theta(aa) + pi)/dd1;
                    djj = jj + rr*sz(2) - (phi(aa) + pi)/dd2;
                    if dii^2 < dii_min^2
                        dii_min = dii;
                    end
                    if djj^2 < djj_min^2
                        djj_min = djj;
                    end
                end
                dii = dii_min;
                djj = djj_min;
                dkk = kk - rs(aa)/dd3;

                values(aa) = kappa(aa, 1)*(dii*dd1*cos(pi/4) - djj*dd2*sin(pi/4))^2 +...
                    kappa(aa, 2)*(dii*dd1*cos(pi/4) + djj*dd2*sin(pi/4))^2 +...
                    kappa(aa, 3)*(dkk*dd3)^2 + ee(aa);

            end
          EEinit(ii, jj, kk) = -log(sum(exp(-values)));
        end
    end
end
EEinit_tet = EEinit * 1;


%% pentamer 
dd1 = 0.05;
dd2 = 0.05;
dd3 = 0.01;

theta = [ ...
    pi/10 3*pi/10 ...
    -pi/10 -3*pi/10];
phi = -theta;
rs = [...
    10.4*rmin 1.9*rmin ...
   10.4*rmin 1.9*rmin];



ee = [...
    -4.5 -24 ...
    -4.5 -24];


kappa = [ ...
    50 50 100; 100 100 200;...
    50 50 100; 100 100 200;];


EEinit = zeros(ceil(2*pi/dd1), ceil(2*pi/dd2),  ceil(2*rbond/dd3));
sz = size(EEinit);
num_well = numel(ee);
values = zeros(num_well, 1);

for ii = 1:sz(1)
    for jj = 1:sz(2)
        for kk = 1:sz(3)
            values = values .* 0;
            for aa = 1:num_well
                dii_min = sz(1);
                djj_min = sz(2);
                for rr = 0
                    dii = ii + rr*sz(1) - (theta(aa) + pi)/dd1;
                    djj = jj + rr*sz(2) - (phi(aa) + pi)/dd2;
                    if dii^2 < dii_min^2
                        dii_min = dii;
                    end
                    if djj^2 < djj_min^2
                        djj_min = djj;
                    end
                end
                dii = dii_min;
                djj = djj_min;
                dkk = kk - rs(aa)/dd3;

                values(aa) = kappa(aa, 1)*(dii*dd1*cos(pi/4) - djj*dd2*sin(pi/4))^2 +...
                    kappa(aa, 2)*(dii*dd1*cos(pi/4) + djj*dd2*sin(pi/4))^2 +...
                    kappa(aa, 3)*(dkk*dd3)^2 + ee(aa);

            end
          EEinit(ii, jj, kk) = -log(sum(exp(-values)));

        end
    end
end
EEinit_pen = EEinit * 1;

%% connection path
dd1 = 0.05;
dd2 = 0.05;
dd3 = 0.01;

theta1 = [0 0];
phi1 = -theta1;
rs1 = [4.4*rmin 4.4*rmin];

theta2 = [pi/10 -pi/10];
phi2 = -theta2;
rs2 = [10.4*rmin 10.4*rmin];

ee0 = [-1 -1];

kappa0 = [50 50 100; ...
    50 50 100];
ll = [10 10];

theta = [];
phi = [];
rs = [];
ee = [];
kappa = [];
for i = 1:numel(ll)
    theta = [theta linspace(theta1(i), theta2(i), ll(i))];
    phi = [phi linspace(phi1(i), phi2(i), ll(i))];
    rs = [rs linspace(rs1(i), rs2(i), ll(i))];
    ee = [ee ones(1, ll(i)).*ee0(i)];
    kappa00 = ones(ll(i), 3);
    for j = 1:ll(i)
        kappa00(j, :) = kappa0(i, :);
    end
    kappa = [kappa; kappa00];
end


EEinit = zeros(ceil(2*pi/dd1), ceil(2*pi/dd2),  ceil(2*rbond/dd3));
sz = size(EEinit);
num_well = numel(ee);
values = zeros(num_well, 1);

for ii = 1:sz(1)
    for jj = 1:sz(2)
        for kk = 1:sz(3)
            values = values .* 0;
            for aa = 1:num_well
                dii_min = sz(1);
                djj_min = sz(2);
                for rr = 0
                    dii = ii + rr*sz(1) - (theta(aa) + pi)/dd1;
                    djj = jj + rr*sz(2) - (phi(aa) + pi)/dd2;
                    if dii^2 < dii_min^2
                        dii_min = dii;
                    end
                    if djj^2 < djj_min^2
                        djj_min = djj;
                    end
                end
                dii = dii_min;
                djj = djj_min;
                dkk = kk - rs(aa)/dd3;

                values(aa) = kappa(aa, 1)*(dii*dd1*cos(pi/4) - djj*dd2*sin(pi/4))^2 +...
                    kappa(aa, 2)*(dii*dd1*cos(pi/4) + djj*dd2*sin(pi/4))^2 +...
                    kappa(aa, 3)*(dkk*dd3)^2 + ee(aa);

            end
           EEinit(ii, jj, kk) = -log(sum(exp(-values)));

        end
    end
end
EEinit_p = EEinit * 1;

%% integrate energy landspace
EEinit = min(EEinit_pen, EEinit_tet);
EEinit = min(EEinit, EEinit_p);
