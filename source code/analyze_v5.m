%%%%
% Simulation analysis and display
% work with analyze_v3_init.m
%%%%

%% create analysis text files

filename_info = filename + '_info.txt';
filename_xya = filename + '_xya.txt';
filename_hb = filename + '_harmonic.txt';
fileID_info = fopen(filename_info,'r');
fileID_xya = fopen(filename_xya,'r');
fileID_hb = fopen(filename_hb,'r');
filename_olig = filename + '_olig.txt';

%% define useful parameters for visualizing simulation
Lx = Lxinit;
Ly = Lyinit;
N = Ninit;

%%% geometry funciton 
% Grr = geometry_trpv3_v2(1000);  % default

%%% simulation analysis setup
jump = 1;   % analysis step
samples = 0:jump*icr:Nsteps;
Nsample = numel(samples);

%%% simulation display setup
scale = 8;   % display scale
jumpf = 10*jump;   % display step
frames = 0:jumpf*icr:Nsteps;
Nframe = numel(frames);

%% define active protomer recorder
% mostly defined by analyze_v3_init.m
% oligomer ID x N protomer, recording active oligomer connectivity
oligomers = zeros(N*2, N);

% oligomer ID x 6 column (start time, dwell time, oligomer state,
% #interfaceial bond, #code, #sub-code), recording active oligomer info
% time is in time-step
oligomers_stat = zeros(N*2, 6);

oligocounter = 0;  % total oligomer counter
activeoligocounter = 0;  % active oligomer counter
protflag = zeros(N, 1);   % oligomer state tracker
f0 = 0;

%% initial condition (analyze_v3_init.m)
f0 = f_init;
oligomers = oligomers_init;
oligomers_stat = oligomers_stat_init;
oligocounter = max(oligomers_stat(: , 5));
activeoligocounter = find(oligomers_stat(:, 2) == 0);
activeoligocounter = activeoligocounter(1);

%% set display matrices
movie_hb_olig = zeros(Lx*scale, Ly*scale, Nframe);   % harmonic bond
movie_hb_olig_nn = zeros(Lx*scale, Ly*scale, Nframe);  % interfacial harmonic
                                                       % bond (between neighbors in oligomers)
movie_xy = zeros(Lx*scale, Ly*scale, Nframe);   % protomer position xy
movie_xy3 = movie_xy.*0;   % protomer contour
movie_xy4 = movie_xy.*0;   % protomer area
movie_a = zeros(Lx*scale, Ly*scale, Nframe);   % protomer position a
movie_hb = zeros(Lx*scale, Ly*scale, Nframe);   % harmonic bond at the patch
stat_olig = zeros(Nsample, 5);   % protomer stats (as in different oligomeic
                                 % states: mono, di, tri, tetra, penta

%% reading files
coors = zeros(N, 6);
bonds = zeros(N, N);

sc = 0;
rho = linspace(-pi, pi, 1001);
rho = rho(2:end);

[nn, ~] = size(oligomers_stat);
[~] = fgetl(fileID_xya);        % first line
[~] = fgetl(fileID_hb);        % first line

byt_hb = 0;

size(oligomers_stat)
fileID_olig = fopen(filename_olig,'w');
fprintf(fileID_olig,'%1s %1s %1s %1s %1s %1s %1s %1s\n','oligo ID',...
    'oligo ID (sub)','start time-step','final time-step',...
    'dwell time (s)', 'oligo state', 'interface bond','particles');

it1 = 1;
it2 = 1;
while ~feof(fileID_xya)

    %%% reading xy files
    s = str2double(fscanf(fileID_xya, '%s', 1));
    p = str2double(fscanf(fileID_xya, '%s', 1));
    x = str2double(fscanf(fileID_xya, '%s', 1));
    y = str2double(fscanf(fileID_xya, '%s', 1));
    a = str2double(fscanf(fileID_xya, '%s', 1));

    %%% target time step found
    %     if ismember(s, samples)
    if it1 > length(samples)
        % skip
    elseif s == samples(it1)
        % calculate site coordinates

        rri = round(numel(Grr)/2);
        Grr0 = Grr(rri);
        xs = x + Grr0*cos(a);
        ys = y + Grr0*sin(a);

        if xs <= 0
            xs = xs + Lx;
        elseif xs > Lx
            xs = xs - Lx;
        end
        if ys <= 0
            ys = ys + Ly;
        elseif ys > Ly
            ys = ys - Ly;
        end

        %%% update coordinate file
        coors(p, 1) = p;
        coors(p, 2) = x;
        coors(p, 3) = y;
        coors(p, 4) = a;
        coors(p, 5) = xs;
        coors(p, 6) = ys;

        %%% coordinate file update finished and start read bond files
        if p == N
            sc = s;
            %             f = find(samples == sc) + f0
            f = it1 + f0;
            f*h*icr*jump
            it1 = it1 + 1;
            wcount = 0;    % counter for writing files
            A = [];
            Aop = [];

            %%% draw particle positions
            %             if ismember(sc, frames)
            if it2 > length(frames)
                % skip
            elseif sc == frames(it2)
                % Do NOT UPDATE it2
                %                 f2 = find(frames == sc);
                f2 = it2;
                for i = 1:N
                    %
                    x = coors(i, 2);
                    y = coors(i, 3);
                    %
                    xx = round(x*scale);
                    yy = round(y*scale);
                    if xx <= 0
                        xx = xx + Lx*scale;
                    elseif xx > Lx*scale
                        xx = xx - Lx*scale;
                    end
                    if yy <= 0
                        yy = yy + Ly*scale;
                    elseif yy > Ly*scale
                        yy = yy - Ly*scale;
                    end

                    movie_xy(xx, yy, f2) = 1;


                    a = coors(i, 4);

                    for j = 1:numel(rho)
                        rri = round((rho(j) - a + pi)/(2*pi/numel(Grr)));
                        if rri <= 0
                            rri = rri + numel(Grr);
                        elseif rri > numel(Grr)
                            rri = rri - numel(Grr);
                        end
                        Grri = Grr(rri);
                        x2 = x + Grri*cos(rho(j));
                        y2 = y + Grri*sin(rho(j));
                        xx2 = round(x2*scale);
                        yy2 = round(y2*scale);
                        if xx2 <= 0
                            xx2 = xx2 + Lx*scale;
                        elseif xx2 > Lx*scale
                            xx2 = xx2 - Lx*scale;
                        end
                        if yy2 <= 0
                            yy2 = yy2 + Ly*scale;
                        elseif yy2 > Ly*scale
                            yy2 = yy2 - Ly*scale;
                        end
                        movie_xy3(xx2, yy2, f2) = 1;

                        x3 = linspace(x, x2, Grri*scale);
                        y3 = linspace(y, y2, Grri*scale);
                        for k = 1:numel(x3)
                            xx3 = round(x3(k)*scale);
                            yy3 = round(y3(k)*scale);
                            if xx3 <= 0
                                xx3 = xx3 + Lx*scale;
                            elseif xx3 > Lx*scale
                                xx3 = xx3 - Lx*scale;
                            end
                            if yy3 <= 0
                                yy3 = yy3 + Ly*scale;
                            elseif yy3 > Ly*scale
                                yy3 = yy3 - Ly*scale;
                            end
                            movie_xy4(xx3, yy3, f2) = 1;
                        end
                    end

                    xs = coors(i, 5);
                    ys = coors(i, 6);
                    xxa = round(xs*scale);
                    yya = round(ys*scale);
                    if xxa <= 0
                        xxa = xxa + Lx*scale;
                    elseif xxa > Lx*scale
                        xxa = xxa - Lx*scale;
                    end
                    if yya <= 0
                        yya = yya + Ly*scale;
                    elseif yya > Ly*scale
                        yya = yya - Ly*scale;
                    end
                    movie_a(xxa, yya, f2) = 1;

                end
            end

            %%% read harmonic bond files
            while ~feof(fileID_hb)
                s = str2double(fscanf(fileID_hb, '%s', 1));
                p1 = str2double(fscanf(fileID_hb, '%s', 1));
                p2 = str2double(fscanf(fileID_hb, '%s', 1));
                if s == sc
                    %%% draw bond positions
                    if it2 > length(frames)
                        % skip
                    elseif sc == frames(it2)
                        % Do NOT UPDATE it2
                        f2 = it2;
                        x1 = coors(p1, 2);
                        x2 = coors(p2, 2);
                        y1 = coors(p1, 3);
                        y2 = coors(p2, 3);
                        a1 = coors(p1, 4);
                        a2 = coors(p2, 4);
                        xs1 = coors(p1, 5);
                        xs2 = coors(p2, 5);
                        ys1 = coors(p1, 6);
                        ys2 = coors(p2, 6);
                        x = 0.5*(xs1 + xs2);
                        y = 0.5*(ys1 + ys2);
                        xx = round(x*scale);
                        yy = round(y*scale);
                        if xx <= 0
                            xx = xx + Lx*scale;
                        elseif xx > Lx*scale
                            xx = xx - Lx*scale;
                        end
                        if yy <= 0
                            yy = yy + Ly*scale;
                        elseif yy > Ly*scale
                            yy = yy - Ly*scale;
                        end
                        % draw lines
                        movie_hb(xx, yy, f2) = 1;
                    end

                    bonds(p1, p2) = 1;
                    bonds(p2, p1) = 1;
                elseif s > sc
                    byt_hb = ftell(fileID_hb) - 30;
                    frewind(fileID_hb);
                    fseek(fileID_hb, byt_hb, "bof");
                    break;
                end
            end

            %%% find bonds holding oligomers
            bonds2 = bonds;      % oligomer bonds, keep updated
            for p1 = 1:N
                pp = find(bonds2(p1, :));
                coors1 = coors(p1, :);
                for i = 1:numel(pp)
                    p2 = pp(i);
                    coors2 = coors(p2, :);
                    xr = coors1(2) - coors2(2);
                    if xr > Lx*0.5
                        xr = xr - Lx;
                    elseif xr < -Lx*0.5
                        xr = xr + Lx;
                    end
                    yr = coors1(3) - coors2(3);
                    if yr > Ly*0.5
                        yr = yr - Ly;
                    elseif yr < -Ly*0.5
                        yr = yr + Ly;
                    end

                    r = sqrt(yr^2 + xr^2);

                    xsr = coors1(5) - coors2(5);
                    if xsr > Lx*0.5
                        xsr = xsr - Lx;
                    elseif xsr < -Lx*0.5
                        xsr = xsr + Lx;
                    end
                    ysr = coors1(6) - coors2(6);
                    if ysr > Ly*0.5
                        ysr = ysr - Ly;
                    elseif ysr < -Ly*0.5
                        ysr = ysr + Ly;
                    end

                    rs = sqrt(ysr^2 + xsr^2);
                    theta = atan2(-yr, -xr) - coors1(4);
                    phi = atan2(yr, xr) - coors2(4);

                    if theta > 2*pi
                        theta = theta - 2*pi*floor(theta/(2*pi));
                    elseif theta < 0
                        theta = theta + 2*pi*(1+floor(-theta/(2*pi)));
                    end

                    if phi > 2*pi
                        phi = phi - 2*pi*floor(phi/(2*pi));
                    elseif phi < 0
                        phi = phi + 2*pi*(1+floor(-phi/(2*pi)));
                    end

                    % theta + phi should be around 2pi for oligomeric bonds
                    % default allowance: 36 degree (10%)

                    % ra should be < r/2. NOTICE: this condition will
                    % likely remove hexameric bond and completely remove
                    % higher oligomeric bond
                    % default allowance: 10%

                    if abs(theta + phi - 2*pi) > pi/4 || ...
                            rs > 3*r/4
                        % condition not met, invalid oligomeric bond
                        bonds2(p1, p2) = 0;
                        bonds2(p2, p1) = 0;
                    end
                end
            end

            %%% finding oligomer connectivity
            conns = bonds2;
            del = numel(conns);
            while del > 0
                conns0 = conns;
                for p1 = 1:N
                    pp = find(conns(p1, :));
                    for i = 1:numel(pp)
                        p2 = pp(i);
                        conns(p1, :) = conns(p1, :) + conns(p2, :);
                    end
                    conns(p1, p1) = 1;
                    conns = conns + conns';
                end
                conns = double(conns > 0);
                del = sum(conns ~= conns0, "all");
            end

            %%% find the interfacial bonds between oligomers
            bonds3 = bonds2.*0;      % interfacial oligomer bonds, keep updated
            for p1 = 1:N
                pp = find(conns(p1, :));
                ppa = coors(pp, 4);
                [ppas, pps] = sort(ppa);    % ordered protomers in oligomer
                p1s = find(pps == find(pp == p1));
                a1 = coors(p1, 4);
                if numel(pps) == 1
                    % no neighboring protomer
                elseif numel(pps) == 2
                    % one neighboring protomer
                    p2s = pps(3-p1s);
                    p2 = pp(p2s);
                    bonds3(p1, p2) = 1;
                    bonds3(p2, p1) = 1;
                else
                    % more neighboring protomer
                    if p1s == 1
                        p2s = pps(p1s+1);
                        p3s = pps(end);
                    elseif p1s == numel(pps)
                        p2s = pps(1);
                        p3s = pps(p1s-1);
                    else
                        p2s = pps(p1s+1);
                        p3s = pps(p1s-1);
                    end
                    p2 = pp(p2s);
                    p3 = pp(p3s);
                    a2 = a1 - ppa(p2s);
                    if a2 > 2*pi
                        a2 = a2 - 2*pi*floor(a2/(2*pi));
                    elseif a2 < 0
                        a2 = a2 + 2*pi*(1+floor(-a2/(2*pi)));
                    end

                    a3 = a1 - ppa(p3s);
                    if a3 > 2*pi
                        a3 = a3 - 2*pi*floor(a3/(2*pi));
                    elseif a3 < 0
                        a3 = a3 + 2*pi*(1+floor(-a3/(2*pi)));
                    end

                    if abs(a2 + a3 - 2*pi) < pi/4
                        bonds3(p1, p2) = 1;
                        bonds3(p2, p1) = 1;
                        bonds3(p1, p3) = 1;
                        bonds3(p3, p1) = 1;
                    end
                end
            end
            bonds3 = bonds3 .* bonds2;

            %%% finding oligomer connectivity
            conns2 = bonds3;
            del = numel(conns2);
            while del > 0
                conns0 = conns2;
                for p1 = 1:N
                    pp = find(conns2(p1, :));
                    for i = 1:numel(pp)
                        p2 = pp(i);
                        conns2(p1, :) = conns2(p1, :) + conns2(p2, :);
                    end
                    conns2(p1, p1) = 1;
                    conns2 = conns2 + conns2';
                end
                conns2 = double(conns2 > 0);
                del = sum(conns2 ~= conns0, "all");
            end

            %%% display oligomer bonds and interfacial oligomer bonds
            if it2 > length(frames)
                % skip
            elseif sc == frames(it2)
                f2 = it2;
                % UPDATE it2 at its last call in the while loop
                it2 = it2 + 1;
                for p1 = 1:N-1
                    for p2 = p1+1:N
                        if bonds2(p1, p2)
                            x1 = coors(p1, 2);
                            x2 = coors(p2, 2);
                            y1 = coors(p1, 3);
                            y2 = coors(p2, 3);
                            %
                            if x1 < x2 && x1 + Lx - x2 < x2 - x1
                                x1 = x1 + Lx;
                            elseif x1 > x2 && x2 + Lx - x1 < x1 - x2
                                x2 = x2 + Lx;
                            end
                            if y1 < y2 && y1 + Ly - y2 < y2 - y1
                                y1 = y1 + Ly;
                            elseif y1 > y2 && y2 + Ly - y1 < y1 - y2
                                y2 = y2 + Ly;
                            end
                            x = 0.5*(x1 + x2);
                            y = 0.5*(y1 + y2);
                            xx = round(x*scale);
                            yy = round(y*scale);
                            if xx <= 0
                                xx = xx + Lx*scale;
                            elseif xx > Lx*scale
                                xx = xx - Lx*scale;
                            end
                            if yy <= 0
                                yy = yy + Ly*scale;
                            elseif yy > Ly*scale
                                yy = yy - Ly*scale;
                            end

                            movie_hb_olig(xx, yy, f2) = 1;

                            if bonds3(p1, p2)
                                movie_hb_olig_nn(xx, yy, f2) = 1;
                            end

                        end

                    end
                end
            end

            %%% analyze oligomer states
            for op = 1:nn
                if oligomers_stat(op, 2) > 0     % active oligomer
                    wflag = 1;
                    pp0 = find(oligomers(op, :));   % active oligomer connectivity
                    oligo0 = oligomers_stat(op, 3);   % active oligomer state
                    if numel(pp0) ~= oligo0
                        warning("oligomer state error");
                    end
                    pp1 = find(conns2(pp0(1), :));   % new oligomer connectivity (first particle)
                    oligo1 = numel(pp1);
                    oligo0b3 = oligomers_stat(op, 4);   % active oligomer interfacial bonds
                    oligo1b3 = numel(find(bonds3(pp0, :)))/2;

                    if oligo1 <= 5
                        stat_olig(f-f0, oligo1) = stat_olig(f-f0, oligo1) + 1;
                    end

                    if sum(ismember(pp0, pp1)) == oligo0 &&...
                            sum(ismember(pp1, pp0)) == oligo1
                        % no change of active oligomer connectivity
                        if oligo0b3 == oligo1b3
                            % no change of active oligomer interfacial bond
                            % update dwell time
                            wflag = 0;
                            oligomers_stat(op, 2) = oligomers_stat(op, 2) + 1;
                            protflag(pp1) = 1;
                        else
                            % change of active oligomer interfacial bond
                            % create new active oligomer with the same code
                            % but different subcode

                            wflag = 0;
                            oligomers_stat(op, 2) = oligomers_stat(op, 2) + 1;
                            protflag(pp1) = 1;

                        end
                    end
                    if wflag
                        % active oligomer connectivity broken
                        % remove broken oligomer connectivity
                        if ~ismember(op, Aop)
                            wcount = wcount + 1;
                            Aop(wcount) = op;
                            % collection of files to be written
                            A{wcount, 1} = [oligomers_stat(op, 5) oligomers_stat(op, 6)...
                                oligomers_stat(op, 1) f jump*icr*h*oligomers_stat(op, 2)...
                                oligomers_stat(op, 3) oligomers_stat(op, 4)];
                            A{wcount, 2} = pp0;
                        end
                    end
                end
            end

            %%% new connectivity
            for p1 = 1:N
                if ~protflag(p1)    % oligomer state not cleared
                    pp1 = find(conns2(p1, :));
                    oligo1 = numel(pp1);
                    if oligo1 == 6
                        warning("hexamer found");
                    end
                    oligo1b3 = numel(find(bonds3(pp1, :)))/2;

                    oligocounter = oligocounter + 1;
                    while oligomers_stat(activeoligocounter, 2) ~= 0
                        % update active oligomer counter
                        activeoligocounter = activeoligocounter + 1;
                        if activeoligocounter > nn
                            activeoligocounter = 1;
                        end
                    end
                    oligomers(activeoligocounter, pp1) = 1;
                    oligomers_stat(activeoligocounter, 1) = f;
                    oligomers_stat(activeoligocounter, 2) = 1;
                    oligomers_stat(activeoligocounter, 3) = oligo1;
                    oligomers_stat(activeoligocounter, 4) = oligo1b3;
                    oligomers_stat(activeoligocounter, 5) = oligocounter;
                    oligomers_stat(activeoligocounter, 6) = 0;
                    protflag(pp1) = 1;
                    if oligo1 <= 5
                        stat_olig(f-f0, oligo1) = stat_olig(f-f0, oligo1) + 1;
                    end

                end
            end


            %%% write files
            if wcount > 0
                for w = 1:wcount

                    %%% write files
                    op = Aop(w);
                    AA = A{w, 1};
                    App = A{w, 2};
                    fprintf(fileID_olig, '%10d %10d %10d %10d %10.6f %10d %10d', AA);
                    for p0 = 1:numel(App)
                        fprintf(fileID_olig, '%10d', App(p0));
                    end
                    fprintf(fileID_olig, '\n');

                    %%% inactivate oligomers
                    oligomers(op, :) = oligomers(op, :) .* 0;
                    oligomers_stat(op, :) = oligomers_stat(op, :) .* 0;
                end
            end

            %%% initilize new coordinate file
            coors = coors.*0;
            bonds = bonds.*0;
            protflag = protflag.*0;
        end
    end
end


fclose(fileID_olig);
%% final processing

sz = size(movie_hb_olig);
movie_hb_olig2 = movie_hb_olig.*0;      % harmonic bond (for display)
movie_hb_olig_nn2 = movie_hb_olig_nn.*0;      % interfacial harmonic bond (for display)
movie_xy2 = movie_xy.*0;      % protomer position xy (for display)
movie_hb2 = movie_hb.*0;      % harmonic bond (for display)
movie_a2 = movie_a.*0;      % protomer position a (for display)
for i = 1:sz(3)
    movie_hb_olig2(:, :, i) = 0.3*scale - bwdist(squeeze(movie_hb_olig(:, :, i)));
    movie_hb_olig_nn2(:, :, i) = 0.3*scale - bwdist(squeeze(movie_hb_olig_nn(:, :, i)));
    movie_xy2(:, :, i) = 0.2*scale - bwdist(squeeze(movie_xy(:, :, i)));
    movie_hb2(:, :, i) = 0.12*scale - bwdist(squeeze(movie_hb(:, :, i)));
    movie_a2(:, :, i) = 0.12*scale - bwdist(squeeze(movie_a(:, :, i)));
end

movie_hb_olig2 = movie_hb_olig2 .* (movie_hb_olig2 > 0);
movie_hb_olig_nn2 = movie_hb_olig_nn2 .* (movie_hb_olig_nn2 > 0);
movie_xy2 = double(movie_xy2 > 0);
movie_hb2 = double(movie_hb2 > 0);
movie_a2 = double(movie_a2 > 0);

%%
fclose(fileID_xya);
fclose(fileID_hb);
fclose(fileID_info);
%% update initial condition

oligomers_stat_init = oligomers_stat;
oligomers_init = oligomers;
f_init = f;

"updated"