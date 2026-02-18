N = Ninit;

% oligomer ID x N protomer, recording active oligomer connectivity
oligomers_init = zeros(N*10, N);   

% oligomer ID x 2 column (start time, dwell time, #bond, #code), recording active oligomer info
% time is in time-step
oligomers_stat_init = zeros(N*2, 6);   

oligocounter = 0;  % oligomer counter
protflag = zeros(N, 1);   % recording if the oligomer state is clear

bonds = bondinit;
% find all connectivity
for i = 1:N
    pp = find(bonds(i, :));
    for j = 1:numel(pp)
        bonds(i, :) = bonds(i, :) + bonds(pp(j), :);
    end
    bonds(i, i) = 1;
end
bonds = double(bonds > 0);

for p = 1:N
    if ~protflag(p)
        pp = find(bonds(p, :) == 1);
        oligocounter = oligocounter + 1;
        oligomers_init(oligocounter, pp) = 1;
        oligomers_stat_init(oligocounter, 1) = 0;
        oligomers_stat_init(oligocounter, 2) = 1;
        oligomers_stat_init(oligocounter, 3) = numel(pp);
        oligomers_stat_init(oligocounter, 4) = numel(pp);
        oligomers_stat_init(oligocounter, 5) = oligocounter;
        oligomers_stat_init(oligocounter, 6) = 0;
        protflag(pp) = 1;
    end
end

f_init = 0;