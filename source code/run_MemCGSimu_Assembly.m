%%% Windows MacOS
%%% template for patch simulations for testing EE for tetramer-pentamer transitions
%%% mem_simu_v10b.m

%% simulation parameters
h = 0.001;      % unit: s, delta time
Nsteps =  100000;       % unit: s, total number of time-steps
                        % for test purpose, note that it takes longer time
                        % for pentamer to emerge （oligomer distribution
                        % stablizes >1000s)
rmin = 0.1;      % minimal distance between sites
rbond = 0.9;     % reaction distance


%% set up simulations
cdadress = cd;
diradress = [];
filenames = [];

%%% energy landscape
files_EE = "init_EE_geo_v2_transition_test63c.mat";   % pre-defined for TRPV3 tetramer/pentamer 
                                                      % transition simulation
run("create_EE_transition_geo.m");
save (files_EE, "EEinit", "EEinit_tet", "EEinit_pen", "EEinit_p"); 

%%% simulation setup
files_coors = "init_coor_geo_v2_N64_L32.mat";   % pre-defined for ~50% membrane coverage

%%% define filename substrings for saving
filessave = "test63c_1";    % energy landspace name
filetitleflags = ["transition_tetra", "geo_v2"];    % other flags

%% run simulations
seed = 1;
it = 1;
%%% load parameters
load(files_EE)
load(files_coors)
%%% core simulations
run("mem_simu_v10b.m")
%%% save intial conditions
filetitlesave = filetitleflags(1) + "_N" + num2str(Ninit/4) + "_L" +...
    num2str(Lxinit) + "_" + filetitleflags(2) + "_" + filessave;
save("init_" + filetitlesave + ".mat", "bondinit", "ainit", "xinit", "yinit", "EEinit");
cd(cdadress)
%%% update initial conditions
run("supplement.m")
cd(dirname)
save("init_" + filetitlesave + "_" + num2str(Nsteps*h) + "s" + ".mat", "bondinit", "ainit", "xinit", "yinit", "EEinit");
diradress{it} = dirname;
filenames{it} = filename;
%%% get ready for the next
cd(cdadress)

filetitlesave = filetitleflags(1) + "_N" + num2str(Ninit/4) + "_L" +...
    num2str(Lxinit) + "_" + filetitleflags(2);
save(filetitlesave + ".mat", "diradress", "filenames");
%% run analysis
cd(cdadress)
run("analyze_v3_init.m")
Grr = geometry_trpv3_v2(1000);
load(files_coors);
dirname = diradress{it};
copyfile("analyze_v5.m", dirname);
cd(dirname)
filename = filenames{it};
run("analyze_v5.m")
%%%
filetitlesave = filetitleflags(1) + "_N" + num2str(Ninit/4) + "_L" +...
    num2str(Lxinit) + "_" + filetitleflags(2)  + "_" + filessave;
save("display_" + filetitlesave + ".mat", "movie_xy", "movie_hb", "movie_a",...
    "movie_xy2", "movie_xy3", "movie_xy4", "movie_hb2", "movie_a2", "-v7.3")
save("display2_" + filetitlesave + ".mat", "movie_hb_olig", ...
    "movie_hb_olig2", "movie_hb_olig_nn", "movie_hb_olig_nn2", "stat_olig", ...
    "-v7.3");
save("init_analyze_" + filetitlesave + ".mat", "oligomers_stat_init", "oligomers_init", "f_init"); 

%%%
delete analyze_v5.m
cd(cdadress)

%% display simulation results
%%% must install MIJ %%%
%%% display protomer (position + contour + area)
MIJ.createImage(movie_xy2 + movie_xy3 + movie_xy4);

%%% display interfacial harmonic bond
MIJ.createImage(movie_hb_olig_nn2);
