cdir = pwd(); ii = findlast("POD", cdir); PODpath = cdir[1:ii[end]] * "/";    
push!(LOAD_PATH, PODpath * "src");
using POD

# path to LAMMPS executive 
lammps = "/Users/ngoccuongnguyen/lammps/build/lmp";

# fit POD potential
POD.fitpod(lammps, "fpod_param.pod", "data.pod");

