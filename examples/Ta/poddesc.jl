cdir = pwd(); ii = findlast("POD", cdir); PODpath = cdir[1:ii[end]] * "/";    
push!(LOAD_PATH, PODpath * "src");
using POD

# path to LAMMPS executive 
lammps = "/Users/ngoccuongnguyen/lammps/build/lmp";
pod = "fpod_param.pod";
datapath = "training";

# compute POD descriptors for the training dataset 
gd, gdd, na = POD.calculatedescriptors(lammps, pod, datapath);



