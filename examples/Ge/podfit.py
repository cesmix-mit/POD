# import external modules
import numpy, os, sys

# Add POD to Python search path
cdir = os.getcwd(); ii = cdir.find("POD");
srcdir = cdir[0:(ii+3)] + "/"  + "/src";
sys.path.append(srcdir);

# import POD module
import POD

# path to LAMMPS executive 
lammps = "/Users/ngoccuongnguyen/lammps/build/lmp";

# fit POD potential
POD.fitpod(lammps, "fpod_param.pod", "data.pod");
