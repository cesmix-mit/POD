import os
import numpy 

def calculatedescriptors(lammps, podfile, datapath):

  # make a compute podfile file
  text = "file_format extxyz";
  text = text + "\n" + "file_extension xyz";
  text = text + "\n" + "path_to_training_data_set " + datapath;
  text = text + "\n" + "compute_pod_descriptors 2";

  fid = open("compute.pod", "w")
  fid.write(text)
  fid.close()

  # make a fit podfile file
  text = "units metal";
  text = text + "\n" + "fitpod " + podfile + " compute.pod";
  fid = open("fit.pod", "w")
  fid.write(text)
  fid.close()

  # run LAMMPS to calculate POD descriptors
  runstr = lammps + " -in fit.pod";  
  os.system(runstr);

  # determine the number of configurations in the dataset
  k = 0;
  while (1 > 0):
    filename = datapath + "/descriptors_config" + str(k+1) + ".bin";
    if os.path.isfile(filename):
      k = k + 1;
    else:
      break;

  num_configs = k;

  gd = []
  gdd = []
  na = numpy.arange(num_configs)
  for i in range(0,num_configs):
    filename = datapath + "/descriptors_config" + str(i+1) + ".bin";
    tmp = numpy.fromfile(open(filename, "r"), dtype=numpy.float64);
    n =  numpy.size(tmp);
    natom = numpy.int64(tmp[0]);
    na[i] = natom
    nd = numpy.int64(tmp[1]);  
    t1 = numpy.reshape(tmp[2:2+nd], (nd, 1), order='F');
    gd.append(t1)
    if (n==(2+nd+3*nd*natom)):
      t2 = numpy.reshape(tmp[(2+nd):], (3, natom, nd), order='F');
    else:
      t2 = [];
    gdd.append(t2)

  return gd, gdd, na  
    
def fitpod(lammps, podfile, datafile):

  # make a fit podfile file
  text = "units metal";
  text = text + "\n" + "fitpod " + podfile + " " + datafile;
  fid = open("fit.pod", "w")
  fid.write(text)
  fid.close()

  # run LAMMPS to calculate POD descriptors
  runstr = lammps + " -in fit.pod";  
  os.system(runstr);

