__precompile__()

module POD

export calculatedescriptors, fitpod

function string2cmd(str::String)

  ind = findall(" ", str);
  n = length(ind);
  cmdstr = Array{String,1}(undef,n+1);
  cmdstr[1] = str[1:(ind[1][1]-1)];
  for i = 2:n
      cmdstr[i] = str[(ind[i-1][1]+1):(ind[i][1]-1)];
  end
  i = n+1;
  cmdstr[i] = str[(ind[i-1][1]+1):end];

  return Cmd(cmdstr)

end

function readbin(filename)
  a = reinterpret(Float64,read(filename));
  return a    
end

function getdescriptors(filename)

  # get data from the binary file
  tmp = readbin(filename);
  n = length(tmp);

  # global POD descriptors
  natom = Int64(tmp[1]);
  nd = Int64(tmp[2]);
  gd = reshape(tmp[3:2+nd], (nd, 1));

  # derivatives of the global POD descriptors wrt atom positions
  if (n==(2+nd+3*nd*natom))
    gdd = reshape(tmp[(2+nd+1):end], (3, natom, nd));
  else
    gdd = [];
  end

  return gd, gdd, natom
end

function calculatedescriptors(lammps, pod, datapath)

  # make a compute pod file
  text = "file_format extxyz";
  text = text * "\n" * "file_extension xyz";
  text = text * "\n" * "path_to_training_data_set " * datapath;
  text = text * "\n" * "compute_pod_descriptors 2";
  open("compute.pod","w") do io
    write(io,text)
  end

  # make a fit pod file
  text = "units metal";
  text = text * "\n" * "fitpod " * pod * " compute.pod";
  open("fit.pod","w") do io
    write(io,text)
  end

  # run LAMMPS to calculate POD descriptors
  runstr = lammps * " -in fit.pod";  
  run(string2cmd(runstr), wait=true);

  # determine the number of configurations in the dataset
  k = 0;
  while (true)
    filename = datapath * "/descriptors_config" * string(k+1) * ".bin";
    if isfile(filename)
      k = k + 1;
    else
      break;
    end
  end
  num_configs = k;

  # get POD descriptors from the binary files
  gd = Array{Any}(undef,num_configs) 
  gdd = Array{Any}(undef,num_configs)  
  na = zeros(Int64, num_configs, 1);
  for i = 1:num_configs
    filename = datapath * "/descriptors_config" * string(i) * ".bin";
    gd[i], gdd[i], na[i] = getdescriptors(filename); 
  end

  return gd, gdd, na
end

function fitpod(lammps, pod, data)

  # make a fit pod file
  text = "units metal";
  text = text * "\n" * "fitpod " * pod * " " * data;
  open("fit.pod","w") do io
    write(io,text)
  end

  # run LAMMPS to calculate POD descriptors
  runstr = lammps * " -in fit.pod";  
  run(string2cmd(runstr), wait=true);
end

end
