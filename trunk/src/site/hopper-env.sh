module swap PrgEnv-pgi PrgEnv-gnu
module load visit/2.5.2
#grep the LD_LIBRARY_PATH, VISITPLUGINDIR and VISIT_MESA_LIB environment variables 
#from VisIt (by executing visit -env)  and set the environment variables accordingly
eval `visit -env | grep LD_LIBRARY_PATH | sed -e 's/^/setenv /' -e 's/=/ /'`
eval `visit -env | grep VISITPLUGINDIR | sed -e 's/^/setenv /' -e 's/=/ /'`
eval `visit -env | grep VISIT_MESA_LIB | sed -e 's/^/setenv /' -e 's/=/ /'`

#~grote/hopper/public/python/bin/python -i ...py
#~grote/hopper/public/python/bin/pyMPI
