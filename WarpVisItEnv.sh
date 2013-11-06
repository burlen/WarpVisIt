#!/bin/bash
# this script uses VisIt to set the environment
# variables that it needs. it must be given the
# path to a VisIt install. On NERSC systems just
# use the module.
if [[ $# != 1 ]] || [[ ! -e $1/bin/visit ]]
then
  echo
  echo "ERROR:"
  echo "  usage"
  echo "    setupenv.sh /path/to/visit-install"
  echo
  exit
fi

# all visit specific stuff can be used directly
VISIT_INSTALL=$1
$VISIT_INSTALL/bin/visit -env | grep -v LD_LIB | sed -e 's/^/export /' > .tmpVisItEnv

# handle this separately because visit left
# alone will destroy the LD_LIBRARY path
ld_vars=`$VISIT_INSTALL/bin/visit -env | grep LD_LIB`
for var in $ld_vars
do
  echo "export $var:$LD_LIBRARY_PATH" >> .tmpVisItEnv
done

# libsim needs yet another variable , with the same info
VISIT=`$VISIT_INSTALL/bin/visit -env | grep VISITARCHHOME | sed -e 's/VISITARCHHOME=//'`
echo "export VISIT=$VISIT" >> .tmpVisItEnv

# the viewer neeeds PATH set
echo "export PATH=$VISIT_INSTALL/bin:$PATH" >> .tmpVisItEnv

# dump for verification
cat .tmpVisItEnv

# bring them in
source .tmpVisItEnv

# clean up
rm .tmpVisItEnv
