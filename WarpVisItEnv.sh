#!/bin/bash
# this script uses VisIt to set the environment
# variables that it needs. it must be given the
# path to a VisIt install.

for i in "$@"
do
    case $i in
        --visit-install=*)
            VISIT_INSTALL=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            if [ ! -d "$VISIT_INSTALL" ]
            then
                echo "ERROR: VisIt install at ($i) doesn't exist."
                exit -1
            fi
            ;;

        --verbose)
            VERBOSE=1
            ;;
    esac
done

if [[ -z "$VISIT_INSTALL" ]]
then
    VISIT_INSTALL=`which visit`
    if [[ -n "$VISIT_INSTALL" ]]
    then
        VISIT_INSTALL=`dirname $VISIT_INSTALL`
    fi
fi

if [[ -z "$VISIT_INSTALL" ]]
then
    echo "ERROR: WarpVisItEnv.sh could not locate the visit install."
    echo
    echo "  usage"
    echo "    WarpVisItEnv.sh --visit-install=/path/to/visit-install [--verbose]"
    echo
    exit -1
fi

# all visit specific stuff can be used directly
$VISIT_INSTALL/bin/visit -env -engine | grep -v LD_LIB | sed -e 's/^/export /' > .tmpVisItEnv

# handled LD_LIBRARY_PATH separately because visit left
# alone will destroy it
ld_vars=`$VISIT_INSTALL/bin/visit -env -engine | grep LD_LIB`
for var in $ld_vars
do
  echo "export $var:$LD_LIBRARY_PATH" >> .tmpVisItEnv
done

# libsim needs yet another variable, with the same info
VISIT=`$VISIT_INSTALL/bin/visit -env -engine | grep VISITARCHHOME | sed -e 's/VISITARCHHOME=//'`
echo "export VISIT=$VISIT" >> .tmpVisItEnv

# the viewer neeeds PATH set
echo "export PATH=$VISIT_INSTALL/bin:$PATH" >> .tmpVisItEnv

# dump what we set for verification
if [[ "$VERBOSE" == "1" ]]
then
    cat .tmpVisItEnv
fi

# bring them in
source .tmpVisItEnv

# clean up
rm .tmpVisItEnv
