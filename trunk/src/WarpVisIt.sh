#!/bin/bash

NUM_PROCS=1
MPI_EXEC=mpiexec

# get the command line arguments
for i in "$@"
do
    case $i in
        -np=*)
            NUM_PROCS=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            ;;

        --mpiexec=*)
            MPI_EXEC=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            echo "MPI_EXEC=$MPI_EXEC"
            ;;

        --warp-script=*)
            WARP_SCRIPT=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            if [ ! -e "$WARP_SCRIPT" ]
            then
                echo "ERROR: $i not found."
                exit
            else
                export WARP_SCRIPT
                echo "WARP_SCRIPT=$WARP_SCRIPT"
            fi
            ;;

        --visit-install=*)
            VISIT_INSTALL=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            if [ ! -d "$VISIT_INSTALL" ]
            then
                echo "ERROR: $i not found."
                exit
            fi
            ;;

        --visit-script=*)
            VISIT_SCRIPT=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            if [ ! -e "$VISIT_SCRIPT" ]
            then
                echo "ERROR: $i not found."
                exit
            else
                export VISIT_SCRIPT
                echo "VISIT_SCRIPT=$VISIT_SCRIPT"
            fi
            ;;

        --sim2-file=*)
            SIM2_FILE=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            export SIM2_FILE
            echo "SIM2_FILE=$SIM2_FILE"
            ;;

        --help)
            echo
            echo "Usage..."
            echo "$0 [-np=X] --warp-script=X --visit-install=X [--sim2-file=X] [--visit-script=X]"
            echo
            echo "    --warp-script   : Warp python code to configure the run"
            echo "    --visit-install : path to the VisIt install directory"
            echo "    --sim2-file     : where the run should place the sim2 file (optional)"
            echo "    --visit-script  : VisIt python code to configure VisIt's rendering pipeline (optional)"
            echo
            exit
            ;;

        *)
            echo
            echo "ERROR: Invalid option $i"
            echo
            exit
            ;;
    esac
done

if [[ -z "$WARP_SCRIPT" ]]
then
    echo
    echo "ERROR: --warp-script=... is a manditory argument."
    echo "       use --help to show usage info."
    echo
    exit
fi

if [[ -n "$VISIT_INSTALL" ]]
then
    # get VisIt's environment
    source WarpVisItEnv.sh $VISIT_INSTALL
fi

# start the run
if [[ $NUM_PROCS < 2 ]]
then
    echo "starting serial run..."
    LD_PRELOAD=$VISIT_MESA_LIB python WarpVisItMain.py
else
    echo "starting $NUM_PROCS way parallel run..."
    $MPI_EXEC -np $NUM_PROCS ./WarpVisItMPIMain.sh
fi
