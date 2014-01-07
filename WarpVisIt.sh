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
            WARPVISIT_WARP_SCRIPT=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            if [ ! -e "$WARPVISIT_WARP_SCRIPT" ]
            then
                echo "ERROR: $i not found."
                exit
            else
                export WARPVISIT_WARP_SCRIPT
                echo "WARPVISIT_WARP_SCRIPT=$WARPVISIT_WARP_SCRIPT"
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

        --interactive)
            echo "WARPVISIT_INTERACTIVE=1"
            export WARPVISIT_INTERACTIVE=1
            ;;

        --script-dir=*)
            WARPVISIT_SCRIPT_DIR=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            if [ ! -e "$WARPVISIT_SCRIPT_DIR" ]
            then
                echo "ERROR: $i not found."
                exit
            else
                export WARPVISIT_SCRIPT_DIR
                echo "VISIT_SCRIPT=$VISIT_SCRIPT"
            fi
            ;;

        --sim2-file=*)
            WARPVISIT_SIM2_FILE=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            export WARPVISIT_SIM2_FILE
            echo "WARPVISIT_SIM2_FILE=$WARPVISIT_SIM2_FILE"
            ;;

        --gdb)
            GDB=1
            ;;

        --help)
            echo
            echo "Usage..."
            echo "$0 [-np=X] --warp-script=X --visit-install=X [--sim2-file=X] [--visit-script=X]"
            echo
            echo "    --warp-script   : Warp python code to configure the run"
            echo "    --visit-install : path to the VisIt install directory"
            echo "    --sim2-file     : where the run should place the sim2 file (optional)"
            echo "    --script-dir    : where to find user supplied scripts (optional)"
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

# validate and set defaults
if [[ -z "$WARPVISIT_WARP_SCRIPT" ]]
then
    echo
    echo "ERROR: --warp-script=... is a manditory argument."
    echo "       use --help to show usage info."
    echo
    exit
fi

if [[ -n "$VISIT_INSTALL" ]]
then
    source WarpVisItEnv.sh $VISIT_INSTALL
fi

if [[ -z "$WARPVISIT_SCRIPT_DIR" ]]
then
    export WARPVISIT_SCRIPT_DIR=`pwd`
fi

# start the run
if [[ $NUM_PROCS < 2 ]]
then
    echo "starting serial run..."
    if [[ $GDB -eq 1 ]]
    then
        echo "starting gdb..."
        #echo "set environment LD_PRELOAD=$VISIT_MESA_LIB"
        echo "run WarpVisItMain.py"
        gdb python
    else
        #LD_PRELOAD=$VISIT_MESA_LIB python WarpVisItMain.py
        python WarpVisItMain.py
    fi
else
    echo "starting $NUM_PROCS way parallel run..."
    #$MPI_EXEC -np $NUM_PROCS ./WarpVisItMPIMain.sh
    $MPI_EXEC -np $NUM_PROCS pyMPI WarpVisItMain.py
fi
