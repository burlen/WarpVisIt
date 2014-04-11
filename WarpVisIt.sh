#!/bin/bash
NUM_PROCS=
MPI_EXEC=mpiexec
# get the command line arguments
ARGS=()
for i in "$@"
do
    case $i in
        -n=*)
            NUM_PROCS=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            ;;

        -np=*)
            NUM_PROCS=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            ;;

        --mpiexec=*)
            MPI_EXEC=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            echo "MPI_EXEC=$MPI_EXEC"
            ;;

        --warpvisit-install=*)
            WARPVISIT_INSTALL=`echo $i | sed 's/[-a-zA-Z0-9]*=//'`
            if [ ! -e "$WARPVISIT_INSTALL" ]
            then
                echo "ERROR: $i not found."
                exit
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


        --gdb)
            GDB=1
            ;;

        --massif)
            MASSIF=1
            ;;

        --memcheck)
            MEMCHECK=1
            ;;

        --cli)
            CLI=1
            ;;

        --help)
            HELP=1
            ;;

        *)
            ARGS+=("$i")
            ;;

    esac
done

if [[ -e "$WARPVISIT_SIM2_FILE" ]]
then
    echo "WARNING: Found sim2 file named $WARPVISIT_SIM2_FILE. deleting it."
    rm $WARPVISIT_SIM2_FILE
fi

if [[ -z "$WARPVISIT_INSTALL" ]]
then
    WARPVISIT_INSTALL=`pwd`
fi

if [[ -z "$VISIT_INSTALL" ]]
then
    VISIT_INSTALL=`which visit`
    if [[ -n "$VISIT_INSTALL" ]]
    then
        VISIT_INSTALL=`dirname $VISIT_INSTALL`
    fi
fi
if [[ -n "$VISIT_INSTALL" ]]
then
    source ${WARPVISIT_INSTALL}/WarpVisItEnv.sh $VISIT_INSTALL
fi

if [[ -n "$HELP" ]]
then
    python ${WARPVISIT_INSTALL}/WarpVisItCLIMain.py --help
    python ${WARPVISIT_INSTALL}/WarpVisItEngineMain.py --help
    echo "WarpVisIt.sh"
    echo "-n or -np : : number of MPI processes to launch"
    echo "--visit-install : : path to visit install, configures the environment for you"
    echo "--mpiexec : : mpiexec command and extra args"
    echo "--gdb : : run the program in gdb"
    echo "--cli : : start a VisIt CLI process"
    echo "--massif : : profile using massif"
    echo "--memcheck : : profile using memcheck"
    echo
    exit
fi


# start the run
if [[ -n "$CLI" ]]
then
    echo -n "starting cli..."
    export WARPVISIT_CLI=1
    python ${WARPVISIT_INSTALL}/WarpVisItCLIMain.py "${ARGS[@]}" &
    CLI_PID=$!
    echo $CLI_PID
    trap "echo stopping cli $CLI_PID...; kill $CLI_PID; exit;" SIGHUP SIGINT SIGTERM EXIT
    unset WARPVISIT_CLI
fi
if [[ -z "$NUM_PROCS" ]]
then
    echo "starting serial run..."
    if [[ $GDB -eq 1 ]]
    then
        echo "starting gdb..."
        echo "run WarpVisItEngineMain.py ${ARGS[@]}"
        gdb python
    elif [[ $MASSIF -eq 1 ]]
    then
        echo "starting with massif..."
        valgrind --tool=massif python ${WARPVISIT_INSTALL}/WarpVisItEngineMain.py "${ARGS[@]}"

    elif [[ $MEMCHECK -eq 1 ]]
    then
        echo "starting with memcheck..."
        valgrind --tool=memcheck --leak-check=full --show-reachable=yes --log-file=memcheck.log python ${WARPVISIT_INSTALL}/WarpVisItEngineMain.py "${ARGS[@]}"

    else
        python ${WARPVISIT_INSTALL}/WarpVisItEngineMain.py "${ARGS[@]}"
    fi
else
    echo "starting $NUM_PROCS way parallel run..."
    $MPI_EXEC -n $NUM_PROCS pyMPI ${WARPVISIT_INSTALL}/WarpVisItEngineMain.py "${ARGS[@]}"
fi
