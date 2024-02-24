#!/bin/bash

# Remote machine details
REMOTE_USERNAME="lsc"
REMOTE_HOST="10.24.6.132"
REMOTE_DIR="Desktop/ptsn4/harald"

# Local file paths
#LOCAL_PATH="~/Desktop/surface_gf/ptsn4/harald_data"

# Function to transfer files and compile and execute Fortran code remotely
run_fortran_remote() {
    # Modify the value of kkm in ptsn4_bulk_fs.f90
    powershell -Command "(Get-Content ptsn4_bulk_fs.f90) -replace '^ *INTEGER *:: *kkm *= *.*$', '  INTEGER :: kkm = $1' | Set-Content ptsn4_bulk_fs.f90"

    # Transfer kmesh.f90 and modified ptsn4_bulk_fs.f90 to remote machine
    scp kmesh.f90 $REMOTE_USERNAME@$REMOTE_HOST:$REMOTE_DIR/kmesh.f90
    scp ptsn4_bulk_fs.f90 $REMOTE_USERNAME@$REMOTE_HOST:$REMOTE_DIR/ptsn4_bulk_fs.f90
    
    echo "File transfer done."

    # Compile and execute kmesh.f90
    ssh $REMOTE_USERNAME@$REMOTE_HOST "cd $REMOTE_DIR && gfortran -o kmesh kmesh.f90 && ./kmesh" 
    
    # Check if compilation and execution of kmesh.f90 were successful
    if [ $? -eq 0 ]; then
        echo "kmesh compilation and execution successful."
    else
        echo "kmesh compilation or execution failed. Exiting."
        exit 1
    fi

    # Compile ptsn4_bulk_fs.f90
    ssh $REMOTE_USERNAME@$REMOTE_HOST "cd $REMOTE_DIR && gfortran -o pt ptsn4_bulk_fs.f90 -llapack"
    
    # Check if compilation was successful
    if [ $? -eq 0 ]; then
        echo "Fortran compilation successful. Running Fortran program remotely..."
        # Execute Fortran program remotely
        ssh $REMOTE_USERNAME@$REMOTE_HOST "cd $REMOTE_DIR && ./pt"
        # Transfer output file back to local machine
        scp $REMOTE_USERNAME@$REMOTE_HOST:$REMOTE_DIR/fermi_surface_0kx.dat .
        echo "Output file transferred successfully."
    else
        echo "Fortran compilation failed. Exiting."
        exit 1
    fi
}

# Main script
echo "Transferring files and running Fortran code on remote machine..."
# Provide the desired value for kkm here (replace 50 with your desired value)
run_fortran_remote 20

