#!/bin/bash

# ENTER R SCRIPT NAME
r_script=

# ENTER PATH TO YOUR R SCRIPT
r_script_path=

# ENTER PATH TO YOUR INPUT FILES
input_dir=

# ENTER PATH TO YOUR cancer_type_info.txt requires BY CNCDriver LIBRARY
cancer_type_file=

# ENTER PATH TO YOUR CNCDriver DATA FILES
cncdriver_data=

# ENTER PATH TO YOUR OUTPUT DIRECTORY
output_dir=

echo "Running CNCDriver script within docker container"

# run docker for the r_script defined above
docker run \
    -v $r_script_path$r_script:/home/$r_script \
    -v $input_dir:/home/input/ \
    -v $cancer_type_file:/home/cancer_type_info.txt \
    -v $cncdriver_data:/home/CNCDriver_data/ \
    -v $output_dir:/home/output/ \
     cnc_driver_docker Rscript /home/$r_script
        
echo "Done running CNCDriver script"
