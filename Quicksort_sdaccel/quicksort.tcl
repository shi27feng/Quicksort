# SDAccel command script

# Define a solution name
create_solution -name solution_fpga -dir . -force

# Define the target platform of the application
add_device -vbnv xilinx:adm-pcie-7v3:1ddr:2.0

# Host source files
add_files "host_quicksort.cpp"

# Kernel definition
create_kernel lqsort_kernel -type clc
add_files -kernel [get_kernels lqsort_kernel] "QuicksortKernels.cl"

create_kernel  gqsort_kernel -type clc
add_files -kernel [get_kernels gqsort_kernel] "QuicksortKernels.cl"


# Define binary containers
create_opencl_binary Quicksort
set_property region "OCL_REGION_0" [get_opencl_binary Quicksort]

create_compute_unit -opencl_binary [get_opencl_binary Quicksort] -kernel [get_kernels lqsort_kernel] -name K1

create_compute_unit -opencl_binary [get_opencl_binary Quicksort] -kernel [get_kernels gqsort_kernel] -name K2

# Compile the design for CPU based emulation
compile_emulation -flow cpu -opencl_binary [get_opencl_binary Quicksort]

# Generate the system estimate report
report_estimate

# Run the design in CPU emulation mode
run_emulation -flow cpu -args "Quicksort.xclbin" 

# Build the application for hardware
#build_system

# Package the results for the card
#package_system
