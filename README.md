
### Compiling
To build the executable, make sure you have loaded the needed modules with
```bash
$ module load gcc-glibc dealii
```
Then run the following commands:
```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
```

This will create the executable in the `build` directory.

### Execution

#### Standard Execution
To run the executable, navigate to the `build` directory and use the following command:

```bash
$ ./solver <dimension> <mesh_file_name> [parameters_file_name]
```

#### MPI Execution
For parallel execution using MPI, run the following command:

```bash
$ mpirun -np <number_of_processes> ./solver <dimension> <mesh_file_name> [parameters_file_name]

```

- `<number_of_processes>`: Specify the number of processes to run.
- `<dimension>`: Specify the problem dimension (e.g., 2 or 3).
- `<mesh_file_name>`: Provide the name of the mesh file.
- `[parameters_file_name]` (optional): Provide the name of the parameters file if needed.

#### Convergence Test Execution
To run the convergence test, navigate to the `build` directory and use the following commands, replacing # with 1 or 2 for specifying the dimension:

```bash
$ ./convergence#D ../input/test#D_convergence.prm
$ ../scripts/plot-convergence.py convergence.csv
```

### Running on a Cluster
To run the program on a cluster, you can use a sample PBS job script in the `scripts` directory. This script submits a job that requests 1 node with 20 processors, and a wall time of 12 hours.

#### Executing the Program

To run the program with 20 processes, distributing them across the node. The output is redirected to stdout.txt.

```bash
$ mpirun -machinefile mpd.nodes -n 20 -npernode 20 solver <dimension> <mesh_file_name> [parameters_file_name] &> stdout.txt
```

#### Submitting the Job
To submit the job script to the cluster, use the following command:
```bash
qsub job_script.sub
```

### Mesh Generation

#### Prerequisites
- Ensure that Gmsh is installed on your system and is available in your system's PATH.

#### Generate the 3D brain mesh
To generate a 3D brain mesh, run:
```bash
gmsh scripts/convert_stl_to_msh.geo -nopopup
```

#### Generate the ellipse mesh
To generate an ellipse mesh, open the `scripts/ellipse.geo` file and adjust the `refinement` parameter to increase or decrease the level of mesh refinement as needed. Then run:

```bash
gmsh scripts/ellipse.geo -nopopup
```


