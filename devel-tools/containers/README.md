1. Install Apptainer (formerly Singularity), either from your Linux
   distribution's package manager or from https://github.com/apptainer/apptainer/releases/

2. Build a container from its definition file using:
```
apptainer build container.sif container.def
```

3. Run a shell inside the container:
```
apptainer run container.sif
```
