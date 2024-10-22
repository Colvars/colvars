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

4. (Optional) Push the container to the GitHub repository's "Packages"
   section; currently we only have one package, called `devel-containers`,
   which comes in different versions depending on the OS and software
   contained.  For example, to update the most commonly used container
   version:
```
apptainer push CentOS9-devel oras://ghcr.io/colvars/devel-containers:CentOS9-devel
```
(note that the above requires having set up an access token for apptainer)
