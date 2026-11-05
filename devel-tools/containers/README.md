1. Install Docker (https://docs.docker.com/engine/install/) or Podman
   (https://podman.io/); with Podman, replace `docker` with `podman` in
   the commands below.

2. Build a container image from its Dockerfile, e.g.:
```
docker build -t colvars/devel-containers:CentOS9-devel -f CentOS9-devel.Dockerfile .
```
   (building the `arm64` image on a non-ARM machine requires QEMU emulation
   to be set up, see https://docs.docker.com/build/building/multi-platform/)

3. Run a shell inside the container, mounting the current directory:
```
docker run --rm -it -v "$PWD":"$PWD" -w "$PWD" colvars/devel-containers:CentOS9-devel
```

4. (Optional) Push the container to the GitHub repository's "Packages"
   section; currently we only have one package, called `devel-containers`,
   which comes in different versions depending on the OS and software
   contained.  The GitHub Actions workflows (via `devel-tools/get_container`)
   download these as Apptainer images, so convert the Docker image first.
   For example, to update the most commonly used container version:
```
apptainer build CentOS9-devel.sif docker-daemon://colvars/devel-containers:CentOS9-devel
apptainer push CentOS9-devel.sif oras://ghcr.io/colvars/devel-containers:CentOS9-devel
```
(note that the above requires having set up an access token for apptainer)
