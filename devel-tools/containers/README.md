## Prerequisites

Install Docker (https://docs.docker.com/engine/install/) or Podman (https://podman.io/); with Podman, replace `docker` with `podman` in all the commands below.  Both packages can create containers from a Dockerfile format.

Podman in particular is best suited for rootless, daemonless operations; if your home is on a network folder, you may want to add these to the file `~/.config/containers/storage.conf`:
```
[storage]
driver = "overlay"
graphroot = "<local_folder>/storage"
runroot = "<local_folder>/runroot"

[storage.options.overlay]
mount_program = "/usr/bin/fuse-overlayfs"
```

For portability, we are no longer using the Apptainer-specific “`.def`” format: however, if Apptainer is the only container software available you can still use it to run the containers produced by Docker or Podman (see below).


## Building containers

Build a container image from its Dockerfile using Docker, e.g.:
```
docker build --tag colvars/devel-containers:CentOS9-devel-x86_64 --file CentOS9-devel-x86_64.Dockerfile .
```
or:
```
podman build --format docker --tag colvars/devel-containers:CentOS9-devel-x86_64 --file CentOS9-devel-x86_64.Dockerfile .
```


## Testing containers

You can open an interactive shell inside a container, optionally mounting the current directory as in the example below using Podman:
```
podman run --rm --interactive --tty -v "$PWD":"$PWD" -w "$PWD" ghcr.io/colvars/devel-containers:CentOS9-devel-x86_64
```
or Apptainer:
```
apptainer run docker://ghcr.io/colvars/devel-containers:CentOS9-devel-x86_64
```


## Push the container to GHCR (optional)

To use the container in workflows, you will need to push it to the GitHub repository's "Packages" section as a Docker image that GitHub Actions can use directly.  All containers used to perform build tests for Colvars are grouped under the `devel-containers` package, with different tags depending on the OS and software contained.

For example, to update the most commonly used container version:
```
printf '%s' "$GHCR_TOKEN" | docker login ghcr.io --username "$GITHUB_USER" --password-stdin
docker tag colvars/devel-containers:CentOS9-devel-x86_64 ghcr.io/colvars/devel-containers:CentOS9-devel-x86_64
docker push ghcr.io/colvars/devel-containers:CentOS9-devel-x86_64
```

Or, with Podman (build with `--format docker` as shown above):
```
printf '%s' "$GHCR_TOKEN" | podman login ghcr.io --username "$GITHUB_USER" --password-stdin
podman tag colvars/devel-containers:CentOS9-devel-x86_64 ghcr.io/colvars/devel-containers:CentOS9-devel-x86_64
podman push --format v2s2 ghcr.io/colvars/devel-containers:CentOS9-devel-x86_64
```

Set `GITHUB_USER` to your GitHub username and `GHCR_TOKEN` to a personal access token (classic) with `write:packages` permission and access to the package.

GitHub Actions can reference the image as
```
container:
  image: ghcr.io/colvars/devel-containers:CentOS9-devel-x86_64
```
