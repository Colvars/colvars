# Create patched Colvars GROMACS versions

Steps to create a release archive of a GROMACS version with Colvars support.
This is done on the GROMACS fork of Colvars organization: https://github.com/Colvars/gromacs

All commands are done within the GROMACS source tree folder.

If a Colvars release branch already exists with a Colvars commit, go to step 3.

## 1. Create a derived branch from a release one

```
git checkout -b release-XXXX-colvars release-XXX
```

## 2. Apply, check and commit Colvars patch

This step implies that the Colvars sources support the desired GROMACS version.

```
# Applying patch
/path/to/colvars/update-colvars-code.sh .
```

Make sure the patch works by compiling GROMACS and eventually, run the Colvars GROMACS tests.

Once verified, commit the patch with its version name within the message:
```
git add <list of files>
git commit -m "Colvars patch versions XXXX-XX-XX"
```

## 3. Create a Colvars-GROMACS tag from a GROMACS minor version tag

Go to the GROMACS minor version to be patched:

```
git checkout tags/vYYYY.Y
```

Cherry pick the Colvars commit from step #2 to be applied on this tag:

```
git cherry-pick <commit SHA1>
```

Make sure the commit works by compiling GROMACS and running Colvars tests.

Create a new tag:
```
git tag -a "vYYYY.Y-colvars" -m "Colvars patch version XXXX-XX-XX for Gromacs YYYY.Y"
```

Return to master and push the tag:
```
git checkout -
git push origin vYYYY.Y-colvars
```
