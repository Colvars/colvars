## Colvars development tools and scripts


### Automatically updating the version strings (`master` branch only)

The functions defined in `version_functions.sh` allow to keep the `COLVARS_VERSION` macro and
related ones up to date.  From the root folder:
```
source devel-tools/version_functions.sh 
update_all_versions <last_commit_before_change>
```
The argument of the `update_all_versions` function is the hash of a commit that *preceded* the
changes added since the last time the macros were updated.  In most cases, this commit is the
`master` branch itself and this task is run by a CI job when a pull request is merged.

The Git commit message describing the version update also contains the tag that triggers a doc
rebuild (see next section).


### Rebuilding the documentation

If a Git commit message contains `[update-doc]`, it will trigger a doc rebuild by a CI job when
when the commit is pushed to GitHub.

The CI job will place the documentation in a branch with the same name as the current branch:
this makes it useful to e.g. automatically build the documentation for the sake of reviewing a
pull request.


### Making updates to other packages

Colvars is primarily distributed as part of other packages.  To make an update, one should first
merge all recent commits from the `master` branch into a branch dedicated to the upcoming release
of that package.  There should be no additional commits in the latter branch compared to
`master`, allowing a fast-forward update.  After that, the `update-colvars-code.sh` script could
be used to patch a worktree of the target package.

For example, the following process is typically used to update LAMMPS:
```
git checkout lammps-develop
git rebase master
./update-colvars-code.sh /path/to/LAMMPS/worktree
cd /path/to/LAMMPS/worktree
git checkout -b colvars-update
git add .
git commit
```

For the purpose of describing the new commit, the `print_pr_list.py` script can be used to
generate a summary of the latest changes in Colvars.  To use the script, take a note of the date
corresponding to the Colvars version currently included in that package, and use that date to
limit the search using the `--since` option.  (This script uses internally the [GitHub
CLI](https://cli.github.com/).)

For example, running the following command from a Colvars worktree generates a summary of pull
requests affecting GROMACS that were merged since 2024-11-14:
```
./print_pr_list.py --backend GROMACS --since 2024-11-14
```
The `--backend` option will automatically filter out PRs that are tagged as specific to another
backend, or PRs tagged with the `maintenance` tag (relevant only to the Colvars repo).

**Note:** we are not always consistent in tagging PRs :-) It is thus a good idea to review the
generated summary, and remove irrelevant entries from it before using it in a commit message to
other packages.  Additionally, the date check may end up also including PRs that had been already
included, but were merged on the same day as the argument to `--since`.


### Managing bugfix updates

GROMACS, LAMMPS and NAMD periodically issue updates to their "stable" releases.  There should be
branches in this repo tracking those releases for each engine.  Because Colvars bugfixes are
generally included first in the `master` branch, they also need to be cherry-picked into those
branches.

For example to find out which PRs fix bugs relevant to NAMD, use the following:
```
./print_pr_list.py --backend NAMD --since <last update> --label bugfix
```
to find all bugfix PRs.  From that list, manually cherry-pick commits from those PRs into the
relevant `namd-3.x` branch (this step may be automated in the future).

Next, manually edit the file `src/colvars_version.h` and manually increment the value of
`COLVARS_PATCH_VERSION`, using e.g. the patch release number of the target package.
The value of `COLVARS_VERSION` should remain unchanged, because no features have been added.

From then on, follow the update procedure is then the same as in the previous section, and submit
a change request to the other package into their stable release branch.
