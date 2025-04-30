## Colvars development tools and scripts

### Automatically updating the version strings

The functions defined in `version_functions.sh` allow to keep the
COLVARS_VERSION macro and related ones up to date.  From the root folder:
```
source devel-tools/version_functions.sh 
update_all_versions <last_commit_before_change>
```
The argument of the `update_all_versions` function is the hash of a commit
that precedes the changes added since the last time that the macros were
updated.  A Git commit is auto-generated, which will also trigger a doc
rebuild when pushed to GitHub.

Note that the new doc will be pushed to the website under a folder named
like the current branch; this means that if the commit described above is
not the last commit in a pull request, direct-push access to master will be
required for updating the doc.


### Git commit hooks

Type `make` in this folder to install a pre-commit hook, which will remove
all trailing whitespace from the staged changes.
