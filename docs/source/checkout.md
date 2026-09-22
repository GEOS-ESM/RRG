# Checking out the model
## Setting up GitHub credentials

Set up a GitHub account. You can check out read-only as a general user, but must be part of the GEOS-ESM group to checkout with commit privileges. Since the GEOS admins are unlikely to know you if you're a new user, you need to ask an existing user to introduce you.

## Checking out GEOS
Much of this is already documented in the
* [GEOSgcm README](https://github.com/GEOS-ESM/GEOSgcm?tab=readme-ov-file)

Only the essential steps are documented below.

### Setting up the environment
First, you need to have the correct modules to check out the code. If you are expert user, make sure you have `git` and `mepo` in your path. If you're not, execute the following to get them (as well as anything else you might need).
```
module purge
module use -a /discover/swdev/gmao_SIteam/modulefiles-SLES15
module load GEOSenv
```
Note the above will only work on `SLES15` nodes. For `SLES12` nodes, replace `SLES15` above with `SLES12`. Currently all login nodes are `SLES12`, all Milan compute nodes are `SLES15`, and Cascade Lake compute nodes are available on both `SLES12` and `SLES15`.

Other notes: There is currently an issue with `cmake` identifying the correct `Python 3` install. If a `mepo` or `parallel_build.csh` run dies, especially near `f2py`, this may be the problem. Working on a fix ...

### Getting the code
Then, decide where you want to set up the model code. While you _can_ check it out in your home directory, we recommend checking it out on a scratch space because it can get pretty big, especially when you are testing multiple different model versions. On NCCS, you typically want to use `$NOBACKUP`. Somewhere in that folder, check out the model with
```
git clone -b v11.10.1 git@github.com:GEOS-ESM/GEOSgcm.git GEOSgcm-v11.10.1
```
Note that `11.10.1` simply happens to be the latest released tag at the time this is being written. There is nothing sacred about that tag. If you want a later or earlier tag, you can find all release versions [here](https://github.com/GEOS-ESM/GEOSgcm/releases).

The model consists of code in several sub-repositories, by default none of which are checked out. There is a file called `components.yaml` in the source tree you just checked out, which contains the tags for each repository that _will be_ checked out.
