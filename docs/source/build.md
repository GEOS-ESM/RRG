# Building the code

## The easy way

```
./parallel_build.csh -mil
```

## The hard way

You must build the code on a compute node of the same architecture as the ones you will be running the model on. For this example we will be building and running on **AMD Milan** nodes.

Run `mepo clone` at the command line to check out all the repositories at the branches/tags in `components.yaml`. Before you do this, it is advisable to run `mepo config set clone.partial blobless` to make checkouts faster. You can read more details [here](https://github.com/GEOS-ESM/GEOSgcm/wiki#recommended-mepo-settings-for-geosgcm).

Then get a terminal on such a compute node with
```
salloc -A s1460 --nodes=1 --constraint=mil --qos=debug -t 60
```
This gets you a terminal on a Milan node under the `debug` queue, which is pretty fast but has a wall clock limit of 1 hour. You could, alternatively, issue this command first thing in the morning without `--qos=debug` and with `-t 480` and get a node for 8 hours. You will need to wait longer to get a node, but once you do, you're set for a day's worth of building and debugging.

Extra details: The `-t 60` command is optional for `--qos=debug` which defaults to and is capped at an hour. The `-A s1460` command is also optional assuming you've run something before. You must, of course, be on the `s1460` compute code use this account. If not, substitute your account charge code. As of 2026-08-28, it seems that Milan nodes are not the default, i.e., you can still land on a Cascade Lake node. So you still need to specify `--constraint=mil`.

Once you get on a compute node, go to the folder where you checked out the source tree, and just to be safe create a clean environment as follows (this assumes you are running a bash/ksh/zsh variant and not C shell):
```
module purge
source @env/g5_modules.sh
```
Now you're ready to build. Since you're already on a compute node, no need to submit a parallel build job. Instead, issue the following commands in order:
```
mkdir build
cd build
cmake .. -DBASEDIR=$BASEDIR/Linux -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install
make -j install
```
This builds the model into `../install`, specifically the GEOS GCM executable is `../install/bin/GEOSgcm.x`. Note that this builds the model with all optimizations turned on. Occasionally you will need to do a *debug* build that checks all array bounds, overflows, etc. Then you substitute `-DCMAKE_BUILD_TYPE=Release` above with `-DCMAKE_BUILD_TYPE=Debug`.

**Important:** When you run `gcm_setup` to set up a new run, this executable is copied over to the run directory. As a result, if you want to fix something in code and recompile, the changes will not be seen in your run unless you copy over the executable again. Therefore, I often symlink `install/bin/GEOSgcm.x` from my run directory.