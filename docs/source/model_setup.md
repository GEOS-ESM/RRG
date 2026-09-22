# Setting up a run

## The easy way

Use `gcm_setup` to clone the run in `/discover/nobackup/bweir/GEOS/runs/carbon-ng_ana`. You can get restarts like this:
```
tar xf /discover/nobackup/projects/gmao/geos_carb_dev/bweir/runs/carbon-ng_ana/restarts/restarts.e20141023_21z.tar
~bweir/bin/striprst.sh
```
The run:
```
sbatch ./gcm_run.j
```
Although it'd probably be better to get an interactive session and try it, e.g.,
```
salloc --time=10:00:00 --constraint=mil --ntasks=1200 --ntasks-per-node=120
```

## The hard way

### Use `gcm_setup` to install the model somewhere
You don't need to be on a compute node for this. On a login node, go into `install/bin` and execute `./gcm_setup`.

* `Experiment ID` is any name you want to give the run. It's a good idea to include the model version and something about which tracers you are running in a short name. Mostly something that you will remember. If you call it (say) `Apple` it's pretty much guaranteed that you won't remember what it is for two years down the line. I'm calling mine `GCM-11.5.2-methane-c180`.
* `Experiment Description` is a short description to help you remember.
* `CLONE` is the ability of a model to copy over someone else's run folder. This is a very useful ability, but for now let's choose `NO`.
* `Atmospheric Horizontal Resolution` depends on what you want to run. I'm choosing `c90`, which is perfectly fine for model development, or for decadal runs just to test out emissions. For science runs, you might want to choose `c180`.
* `Vertical Resolution` of 72 layers is fine, but if you chose `c180` or finer above, your default might be `181`.
* Default `Microphysics` of `BACM_1M` is fine
* Default `TRUE` for `Hydrostatic Atmosphere` is fine
* Use `IOSERVER` if you're running `c180` or higher
* Default processor type of `mil` is fine
* Default `NO` to `COUPLED Ocean/Sea-Ice Model` is fine
* Choose `CS` (cubed sphere) for `Data_Ocean Horizontal Resolution`
* Default choice `Icarus-NLv3` for land surface boundary conditions is fine
* Default choice `Catchment` for land surface model is fine
* Accept default choice to run GOCART with `Actual` aerosols
* Choose to use `OPS` emission files for GOCART, because the `AMIP` emission files do not exist for recent years
* For `c180`, a `HEARTBEAT_DT` of `450` is fine
* Don't worry about the `HISTORY template`, you are going to change the history file anyway
* The `HOME Directory` is where the run folder will be created. Just make sure it's created somewhere inside one of the large partitions you have allocation on. E.g., `/discover/nobackup/${USER}` or `/discover/nobackup/projects/gmao/geos_carb/${USER}`. If you're starting out with the GCM, try the first one. Eventually when you do production runs you will switch to the latter.
* In theory `EXPERIMENT Directory` can be different from `HOME Directory`, but no one has ever tried it. Either set it to be the same, or try at your own risk and don't expect any sympathy if you break something.
* The `Build directory` should already be correct
* Our `GROUP ID` is `s1460`

### `TMPDIR` issue on discover

Every so often `gcm_setup` will fail with errors like
```
/tmp/tmp.sVmzWQGKy5: Permission denied.
/bin/mv: cannot stat '/discover/nobackup/projects/gmao/geos_carb/sbasu1/runs/GCM/test_restarts/AGCM.rc.tmpl': No such file or directory
/tmp/tmp.VXxsGkzxBY: Permission denied.
/bin/mv: cannot stat '/discover/nobackup/projects/gmao/geos_carb/sbasu1/runs/GCM/test_restarts/AGCM.rc.tmpl': No such file or directory
/tmp/tmp.kAg9MZYON8: Permission denied.
/bin/mv: cannot stat '/discover/nobackup/projects/gmao/geos_carb/sbasu1/runs/GCM/test_restarts/AGCM.rc.tmpl': No such file or directory
/tmp/tmp.Xced0YAEi2: Permission denied.
/bin/mv: cannot stat '/discover/nobackup/projects/gmao/geos_carb/sbasu1/runs/GCM/test_restarts/AGCM.rc.tmpl': No such file or directory
/tmp/tmp.igPE8a3dYw: Permission denied.
/bin/mv: cannot stat '/discover/nobackup/projects/gmao/geos_carb/sbasu1/runs/GCM/test_restarts/AGCM.rc.tmpl': No such file or directory
/tmp/tmp.OL9WF9FmLj: Permission denied.
/bin/mv: cannot stat '/discover/nobackup/projects/gmao/geos_carb/sbasu1/runs/GCM/test_restarts/AGCM.rc.tmpl': No such file or directory
/tmp/tmp.73LABTUPJy: Permission denied.
/bin/mv: cannot stat '/discover/nobackup/projects/gmao/geos_carb/sbasu1/runs/GCM/test_restarts/AGCM.rc.tmpl': No such file or directory
/tmp/tmp.rdwOgZdlbP: Permission denied.
/bin/mv: cannot stat '/discover/nobackup/projects/gmao/geos_carb/sbasu1/runs/GCM/test_restarts/AGCM.rc.tmpl': No such file or directory
/tmp/tmp.Qak33WKR9E: Permission denied.
/bin/mv: cannot stat '/discover/nobackup/projects/gmao/geos_carb/sbasu1/runs/GCM/test_restarts/AGCM.rc.tmpl': No such file or directory
/tmp/tmp.vDLxsZWrwJ: Permission denied.
/bin/mv: cannot stat '/discover/nobackup/projects/gmao/geos_carb/sbasu1/runs/GCM/test_restarts/AGCM.rc.tmpl': No such file or directory
cat: /discover/nobackup/projects/gmao/geos_carb/sbasu1/runs/GCM/test_restarts/AGCM.rc.tmpl.tmp: No such file or directory
```
For some unknown reason, `/tmp` on discover acts up with denied permissions. Probably because it's mounted with `noexec`. To solve, do
```
export TMPDIR=/discover/nobackup/$USER/tmp
mkdir -p $TMPDIR
```
before executing `gcm_setup`.
## Create restart files

This is a dark art. Remembering Robert the Bruce before embarking on this endeavor would be well advised.

GEOS restart files are called `*_rst`, even though they're really netcdf files. Ours not to reason why, ours but to do and die. You will see some `*_import_rst` and some `*_internal_rst`. Ignore the first kind, you will only need to supply the second kind for a new run. There are two types of `*_internal_rst` restart files, upper air (3D) restarts and surface (2D) restarts. Upper air restarts are defined on the cube, contains variables with shape `levels x N x 6N`, and are fairly easily created by the provided scripts for creating/remapping restarts (more below). There are very few ways in which these can "go wrong". Surface restarts can also be created by the provided remapping scripts. However, these will very likely make you weep. Instead of being on grids, surface restarts are provided as a list of tiles (my theory is that whoever made that decision was trying to save disk space and reinvented the wheel instead of relying on compression algorithms). **Every single land model** has a different ordering of these tiles, and understanding what your land model is requires a fair amount of expert knowledge. Worse, the choice of a land model makes pretty much zero difference in a replay run, yet your model will crash unless you do this correctly. In moments of frustration, remember Robert the Bruce.

### Creating restarts using GEOS-provided scripts

The script to create restarts is called `install/bin/remap_restarts.py`. **Do not run this on a compute node** if you are in the `debug` queue because it will try to submit a regrid run to the debug queue, which will fail. Better to just run this on a login node, as follows:
```
module purge
source @env/g5_modules.sh
install/bin/remap_restarts.py
```
This will present you with a series of questions, answer as follows.
* Remap from archived GEOS-IT restarts? `Yes`
* Enter the restart date and hour:  Enter YYYYMMDDHH, where, for GEOS-IT, DD is either 14 or 28, and HH is 21
* Is the upper air input hydrostatic? `Yes`
* Enter output directory for new restarts: Make sure this is a unique folder which is *not* your run folder, you can later copy them over
* Remap to a stretched cubed-sphere grid? `No`
* Enter atmospheric grid for new restarts: Enter the same atmospheric resolution you entered for `gcm_setup`
* Select ocean model for new restarts: `data`
* Select data ocean grid/resolution for new restarts: `CS`
* Enter number of atmospheric levels for new restarts: Choose what you chose for `gcm_setup`
* Select boundary conditions (BCs) version for new restarts: This depends on what you chose for the land boundary condition in `gcm_setup`. If you chose `Icarus-NLv3` there, choose `NL3` here.
* Land BCs for input restarts: You will be presented a folder choice, accept it
* Select BCs base directory for new restarts: Select what you are given
* Land BCs for output restarts: Select what you are given
* Remap upper air restarts? `Yes`
* Remap agcm_import_rst (a.k.a. IAU) file needed for REPLAY runs? `No`
* Remap surface restarts? `Yes`
* Remap bkg files? `No`
* Write lcv file? `No`
* Enter value of WEMIN. No idea what this is, just choose what you are given.
* Enter value of zoom parameter for surface restarts [1-8]? No idea what this is, just choose what you are given.
* Enter experiment ID for new restarts: Fine to leave this blank.
* Add labels for BCs version and atm/ocean resolutions to restart file names? `No`
* SLURM or PBS quality-of-service (qos)? `debug`
* ('Select/enter SLURM or PBS account:\n',) `s1460`
* ('Enter SLURM or PBS partition: (If desired; can leave blank.)\n',) Leave blank.

After entering all the questions, it will write a YAML file with your choices and, if needed, submit a job to the queue to regrid the upper air restarts. It will _make you wait while it does_, i.e., the `sbatch` command won't exit. Don't close the terminal or quit at this point, hopefully the debug queue will be quick enough. Once the job is done, you need to copy over the `*_rst.nc4` files from the ouput folder (above) to your run directory and remove the extension `.nc4`.

## Set up a minimal model to run just `PCHEM`

Before running with GOCART, RRG etc., it's good practice to run a "minimal" replay configuration. Once that works, you can add components and tracers. For this example, we will run with just PCHEM (I believe this stands for Parameterized CHEMistry) for chemistry, and specify that for the source of radiative forcing, aerosols, etc. This obviates the need for tracer restart files and emissions and gets you running the GCM. To do that, make the following modifications.

1. Comment out the following lines in `AGCM.rc`
    ```
    # Enable wet scavenging
    #MCHEMTRI_increments::
    #DU::DU default
    #SS::SS default
    #SU::SO4 default
    #CA.bc::CA.bcphilic default
    #CA.br::CA.brphilic default
    #CA.oc::CA.ocphilic default
    #NI::NO3an1 "NI::NO3an2,NI::NO3an3"
    #PCHEM::OX default
    #::
    ```
The lines might not be exactly as above, so basically the table of `MCHEMTRI_increments`. Since you're running `PCHEM`, you can leave the table uncommented with just the `PCHEM::` components. Also, if your `AGCM.rc` has `MTRI_increments`, comment that table out as well.

2. In `AGCM.rc`, set `NX` to 6 and `NY` to 36 for a `c90` run. This runs the GCM on 216 cores. `NY` must always be `6*NX`, and `NX` has to be chosen so that `AGCM_IM` is at least a factor of 3 larger. Change `NUM_BACKEND_PES` to 32. Each Milan node has 126 cores available, so if you're going to keep a node separate for I/O, no point in using just 16 of its cores.

3. In `RC/GEOS_ChemGridComp.rc`, set everything to `FALSE`, except `ENABLE_PCHEM`. Set that to `TRUE`.

4. In `AGCM.rc`, set the appropriate `RATS` and `AERO` providers,
    ```
    RATS_PROVIDER: PCHEM   # options: PCHEM, GMICHEM, STRATCHEM (Radiatively active tracers)
    AERO_PROVIDER: none    # options: GOCART2G, MAM, none  (Radiatively active aerosols)
    ANALYSIS_OX_PROVIDER: PCHEM   # options: PCHEM, GMICHEM, STRATCHEM, GOCART
    ```

5. In `AGCM.rc` set `USE_AEROSOL_NN: .false.`. That key may not exist in recent model tags, add it.

6. In `RC/GOCART2G_GridComp.rc`, keep all `ACTIVE_INSTANCES_*` and `PASSIVE_INSTANCES_*` blank, e.g.,
    ```
    ACTIVE_INSTANCES_DU:
    PASSIVE_INSTANCES_DU:

    ACTIVE_INSTANCES_SS:
    PASSIVE_INSTANCES_SS:

    ACTIVE_INSTANCES_SU:
    PASSIVE_INSTANCES_SU:

    ACTIVE_INSTANCES_CA:
    PASSIVE_INSTANCES_CA:

    ACTIVE_INSTANCES_NI:
    PASSIVE_INSTANCES_NI:
    ```

7. In `HISTORY.rc` do not ask for any collection to be written, i.e.,
    ```
    COLLECTIONS:
    ::
    ```

## Tell the model when to run and for how long

In the run directory, create a file called `cap_restart`. This tells the model when to start. The format is `YYYYMMDD HHMMSS` with exactly one newline at the end. This should match the time you specified when creating the restart files above. E.g., if you created restart files for 2015-09-28 21z, the `cap_restart` should contain a single string, `20150928 210000`, with a newline at the end.

Also in the run directory, there is a file called `CAP.rc`. Set the `END_DATE:` in there to be (say) 18 days after `cap_restart`. Select `JOB_SGMT` to be `00000005 000000` to run for 5 days at a time. Select `NUM_SGMT` to be 2 to run two of those segments in each job. This _should_ result in two consecutive jobs, one running for 10 days (5+5), the other running for 8 days (5+3), before hitting the `END_DATE` you specified.
