# CHyM_cpl

CHyM code for RegESM coupling. The code in this directory contains the code
to be used in the RegCM-ES1-1 coupled model. The model description is
identical to the model available in stand-alone version in [CHyM-roff](https://github.com/graziano-giuliani/CHyM-roff).

The difference is the removal of the MPI parallelization code, the absence
of any pre processing program, the lack of input data part (coming in from the
[RegESM](https://github.com/graziano-giuliani/RegESM) coupling) and the
very simplified output module.

The eventually binary created program is just a placeholder to verify in debug
all is correctly working. The model can actually work ONLY in the coupled
settings.

Refer to the "parent" code CHyM-roff on how to prepare the static dataset
needed for the model use. Assuming the model setting are the same, the model
result is expected to be the same as the parent CHyM-roff.

The model needs a similar namelist as the one used for the CHyM-roff at runtime.

The Makefile in this directory is the same used on the ICTP MED-cordex target
platform, the CINECA leonardo pre-exascale system.
