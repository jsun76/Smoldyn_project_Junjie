# Actin load adaptation project

## Starting point

This work starts from `MatsulabUW/Smoldyn` branch
`akamatsu/filament-capabilities`. The first target is to reproduce the lab's
branched actin compression baseline before changing model parameters or code.
The experimental comparison is Bieling et al., *Cell* (2016),
doi:10.1016/j.cell.2015.11.057.

## Baseline inputs

The lab's `Shared data / Smoldyn` folder contains the model and reference
artifacts: `compressionPoC3D_baseline.txt`, `baseline_run.log`,
`baseline.png`, `analyze_compression.py`, and
`baseline_stride10.simularium`. Keep generated outputs and large trajectories
out of Git unless the lab chooses a storage policy for them.

## Reproduction order

1. Build Smoldyn from the lab branch and record the exact commit and build
   settings.
2. Run `compressionPoC3D_baseline.txt` with the recorded random seed and
   save the command, standard output, and generated files.
3. Calculate SHA-256 hashes for the model, executable, and outputs; compare
   output hashes with the lab's expected values when provided.
4. Run `analyze_compression.py` and compare its figure with `baseline.png`.
5. Convert the trajectory to Simularium format and inspect it alongside
   `baseline_stride10.simularium` in the lab viewer.

The first paper-facing comparisons should use like-for-like definitions and
units for network growth velocity, actin density, and applied stress. Record
each parameter's source before calibrating it.

## Commit and notes routine

- Make one focused commit for each coherent change to model inputs, analysis,
  code, or documentation. Include the reason and validation in its message.
- Before each commit, review `git status` and `git diff`, then run the relevant
  example or check. Do not commit generated trajectories or local build files.
- Put the commit SHA, run command, input/output hashes, and result summary in
  the corresponding Roam note. Push the personal branch after each validated
  milestone once a writable remote is available.
- Fetch the lab branch regularly and review its changes before merging them
  into the personal branch.
