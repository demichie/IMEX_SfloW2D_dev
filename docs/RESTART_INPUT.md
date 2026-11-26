# Initial-state and restart input groups

The `RESTART` run flag selects initialization from files. The file lists
described below are for initial-state input, **not** successive times of a
binary checkpoint series. The binary `restart.bin` is one complete snapshot
of the whole simulation, including all layers and runtime state.

## Single layer

Exactly one of `RESTART_PARAMETERS` or `RESTART_L1_PARAMETERS` is required
when `RESTART=T`. Both spellings bind the same L1 file list and release
schedule; do not provide both. A single-layer run uses the corresponding
`INITIAL_PARAMETERS` or `INITIAL_L1_PARAMETERS` for primitive values.

```fortran
&RESTART_PARAMETERS
 N_RESTART_FILES=2,
 RESTART_FILES='first.asc','second.asc',
 RELEASE_TIME=0.D0,10.D0,
/
```

The existing single-layer `RELEASE_TIME` schedule remains unchanged.

## Two layers: text initial states

Provide `RESTART_L1_PARAMETERS` and `RESTART_L2_PARAMETERS` separately.
Each group has its own list and count; the L2 members carry the `_L2`
suffix, matching the `INITIAL_L2_PARAMETERS` input convention.

```fortran
&RESTART_L1_PARAMETERS
 N_RESTART_FILES=2,
 RESTART_FILES='first_l1.asc','second_l1.asc',
 RELEASE_TIME=0.D0,0.D0,
/
&RESTART_L2_PARAMETERS
 N_RESTART_FILES_L2=1,
 RESTART_FILES_L2='NONE',
 RELEASE_TIME_L2=0.D0,
/
```

The example accumulates two same-time thickness rasters into L1 and leaves L2
dry. A layer with no initial material has a single `'NONE'` entry. Two
`NONE` entries, one per group, initialize both layers dry. All nonempty
files in one primitive initialization must use `.asc` format. Composition,
temperature and initial velocity are defined in `INITIAL_L1_PARAMETERS`
and `INITIAL_L2_PARAMETERS`.

**Two-layer delayed releases are not yet supported.** The solver rejects
nonzero `RELEASE_TIME` or `RELEASE_TIME_L2` for two-layer primitive
initialization. The current coupled state does not defer material activation
independently per layer. Accepting such input without a separate activation
mechanism could generate nonphysical interlayer coupling before release.
This explicit restriction will be removed only after scientifically
validated layer-specific activation is implemented.

For conservative `.q_2d` input, provide exactly one complete snapshot file
for each layer, at the same output index. Do not combine `.q_2d` files with
`.asc` rasters. These snapshots are initial states, not multiple releases.

## Binary checkpoint

A `.bin` checkpoint contains both layer states in one file. It is supplied
via `RESTART_PARAMETERS` or `RESTART_L1_PARAMETERS`; in this case
`RESTART_L2_PARAMETERS` must be absent because the file already contains
the complete canonical multi-layer state and restart metadata.

The input reader rejects duplicate L1 group spellings, missing entries,
unused extra file entries and mixed input formats. This is a development
input contract; there are no implicit legacy format fallbacks.
