# this is example-file: Makefile.am

# the subdirectories of the project to go into
SUBDIRS =   \
            src

# Test fixtures keep one lowercase input file; test copies use the solver's
# required uppercase filename to avoid case-only aliases in Git.
EXTRA_DIST = \
  COPYING \
  CITATION.cff \
  CONTRIBUTING.md \
  README.md \
  TESTS/run_tests.sh \
  TESTS/TEST_SOURCE_REGISTRY/README.md \
  TESTS/TEST_SOURCE_REGISTRY/run_test.sh \
  TESTS/TEST_SOURCE_REGISTRY/test_source_registry.f90 \
  TESTS/TEST_RHEOLOGY_CELL_LOCAL/README.md \
  TESTS/TEST_RHEOLOGY_CELL_LOCAL/run_test.sh \
  TESTS/TEST_RHEOLOGY_CELL_LOCAL/test_source_modules.f90 \
  TESTS/TEST_RHEOLOGY_CELL_LOCAL/test_analytic_velocity.f90 \
  TESTS/TEST_RHEOLOGY_CELL_LOCAL/test_lava_cooling.f90 \
  TESTS/TEST_RHEOLOGY_CELL_LOCAL/test_pore_pressure_diffusion.f90 \
  TESTS/TEST_RHEOLOGY_CELL_LOCAL/test_model_names.f90 \
  TESTS/TEST_RHEOLOGY_CELL_LOCAL/check_inputs.py \
  TESTS/TEST_RHEOLOGY_CELL_LOCAL/check_stages.py \
  TESTS/TEST_SOURCE_MASS/README.md \
  TESTS/TEST_SOURCE_MASS/run_test.sh \
  TESTS/TEST_SOURCE_MASS/make_cases.py \
  TESTS/TEST_SOURCE_MASS/check_mass.py \
  TESTS/TEST_OUTPUT_PATHS/run_test.sh \
  TESTS/TEST_OUTPUT_PATHS/test_snapshot_index.f90 \
  TESTS/TEST_2D/IMEX_sfloW2D.inp \
  TESTS/TEST_2D/INIT_2D_0000.q_2d \
  TESTS/TEST_2D/cleanFolder.sh \
  TESTS/TEST_2D/topography_dem.asc \
  TESTS/TEST_BOX/IMEX_sfloW2D.inp \
  TESTS/TEST_BOX/INIT_BOX_0000.q_2d \
  TESTS/TEST_BOX/cleanFolder.sh \
  TESTS/TEST_BOX/topography_dem.asc \
  TESTS/TEST_BW/IMEX_sfloW2D.inp \
  TESTS/TEST_BW/INIT_BW_0000.q_2d \
  TESTS/TEST_BW/cleanFolder.sh \
  TESTS/TEST_BW/topography_dem.asc \
  TESTS/TEST_CF/IMEX_sfloW2D.inp \
  TESTS/TEST_CF/cleanFolder.sh \
  TESTS/TEST_TAAL/IMEX_sfloW2D.inp \
  TESTS/TEST_TAAL/cleanFolder.sh \
  TESTS/TEST_DAMBREAK/TEST_DYNAMIC_FINH/IMEX_SfloW2D.inp \
  TESTS/TEST_DAMBREAK/TEST_DYNAMIC_FINH/cleanFolder.sh \
  TESTS/TEST_DAMBREAK/TEST_DYNAMIC_FINH/example_2D_0000.q_2d \
  TESTS/TEST_DAMBREAK/TEST_DYNAMIC_FINH/topography_dem.asc \
  TESTS/TEST_DAMBREAK/TEST_OFF_FINH/IMEX_SfloW2D.inp \
  TESTS/TEST_DAMBREAK/TEST_OFF_FINH/cleanFolder.sh \
  TESTS/TEST_DAMBREAK/TEST_OFF_FINH/example_2D_0000.q_2d \
  TESTS/TEST_DAMBREAK/TEST_OFF_FINH/topography_dem.asc \
  TESTS/TEST_DAMBREAK/TEST_STATIC_FINH/IMEX_SfloW2D.inp \
  TESTS/TEST_DAMBREAK/TEST_STATIC_FINH/cleanFolder.sh \
  TESTS/TEST_DAMBREAK/TEST_STATIC_FINH/example_2D_0000.q_2d \
  TESTS/TEST_DAMBREAK/TEST_STATIC_FINH/topography_dem.asc \
  TESTS/TEST_DAMBREAK/TEST_STEP_FINH/IMEX_SfloW2D.inp \
  TESTS/TEST_DAMBREAK/TEST_STEP_FINH/cleanFolder.sh \
  TESTS/TEST_DAMBREAK/TEST_STEP_FINH/example_2D_0000.q_2d \
  TESTS/TEST_DAMBREAK/TEST_STEP_FINH/topography_dem.asc \
  EXAMPLES/EXAMPLE_CF/topo.zip \
  EXAMPLES/EXAMPLE_TAAL/topography_dem.asc \
  TESTS/DIAGNOSTICS/README.md \
  TESTS/DIAGNOSTICS/HP_PCCU_1D/README.md \
  TESTS/DIAGNOSTICS/HP_PCCU_1D/hp_pccu_1d_core.py \
  TESTS/DIAGNOSTICS/HP_PCCU_1D/kt_v2_1d.py \
  TESTS/DIAGNOSTICS/HP_PCCU_1D/run_inclined_1d.py \
  TESTS/DIAGNOSTICS/HP_PCCU_1D/smooth_bump_reconstruction_1d.py \
  TESTS/DIAGNOSTICS/HP_PCCU_1D/fortran_inclined_1d/README.md \
  TESTS/DIAGNOSTICS/HP_PCCU_1D/fortran_inclined_1d/create_case.py \
  TESTS/DIAGNOSTICS/HP_PCCU_1D/fortran_inclined_1d/IMEX_SfloW2D.inp \
  TESTS/DIAGNOSTICS/HP_PCCU_1D/fortran_inclined_1d/inclined_excavation_1D_0000.q_2d \
  TESTS/DIAGNOSTICS/HP_PCCU_1D/fortran_inclined_1d/topography_dem.asc \
  TESTS/DIAGNOSTICS/ARCHIVE/README.md \
  TESTS/DIAGNOSTICS/ARCHIVE/kt_v2_original_1d.py \
  TESTS/DIAGNOSTICS/ARCHIVE/kt_ritter_1d.f90

check-local:
	$(SHELL) $(top_srcdir)/TESTS/TEST_OUTPUT_PATHS/run_test.sh \
	  "$(FC)" "$(OPENMP_FCFLAGS)" "$(abs_top_builddir)/src"
	$(SHELL) $(top_srcdir)/TESTS/run_tests.sh \
	  "$(abs_top_builddir)/src/IMEX_SfloW2D"
