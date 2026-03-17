✦ I've successfully integrated checkpoint captures into the Metabolomics and Lipidomics workflows, mirroring the Proteomics architecture. I've added a "Capture unit test checkpoints"
  checkbox to the import UIs and inserted .capture_checkpoint calls at key processing steps (cp01, cp02, cp05, cp06, cp07, cp08).


  To make "vibe coding" easier for linking these fixtures to tests, I recommend using a convention-based directory structure like tests/testdata/<dataset>/<omics>/cp<NN>_<label>.rds. I
  also suggest creating a tests/testthat/helper-checkpoints.R file with a get_checkpoint function to dynamically resolve paths, and including a manifest.json in each checkpoint directory
  to provide context-aware metadata for an AI to use during test generation.