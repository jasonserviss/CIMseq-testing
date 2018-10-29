
packages <- c("CIMseq", "sp.scRNAseqData", "tidyverse", "future", "future.apply")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

currPath <- getwd()

if(file.exists(file.path(currPath, '../mouseAnalysis_engeOnly/data/CIMseqData.rda'))) {
  load(file.path(currPath, '../mouseAnalysis_engeOnly/data/CIMseqData.rda'))
}

print("Starting permutation")
nPerms <- 5
perms <- map(1:nPerms, function(x) {
  print(x)
  #plan(multiprocess)

  future::plan(
    future.batchtools::batchtools_slurm,
    template = "/crex/proj/snic2018-8-151/private/batchtools.slurm.tmpl",
    resources = list(
      account = "snic2018-8-151", partition = "core", ntasks = 1L,
      time = "24:00:00", jobname = "SCM.permute",
      modules = "R_packages/3.5.0", R = "R/3.5.0", log.file = file.path(currPath, "logs/slurm.txt")
    ),
    workers = 100
  )

  spSwarm(
    cObjSng, cObjMul, uObj, maxiter = 100, swarmsize = 500,
    nSyntheticMultiplets = 400, permute = TRUE, seed = x
  )
})

print("Finished permutation")
save(perms, file = file.path(currPath, "data/permutations.rda"))
writeLines(capture.output(sessionInfo()), file.path(currPath, "logs/sessionInfo_spSwarm.txt"))
print("Finished")
