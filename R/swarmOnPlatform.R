#' runSwarmMultiprocess
#'
#'
#' @name runSwarmMultiprocess
#' @rdname runSwarmMultiprocess
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param swarmInit See \code{\link[CIMseq]{CIMseqSwarm}}.
#' @param maxiter See \code{\link[CIMseq]{CIMseqSwarm}}.
#' @param nSyntheticMultiplets See \code{\link[CIMseq]{CIMseqSwarm}}.
#' @param eps.stagnate See \code{\link[CIMseq]{pso.2.0}}.
#' @param maxit.stagnate See \code{\link[CIMseq]{pso.2.0}}.
#' @param currPath The analysis directory.
#' @param ... additional arguments to pass on
#' @return Saves swarm output to data folder in current directory.
#' @author Jason T. Serviss
#'
#'
#'
NULL
#' @export
#' @import CIMseq
#' @import future.apply

runSwarmMultiprocess <- function(
  singlets, multiplets, 
  swarmInit = NULL, 
  maxiter = 100,
  nSyntheticMultiplets = 400,
  eps.stagnate = 1,
  maxit.stagnate = 5,
  currPath = getwd(),
  ...
){
  if(is.null(swarmInit)) swarmInit <- swarmInit(singlets, 2)
  future::plan(multiprocess)
  print(paste0("Starting deconvolution at ", Sys.time()))
  sObj <- CIMseqSwarm(
    singlets, multiplets, maxiter = maxiter, swarmsize = ncol(swarmInit), 
    nSyntheticMultiplets = nSyntheticMultiplets, swarmInit = swarmInit,
    psoControl = list(eps.stagnate = eps.stagnate, maxit.stagnate = maxit.stagnate)
  )
  print(paste0("Finished deconvolution at ", Sys.time()))
  save(sObj, file = file.path(currPath, "data/sObj.rda"))
}

.runSwarmBatch <- function(
  singlets, multiplets, 
  swarmInit = NULL, 
  maxiter = 100,
  nSyntheticMultiplets = 400,
  eps.stagnate = 1,
  maxit.stagnate = 5,
  currPath = getwd(),
  time = "24:00:00",
  ...
){
  if(is.null(swarmInit)) swarmInit <- swarmInit(singlets, 2)
  options(future.wait.interval = 10000.0)
  options(future.wait.timeout = 1e9)
  future::plan(
    future.batchtools::batchtools_slurm,
    template = "/crex/proj/snic2018-8-151/private/batchtools.slurm.tmpl",
    resources = list(
      account = "snic2018-8-151", partition = "core", ntasks = 1L,
      time = time, jobname = "CIMseq",
      modules = "R_packages/3.5.0", R = "R/3.5.0", 
      log.file = file.path(currPath, "logs/slurm.txt")
    ),
    workers = 100
  )
  
  #run deconvolution
  print(paste0("Starting deconvolution at ", Sys.time()))
  sObj <- CIMseqSwarm(
    cObjSng, cObjMul, maxiter = 100, swarmsize = 500, nSyntheticMultiplets = 400
  )
  print(paste0("Finished deconvolution at ", Sys.time()))
  save(sObj, file = file.path(currPath, "data/sObj.rda"))
}

#' runSwarmUppmax
#'
#'
#' @name runSwarmUppmax
#' @rdname runSwarmUppmax
#' @param singlets CIMseqSinglets; A CIMseqSinglets object.
#' @param multiplets CIMseqMultiplets; A CIMseqMultiplets object.
#' @param swarmInit See \code{\link[CIMseq]{CIMseqSwarm}}.
#' @param maxiter See \code{\link[CIMseq]{CIMseqSwarm}}.
#' @param nSyntheticMultiplets See \code{\link[CIMseq]{CIMseqSwarm}}.
#' @param eps.stagnate See \code{\link[CIMseq]{pso.2.0}}.
#' @param maxit.stagnate See \code{\link[CIMseq]{pso.2.0}}.
#' @param currPath The analysis directory.
#' @param args List of arguments passed to the script. Should have sample in
#'  \code{args[1]}.
#' @param ... additional arguments to pass on
#' @return Saves swarm output to data folder in current directory.
#' @author Jason T. Serviss
#'
#'
#'
NULL
#' @export
#' @import CIMseq
#' @import future.apply

runSwarmUppmax <- function(
  singlets, multiplets, 
  swarmInit = NULL, 
  maxiter = 100,
  nSyntheticMultiplets = 400,
  eps.stagnate = 1,
  maxit.stagnate = 5,
  currPath = getwd(),
  args,
  ...
){
  plan(sequential)
  sample <- as.character(args[1])
  print(paste0("Running ", sample))
  counts <- getData(multiplets, "counts")
  counts.ercc <- getData(multiplets, "counts.ercc")
  multiplets.2 <- CIMseqMultiplets(
    matrix(counts[, sample], ncol = 1, dimnames = list(rownames(counts), sample)),
    matrix(counts.ercc[, sample], ncol = 1, dimnames = list(rownames(counts.ercc), sample)),
    getData(multiplets, "features")
  )
  
  #run deconvolution
  if(is.null(swarmInit)) swarmInit <- swarmInit(singlets, 2)
  print(paste0("Starting deconvolution at ", Sys.time()))
  sObj <- CIMseqSwarm(
    singlets, multiplets.2, maxiter = maxiter, swarmsize = ncol(swarmInit), 
    nSyntheticMultiplets = nSyntheticMultiplets, swarmInit = swarmInit,
    psoControl = list(eps.stagnate = eps.stagnate, maxit.stagnate = maxit.stagnate)
  )
  print(paste0("Finished deconvolution at ", Sys.time()))
  save(sObj, file = file.path(currPath, paste0("tmp/sObj_", sample, "_uppmax.rda")))
}