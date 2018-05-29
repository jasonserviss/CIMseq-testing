
# Function to create simulated derivative new cell types based on a rofile, by switching genes and adding some random noise.
createCelltypeFromProfile <- function(origProfile, nSwitched, mulNoise = 0, linNoise = 0) {
    newProfile <- origProfile
    if(nSwitched > 0) {
        switchIdx <- sample(1:(length(origProfile)), nSwitched*2, replace=F)
        from <- switchIdx[1:(length(switchIdx)/2)]
        to <- switchIdx[-1:-(length(switchIdx)/2)]
        newProfile[from] <- origProfile[to]
        newProfile[to] <- origProfile[from]
    }
    newProfile <- newProfile * (runif(n=length(newProfile), min=1-mulNoise, max=1+mulNoise)) + rbinom(n=length(newProfile), size=linNoise, prob=0.5)-linNoise/2# Add some noise - multiplicative and additive.
    newProfile[newProfile < 0] <- 0 # Make sure we have no negative values. 
    newProfile
}

# Simulate a cell from a cell type (probability profile) by randomly picking reads with a prob. based on the profile
cellFromProfile <- function(profile, nreads) {
    tbl <- table(sample(1:length(profile), prob=profile, size=nreads, replace=T))
    ret <- profile * 0
    ret[as.integer(names(tbl))] <- tbl
    ret
}

# Merge two profiles with specified fractions.
makeDoubletProfile <- function(prof1, prof2, frac1, frac2) {
    frac1 <- frac1/sum(frac1)
    frac2 <- frac2/sum(frac2)
    rowSums(cbind(prof1*frac1, prof2*frac2))
}


# Make cell type profiles based on single cell data of cell lines
load("counts_s.rda")
rows <- as.integer(substr(sapply(strsplit(colnames(counts.s), "\\."), function(x) {x[2]}), 2, 100))

counts.s.norm <- apply(counts.s, 2, function(x) {x/sum(x) * 1000000}) # Normalize for downstream

ngenes <- 1000
selGenes <- order(apply(counts.s.norm, 1, function(x) {max(x)}), decreasing=T)[1:1000] # Select top 1000 genes

# Create cell type profiles
Celltype1 <- rowMeans(counts.s.norm[selGenes, rows %in% 1:4])
Celltype2 <- rowMeans(counts.s.norm[selGenes, rows %in% 5:8])
Celltype3 <- rowMeans(counts.s.norm[selGenes, rows %in% 9:12])

targetDistr <- sort(Celltype1 - Celltype2) # This is the distribution of diffs that we are aiming for


syntCellTypes <- lapply(1:10, function(x) {createCelltypeFromProfile(Celltype1, 100, mulNoise=0.2, linNoise=200)})

plot(targetDistr, col="red", type='n', ylim=c(-4000, 4000))
for(i in 1:(length(syntCellTypes)-1) ) {
    for(j in (i+1):length(syntCellTypes)) {
        lines(sort(syntCellTypes[[i]] - syntCellTypes[[j]]), col="grey")
    }
}
lines(sort(Celltype2 - Celltype3), col="green")
lines(sort(Celltype1 - Celltype3), col="red")
lines(targetDistr, col="black")
legend("topleft", fill=c("black", "grey", "red", "green"), legend=c("real", "simulated", "real", "real"))
dev2bitmap("Distr_of_differences.pdf", type="pdfwrite")

# Create a dataset with simulated cell types
newPops <- lapply(1:10, function(x) {
    sapply(1:30, function(y) {
        cellIdx <- sample(1:(ncol(counts.s)),1)
        rawProf <- cellFromProfile(syntCellTypes[[x]], nreads=sum(counts.s[selGenes,cellIdx])) # generate profile with sensible number of reads
        rawProf * sum(counts.s.norm[selGenes,cellIdx]) / sum(counts.s[selGenes,cellIdx]) # Normalize
    })
})
newPops <- do.call(cbind, newPops) # newPops now has the simulated count data.


a <- tsne(dist(t(cbind(counts.s.norm[selGenes,], newPops))), max_iter=600, epoch_callback=plot)
plot(a, col=c(rep("red", ncol(counts.s.norm)), rep("black", nrow(a)-ncol(counts.s.norm))), pch=19)
legend("topleft", fill=c("red", "black"), legend=c("real", "simulated"))
dev2bitmap("tSNE_of_original_and_simulated.pdf", type="pdfwrite")



###################################################################
# Create doublets based on the above cell types - one real connection and the rest random.
minFrac <- 0.3
maxFrac <- 0.7
nSpecDbl <- 50
realCon.doublets <- sapply(1:nSpecDbl, function(y) { # Make specific doublets
    frac1 <- runif(n=1, min=minFrac, max=maxFrac)
    frac2 <- 1-frac1
    realCon <- makeDoubletProfile(prof1 = syntCellTypes[[1]], prof2=syntCellTypes[[2]], frac1=frac1, frac2=frac2)
    cellIdx <- sample(1:(ncol(counts.s)),1)
    rawProf <- cellFromProfile(realCon, nreads=sum(counts.s[selGenes,cellIdx])) # generate profile with sensible number of reads
    rawProf * sum(counts.s.norm[selGenes,cellIdx]) / sum(counts.s[selGenes,cellIdx]) # Normalize
})


# Make random doublets
nRandDbl <- 200
randomCon.doublets <- sapply(1:nRandDbl, function(x) {
    cellIdx <- sample(1:(ncol(counts.s)),1)
    Celltype1 <- syntCellTypes[[sample(1:length(syntCellTypes), 1)]]
    Celltype2 <- syntCellTypes[[sample(1:length(syntCellTypes), 1)]]
    frac1 <- runif(n=1, min=minFrac, max=maxFrac)
    frac2 <- 1-frac1
    dblProf <- makeDoubletProfile(prof1=Celltype1, prof2=Celltype2, frac1=frac1, frac2=frac2)
    rawProf <- cellFromProfile(dblProf, nreads=sum(counts.s[selGenes,cellIdx])) # generate profile with sensible number of reads
    rawProf * sum(counts.s.norm[selGenes,cellIdx]) / sum(counts.s[selGenes,cellIdx]) # Normalize
})

all.doublets <- cbind(realCon.doublets, randomCon.doublets)

aa <- tsne(dist(t(cbind(all.doublets))), max_iter=1000, epoch_callback=function(x) {plot(x, col=c(rep("black", nSpecDbl), rep("red", nRandDbl)))})
plot(aa, col=c(rep("red", nSpecDbl), rep("black", nRandDbl)), pch=19)
legend("topleft", fill=c("red", "black"), legend=c("Specific doublets", "Random doublets"))
dev2bitmap("Doublets_tSNE.pdf", type="pdfwrite")
 
