

#############################################################################
#############################################################################
##### These are the old functions for simulating a clade and turning it into
##### a phylogeny from paleotree. The new ones are far too slow for the size
##### of data we are dealing with here.
#############################################################################
#############################################################################

simFossilTaxa<-function (p, q, anag.rate = 0, prop.bifurc = 0, prop.cryptic = 0, 
    nruns = 1, mintaxa = 1, maxtaxa = 1000, mintime = 1, maxtime = 1000, 
    minExtant = 0, maxExtant = NULL, min.cond = TRUE, count.cryptic = FALSE, 
    print.runs = FALSE, sortNames = FALSE, plot = FALSE) 
{
    if (any(c(p, q, anag.rate, prop.bifurc, prop.cryptic) < 0)) {
        stop("Error: bad parameters input, p, q, anag.rate, prop.bifurc or prop.cryptic are less than 0")
    }
    if (prop.bifurc > 0 & prop.cryptic == 1) {
        stop("Error: Prop.bifurc greater than 0 even though cryptic cladogenesis = 1??")
    }
    if (nruns < 1) {
        stop("Error: nruns<1")
    }
    if (maxtaxa < 0) {
        stop("Error: maxtaxa<0")
    }
    if (mintaxa < 1) {
        stop("Error: mintaxa<1")
    }
    if (mintime < 1) {
        stop("Error: mintime<1")
    }
    if (maxtime < mintime) {
        stop("Error: maxtime<mintime")
    }
    if (mintaxa > maxtaxa) {
        stop("Error: mintaxa > maxtaxa")
    }
    if (maxtaxa > 10000 & maxtime > 10000) {
        warning("Warning: Unrealistic limits for maxtaxa or maxtime")
    }
    if (minExtant < 0) {
        stop("Error: minExtant<0")
    }
    if (minExtant > mintaxa) {
        mintaxa <- minExtant
    }
    if (!is.null(maxExtant)) {
        if (maxExtant < 0) {
            stop("Error: maxExtant<0")
        }
        if (maxExtant > maxtaxa) {
            maxtaxa <- maxExtant
        }
        if (minExtant > maxExtant) {
            stop("Error: maxExtant is set higher than minExtant")
        }
    }
    if (!min.cond) {
        message("No conditioning during simulation; run until max limits or total extinction")
    }
    minExtant1 <- ifelse(minExtant == 0, 0, minExtant + 1)
    results <- list()
    ntries <- 0
    for (i in 1:nruns) {
        ntries <- ntries + 1
        taxad <- matrix(c(1, NA, 0, NA, 1), 1, )
        pqw <- p + q + anag.rate
        maxtime1 <- maxtime
        continue <- TRUE
        eval <- FALSE
        while (any(is.na(taxad[, 4])) & continue) {
            tpot <- is.na(taxad[, 4])
            tpot2 <- min(taxad[tpot, 3]) == taxad[, 3]
            tpick <- which(tpot & tpot2)[1]
            tpick_FO <- taxad[tpick, 3]
            wait <- 0
            while (is.na(taxad[tpick, 4]) & continue) {
                wait <- rexp(1, rate = pqw) + wait
                type <- sample(1:3, 1, prob = c(p/pqw, q/pqw, 
                  anag.rate/pqw))
                if (type == 1) {
                  type1 <- sample(1:3, 1, prob = c((1 - prop.cryptic) - 
                    ((1 - prop.cryptic) * prop.bifurc), (1 - 
                    prop.cryptic) * prop.bifurc, prop.cryptic))
                  if (type1 == 1) {
                    taxad <- rbind(taxad, c(max(taxad[, 1]) + 
                      1, taxad[tpick, 1], wait + tpick_FO, NA, 
                      max(taxad[, 1]) + 1))
                  }
                  if (type1 == 2) {
                    taxad[tpick, 4] <- wait + tpick_FO
                    taxad <- rbind(taxad, c(max(taxad[, 1]) + 
                      1, taxad[tpick, 1], wait + tpick_FO, NA, 
                      max(taxad[, 1]) + 1))
                    taxad <- rbind(taxad, c(max(taxad[, 1]) + 
                      1, taxad[tpick, 1], wait + tpick_FO, NA, 
                      max(taxad[, 1]) + 1))
                  }
                  if (type1 == 3) {
                    taxad <- rbind(taxad, c(max(taxad[, 1]) + 
                      1, taxad[tpick, 1], wait + tpick_FO, NA, 
                      taxad[tpick, 5]))
                  }
                }
                if (type == 2) {
                  taxad[tpick, 4] <- wait + tpick_FO
                }
                if (type == 3) {
                  taxad[tpick, 4] <- wait + tpick_FO
                  taxad <- rbind(taxad, c(max(taxad[, 1]) + 1, 
                    tpick, wait + tpick_FO, NA, max(taxad[, 1]) + 
                      1))
                }
                if (count.cryptic) {
                  numtax <- nrow(taxad)
                }
                else {
                  numtax <- length(unique(taxad[, 5]))
                }
                if (numtax > maxtaxa) {
                  maxtime1 <- min(c(maxtime1, taxad[maxtaxa + 
                    1, 3]))
                  if (!min.cond) {
                    eval <- TRUE
                  }
                }
                if (wait > maxtime1) {
                  taxad[tpick, 4] <- wait
                }
                if (sum(taxad[, 3] < maxtime1) > (maxtaxa * 2)) {
                  taxad[is.na(taxad[, 4]), 4] <- maxtime + 1
                }
                if (count.cryptic) {
                  numtax <- nrow(taxad)
                  numext <- sum(is.na(taxad[, 4])) + sum(taxad[!is.na(taxad[, 
                    4]), 4] >= maxtime1)
                }
                else {
                  numtax <- length(unique(taxad[, 5]))
                  numext <- length(unique(taxad[(is.na(taxad[, 
                    4]) | taxad[, 4] >= maxtime1), 5]))
                }
                if (max(taxad[, 3:4], na.rm = TRUE) >= mintime & 
                  ifelse(numext > 0, numtax > mintaxa, numtax >= 
                    mintaxa) & numext >= minExtant1 & ifelse(is.null(maxExtant), 
                  TRUE, numext <= maxExtant) & min.cond) {
                  if (any(is.na(taxad[, 4]))) {
                    maxtime2 <- min(c(maxtime1, max(taxad[is.na(taxad[, 
                      4]) | taxad[, 4] >= maxtime1, 3])))
                    if (maxtime2 > mintime) {
                      maxtime1 <- maxtime2
                    }
                  }
                  eval <- TRUE
                }
                continue <- ifelse(any(is.na(taxad[, 4])), any(taxad[is.na(taxad[, 
                  4]), 3] <= maxtime1), FALSE)
                if (!continue & !min.cond) {
                  eval <- TRUE
                }
                taxad_save <- taxad
            }
            if (!continue & eval) {
                taxad <- matrix(taxad[taxad[, 3] < maxtime1, 
                  ], sum(taxad[, 3] < maxtime1), )
                if (any(is.na(taxad[, 4]))) {
                  stop("Error: Live creatures escaping simulation! Get out now while you still have time!")
                }
                posstimes <- sort(unique(c(taxad[, 3:4], maxtime1)))
                maxtimes <- posstimes[posstimes >= mintime & 
                  posstimes <= maxtime1]
                if (length(maxtimes) == 0) {
                  eval <- FALSE
                }
                else {
                  mtds <- lapply(maxtimes, function(x) matrix(taxad[taxad[, 
                    3] < x, ], sum(taxad[, 3] < x), ))
                  if (count.cryptic) {
                    numtaxa <- sapply(mtds, function(x) nrow(x))
                    numexta <- sapply(1:length(mtds), function(x) sum(mtds[[x]][, 
                      4] >= maxtimes[x]))
                  }
                  else {
                    numtaxa <- sapply(mtds, function(x) length(unique(x[, 
                      5])))
                    numexta <- sapply(1:length(mtds), function(x) length(unique(mtds[[x]][mtds[[x]][, 
                      4] >= maxtimes[x], 5])))
                  }
                  minta <- numtaxa >= mintaxa
                  maxta <- numtaxa <= maxtaxa
                  minti <- maxtimes >= mintime
                  maxti <- maxtimes <= maxtime
                  maxext <- if (!is.null(maxExtant)) {
                    maxExtant >= numexta
                  }
                  else {
                    TRUE
                  }
                  minext <- minExtant <= numexta
                  evalcond <- maxext & minext & minta & maxta & 
                    minti & maxti
                  if (any(evalcond)) {
                    chosen <- rev(which(evalcond))[1]
                    taxad <- mtds[[chosen]]
                    maxtime1 <- maxtimes[chosen]
                  }
                  else {
                    eval <- FALSE
                  }
                }
            }
            if (!continue & !eval) {
                taxad <- matrix(c(1, NA, 0, NA, 1), 1, )
                ntries <- ntries + 1
                continue <- TRUE
                eval <- FALSE
                evalcond <- NULL
                maxtime1 <- maxtime
            }
        }
        taxad1 <- cbind(taxad[, 1:4, drop = FALSE], (taxad[, 
            4] >= maxtime1), taxad[, 5])
        taxad1[, 3:4] <- maxtime1 - taxad1[, 3:4]
        taxad1[, 3] <- round(taxad1[, 3], digits = 4)
        taxad1[, 4] <- round(taxad1[, 4], digits = 4)
        taxad1[taxad1[, 4] < 0, 4] <- 0
        taxad2 <- cbind(matrix(1:nrow(taxad1), , 1), matrix(match(taxad1[, 
            2], taxad1[, 1]), , 1), matrix(taxad1[, 3:5], , 3), 
            matrix(match(taxad1[, 6], taxad1[, 1]), , 1))
        colnames(taxad2) <- c("taxon.id", "ancestor.id", "orig.time", 
            "ext.time", "still.alive", "looks.like")
        names <- paste("t", taxad2[, 1], sep = "")
        if (any(taxad2[, 6] != taxad2[, 1])) {
            for (cry in which(sapply(taxad2[, 6], function(x) sum(x == 
                taxad2[, 6]) > 1))) {
                names[cry] <- paste("t", taxad2[cry, 6], ".", 
                  sum(taxad2[1:cry, 6] == taxad2[cry, 6]), sep = "")
            }
        }
        rownames(taxad2) <- names
        if (sortNames) {
            taxad2 <- taxad2[order(as.numeric(substring(rownames(taxad2), 
                2))), ]
        }
        results[[i]] <- taxad2
        if (plot) {
            taxicDivCont(results[[i]], int.length = 0.2)
            if (nruns > 1) {
                title(paste("Run Num.", i, " of ", nruns, sep = ""))
            }
        }
    }
    if (print.runs) {
        message(paste(nruns, " runs accepted from ", ntries, 
            " total runs (", signif(nruns/ntries, 2), " Acceptance Probability)", 
            sep = ""))
    }
    if (nruns == 1) {
        results <- results[[1]]
    }
    return(results)
}


taxa2cladogram<-function (taxad, drop.cryptic = FALSE, plot = FALSE) 
{
    if (any(taxad[, 6] != taxad[, 1])) {
        for (i in which(taxad[, 1] != taxad[, 6])) {
            taxad[taxad[, 2] == taxad[i, 1], 2] <- taxad[i, 6]
        }
    }
    tlabs <- rownames(taxad)
    desc <- lapply(taxad[, 1], function(x) (taxad[taxad[, 2] == 
        x, 1])[!is.na(taxad[taxad[, 2] == x, 1])])
    ndesc <- sapply(desc, length)
    rank <- numeric(length(ndesc))
    rank[ndesc == 0] <- 1
    rank[rank == 0] <- NA
    while (any(is.na(rank))) {
        rank <- sapply(1:length(rank), function(x) ifelse(!is.na(rank[x]), 
            rank[x], 1 + max(rank[sapply(desc[[x]], function(y) which(y == 
                taxad[, 1]))])))
    }
    comp <- numeric(length(ndesc))
    lab <- list()
    lab[rank == 1] <- tlabs[rank == 1]
    comp[rank == 1] <- 1
    while (any(comp == 0)) {
        tpot <- comp == 0
        tpot2 <- rank == min(rank[tpot])
        tpick <- which(tpot & tpot2)[1]
        dlab <- paste(unlist(lab[desc[[tpick]]]), ",", sep = "", 
            collapse = "")
        lab[[tpick]] <- paste("(", dlab, tlabs[tpick], ")", sep = "")
        comp[tpick] <- 1
    }
    tree1 <- paste(lab[[1]], ";", sep = "")
    tree2 <- read.tree(text = tree1)
    if (drop.cryptic & any(taxad[, 6] != taxad[, 1])) {
        tree2 <- drop.tip(tree2, tlabs[taxad[, 6] != taxad[, 
            1]])
        tree2 <- collapse.singles(tree2)
    }
    if (plot) {
        plot(ladderize(tree2), show.tip.label = FALSE)
    }
    return(tree2)
}