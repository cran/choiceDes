
##Algorithms for producing DCM designs
##For use in choiceDes_0.9-3
##Written for [R]-language (tested/built on version 3.5.0)
##Requires: base, stat, AlgDesign
##Revised 2018-06-18

write.tab <- function(x, f) {
	try(
		write.table(x, f, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
	)
}

dcm.design <- function(cand, nb, sets, alts, fname=NULL, Rd=20, print=TRUE) {
	
	## Initialise
	ta <- Sys.time()
	options(contrasts=c("contr.sum","contr.poly"))
	if(!is.list(cand)) { cand <- list(cand) }
	if(print) {
		cat("dcm.design 0.9-2, ", date(), "\n")
		flush.console()
	}
	
	##ITERATE THROUGH THE CANDIDATE SETS
	drows <- nb * sets * alts
	for(i in 1:length(cand)) {
		##GENERATE FULL-FACTORIAL CANDIDATE SET AND DO EFFECTS CODING
		cand.lev <- gen.factorial(cand[[i]], factors="all")
		cand.eff <- model.matrix(~., cand.lev)
		cand.eff <- cand.eff[,-1]
			
		##LENGTHEN DESIGN SETS, AS REQUIRED, SO THAT THEY CONTAIN
		##AT LEAST AS MANY ROWS AS NEEDED DESIGN
		dd <- cand.lev
		ee <- cand.eff
		while(nrow(cand.lev) < drows) {
			cand.lev <- rbind(cand.lev, dd)
			cand.eff <- rbind(cand.eff, ee)
		}
		rownames(cand.lev) <- 1:nrow(cand.lev)
		rownames(cand.eff) <- 1:nrow(cand.eff)
		
			
		##GENERATE DESIGN
		desa.eff <- optFederovC(cand.eff, drows, nRepeats=Rd)
		desa.lev <- cand.lev[desa.eff$rows,]
			
		rownames(desa.eff$design) <- 1:drows
		rownames(desa.lev) <- 1:drows
			
		##MERGE DESIGNS IN THE CASE OF MULTIPLE CANDIDATE SETS
		if(i == 1) { 
			desb.eff <- desa.eff$design
			desb.lev <- desa.lev
		} else {
			desb.eff <- cbind(desb.eff, desa.eff$design)
			desb.lev <- cbind(desb.lev, desa.lev)
		}
	}
		
	##BLOCK THE FINAL DESIGN
	colnames(desb.eff) <- 1:ncol(desb.eff)
	desx.eff <- optBlockC(desb.eff, rep(sets * alts, nb), nRepeats=Rd)
	desx.lev <- desb.lev[desx.eff$rows,]
	desf.lev <- dcm.design.sort(desx.lev, nb, sets, alts)
	
	##CALCULATE EFFICIENCY, RETURN RESULTS
	d.eff <- dcm.design.effcy(desx.eff$design)
	c.out <- list(levels=desf.lev, effects=desx.eff, d.eff=d.eff)
	if(!is.null(fname)) { write.tab(c.out$levels, fname) }
	
	tb <- Sys.time() - ta
	if(print) {
		cat("parameters: ", ncol(d.eff$V), ", D-efficiency: ", d.eff$D, "\n", sep="")
		cat("time elapsed: ", tb, " secs\n", sep="")
	}
	
	return(c.out)

}

dcm.design.cand <- function(cand, nb, sets, alts, fname=NULL, Rd=20, print=TRUE) {
		
	##THIS FUNCTION TAKES A LEVELS-CODED CANDIDATE SET GENERATED INDEPENDENT
	##FROM THIS FUNCTION (CAN BE RESTRICTED IF DESIRED) AND GENERATES A
	##DESIGN FROM THAT SET
	
	ta = Sys.time()
	options(contrasts=c("contr.sum","contr.poly"))
	if(print) {
		cat("dcm.design.cand 0.9-2, ", date(), "\n")
		flush.console()
	}
	
	cand.lev <- cand
	if(!is.data.frame(cand.lev)) { cand.lev <- read.table(cand.lev, header=TRUE, sep="\t") }
	
	drows = nb * sets * alts
	for(i in 1:ncol(cand.lev)) {
		cand.lev[,i] <- as.factor(cand.lev[,i])
	}
	colnames(cand.lev) <- paste("x",1:ncol(cand.lev),sep="")
		
	cand.eff <- model.matrix(~., cand.lev)
	cand.eff <- cand.eff[,-1]
	
	##LENGTHEN DESIGN SETS, AS REQUIRED, SO THAT THEY CONTAIN
	##AT LEAST AS MANY ROWS AS NEEDED DESIGN
	dd <- cand.lev
	ee <- cand.eff
	while(nrow(cand.lev) < drows) {
		cand.lev <- rbind(cand.lev, dd)
		cand.eff <- rbind(cand.eff, ee)
	}
	rownames(cand.lev) <- 1:nrow(cand.lev)
	rownames(cand.eff) <- 1:nrow(cand.eff)
	
	##GENERATE AND BLOCK DESIGN
	desa.eff <- optFederovC(cand.eff, drows, nRepeats=Rd)
	desx.eff <- optBlockC(desa.eff$design, rep(sets * alts, nb), nRepeats=Rd)
	desx.lev <- cand.lev[desx.eff$rows,]
	desf.lev <- dcm.design.sort(desx.lev, nb, sets, alts)
	
	##CALCULATE EFFICIENCY, RETURN RESULTS
	d.eff <- dcm.design.effcy(desx.eff$design)
	c.out <- list(levels=desf.lev, effects=desx.eff, d.eff=d.eff)
	if(!is.null(fname)) { write.tab(c.out$levels, fname) }
	
	tb <- Sys.time() - ta
	if(print) {
		cat("parameters: ", ncol(d.eff$V), ", D-efficiency: ", d.eff$D, "\n", sep="")
		cat("time elapsed: ", tb, " secs\n", sep="")
	}
	
	return(c.out)	
}

dcm.design.effcy <- function(des) {

	##INTERNAL
	##CALCULATE LINEAR EFFICIENCIES IN AN EFFECTS-CODED DESIGN
	
	des <- as.matrix(des)
	
	M <- t(des) %*% des
	D <- det(M)^(-1/ncol(des))
	V <- solve(M)
	s <- sqrt(diag(V))
	
	return(list(D=D, V=V, s=s))
}


dcm.design.sort <- function(design, nb, sets, alts) {
	
	##INTERNAL
	##ADD CARD, VERSION, TASK VARIABLES TO LEVELS-CODED DESIGN
	##AND DO FINAL SORT
	
	rand <- runif(nb * sets * alts)
	
	card <- rep(1:(nb * sets), each=alts)
	vers <- rep(1:nb, each=(sets * alts))
	task <- rep(rep(1:sets, each=alts), nb)

	return(cbind(card, vers, task, design[order(vers, rand),]))
}

