
##Algorithms for producing MaxDiff and other tradeoff designs
##For use in choiceDes_0.9-3
##Written for [R]-language (tested/built on version 3.5.0)
##Requires: base, stat, AlgDesign
##Revised 2018-06-15

tradeoff.des <- function(items, shown, vers, tasks, fname=NULL, Rd=20, Rc=NULL, print=TRUE) {

	## Initialise
	ta <- Sys.time()
	if(print) {
		cat("tradeoff.des 0.9-2, ", date(), "\n")
		flush.console()
	}

	drows <- vers * tasks
	if(is.null(Rc)) { Rc <- max(1000, 10 * drows) }

	## Find optimal design
	des.c <- tryCatch(
		optBlockC(factor(1:items), rep(shown, drows), nRepeats=Rd),
			error=function(e) return(NULL)
		)
	if(is.null(des.c)) { stop("PROCEDURE STOPPED: Insufficient design space") }
	
	des.m <- t(matrix(des.c$rows, shown, drows))
	
	## Optimize positional balance
	crit <- drows * shown
	crit.vec <- vector("numeric", Rc)
	
	#reorder all rows in design as a starting place
	des.p <- des.m
	for(d in 1:drows) { des.p[d,] <- sample(des.p[d,], shown) }
	
	#reorder one row at a time in candidate
	#and evaluate positional balance
	cc <- 0
	for(i in 1:Rc) {
		ii <- sample(1:drows, 1)
		des.try <- des.p
		des.try[ii,] <- sample(des.try[ii,], shown)
			
		da <- as.vector(unlist(apply(des.try, 2, function(x) table(factor(x, levels=1:items)))))
		if(length(da) == items * shown) {
			crit.try <- abs(sqrt(items * shown) - sum(svd(da)$u))
			if(crit.try < crit) {
				des.p <- des.try
				crit <- crit.try
			}
		}
		
		#count number of iterations that crit has been stable
		#reset counter if crit changes in this iteration
		crit.vec[i] <- crit
		if(i > 1) {
			cc <- cc + 1
			if(crit.vec[i] < crit.vec[i-1]) { cc <- 0 }
		}
	}
	
	## Block reordered design into versions
	des.d <- as.data.frame(des.p)
	for(g in 1:shown) { des.d[,g] <- as.factor(des.d[,g]) }
	
	#error-check block to versions
	des.f <- tryCatch(
		optBlockC(des.d, rep(tasks, vers), nRepeats=Rd), 
			error=function(e) return(NULL), 
			warning=function(e) return(NULL)
		)
	
	if(is.null(des.f)) {
		des.x <- des.p
		warning("Insufficient design space to block into versions; using random assignment instead", 
			call.=FALSE, immediate.=TRUE
		)		
	} else { des.x <- des.p[des.f$rows,] }
	
	
	## Calculate one-way and positional balance
	bal.1w <- table(factor(des.x, levels=1:items))
	bal.1w.mn <- mean(bal.1w)
	bal.1w.sd <- sd(bal.1w)
	
	bal.ps <- apply(des.x, 2, function(x) table(factor(x, levels=1:items)))
	bal.ps.mn <- mean(bal.ps)
	bal.ps.sd <- sd(as.vector(bal.ps))
	
	## Calculate two-way (pairwise) balance
	bal.2w <- pw.eval(items, shown, drows, des.x)
	
	bal.final <- list(
		one.way = list(tbl=bal.1w, mean=bal.1w.mn, stdv=bal.1w.sd),
		col.pos = list(tbl=bal.ps, mean=bal.ps.mn, stdv=bal.ps.sd),
		two.way = bal.2w
	)
	
	Rc.crit = list(Rc=Rc, crit.stable=cc, crit.vec=crit.vec)
	
	## Add card, vers and task columns to design
	card <- 1:drows
	ver <- rep(1:vers, rep(tasks, vers))
	tsk <- rep(1:tasks, vers)
	design <- cbind(card, ver, tsk, des.x)
	colnames(design) <- c("card", "version", "task", paste("item", 1:shown, sep=""))
	
	## Check for zeroes in bal.ps
	zvec <- which(bal.ps == 0)
	if(length(zvec) > 0) {
		warning("SPARSE DESIGN: Not all items appear in all column positions", 
			call.=FALSE, immediate.=TRUE
		)
	}
	
	## Return results
	tb <- Sys.time() - ta
	if(!is.null(fname)) { write.tab(design, fname) }
	if(print) {
		cat("total iterations: ", Rc, "\n", sep="")
		cat("critical D: ", crit, ", stable: ", cc, " iterations\n", sep="")
		cat("time elapsed: ", tb, " secs\n", sep="")
	}
	
	return(list(design=design, balance=bal.final, Rc.crit=Rc.crit, time.elapsed=tb))
}

pw.eval <- function(items, shown, drows, des) {

	pw.mat <- matrix(0, items, items)
	for(i in 1:drows) {
		for(j in 1:shown) {
			a <- des[i,j]
			for(k in 1:shown) {
				b <- des[i,k]
				pw.mat[a, b] <- pw.mat[a, b] + 1
			}
		}
	}
	
	pw.vec <- as.vector(pw.mat)
	pw.vec.od <- pw.vec[-which(pw.vec==diag(pw.mat))]
	pw.mean.od <- mean(pw.vec.od)
	pw.sdev.od <- sd(pw.vec.od)
	
	return(list(tbl=pw.mat, od.mean=pw.mean.od, od.stdv=pw.sdev.od))

}

cp.scree <- function(des) {
	Rc <- des$Rc.crit$Rc
	crit.val <- des$Rc.crit$crit.vec
	
	plot(1:Rc, crit.val, type="l") 
}