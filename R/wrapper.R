
##Slimmed down blocking for optimal designs
##For use in choiceDes_0.9-3
##Written for [R]-language (tested/built on version 3.5.0)
##Requires: base, stat, AlgDesign
##Revised 2018-06-18

optBlockC <- function(withinData, blocksizes, nRepeats = 5) {

	out <- optBlock(~., withinData, blocksizes, nRepeats=nRepeats)
    return(out)
	
}

##Slimmed down function for generating optimal designs
##For use in choiceDes_0.9-3
##Written for [R]-language (tested/built on version 3.5.0)
##Requires: base, stat, AlgDesign
##Revised 2018-06-15

optFederovC <- function(modelData, nTrials, nRepeats = 5) { 

    out <- optFederov(~., modelData, nTrials, nRepeats=nRepeats)
    return(out)
	
}










