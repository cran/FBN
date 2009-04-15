##rm(list = ls())
#
#setwd("D:/Adi/R/FBNormalization")
#memory.limit(size=4000)
#
#rawDataFileName = "CBS_25HMLCs_12022008.txt"
#fishProbesFileName = "FISHprobes_and_FISHcns_25HMLCs_4probes.txt"
#normDataFileName= "CBS_25HMLCs_FBN_13112008.txt"

#rawDataFileName = "CBS_45MM_11022008.txt"
#fishProbesFileName = "FISHprobes_and_FISHcns_45MM_30092008.txt"
#normDataFileName= "CBS_45MM_FBN_13112008.txt"

#############################################################################################
#                      FBN FUNCTION                      									#
#                                                        									#
# FBNormalization(rawDataFileName, fishProbesFileName, normDataFileName, debugFlag = FALSE) #
#############################################################################################
FBNormalization <- function (rawDataFileName = NULL, fishProbesFileName = NULL, normDataFileName = NULL, debugFlag = FALSE)
{
	k = list.files()
	########################################
	# check the existence of the input files
	if((length(which(k == rawDataFileName)) == 0) || is.null(rawDataFileName)){
		cat("ERROR - can not find rawDataFileName: ", rawDataFileName, "\n")
		return(NULL)	
	}
 	if((length(which(k == fishProbesFileName)) == 0) || is.null(fishProbesFileName)){
		cat("ERROR - can not find fishProbesFileName:" ,  fishProbesFileName, "\n")
		return(NULL)
	}
	##############################################################################
	# check if the output file exists or not, and ask if necessary for overwriting
	while((length(which(k == normDataFileName)) != 0) || is.null(normDataFileName) || (normDataFileName == "")) {
		if(is.null(normDataFileName) || (normDataFileName == "")){
			cat("WARNING: normDataFileName is NULL.\n")
			r = "n"
		}
		else{
			cat("WARNING: ", normDataFileName, " already exists. Overwrite Yes/No/Stop?(y/n/s):")
			r = ""
			while((r != "y") && (r != "n")&& (r != "s")) {
				r = readline()
			}
		}
		if(r == "s") return()
		if(r != "y") {
			cat("Please input a new output file name: ")
			normDataFileName = readline()
		}
		if(r == "y") break
	}	
	rm(k)
	if(length(grep(".txt", normDataFileName)) == 0)
		normDataFileName = paste(normDataFileName, ".txt", sep = "")
	######################
	# read the input files
	rawData = read.table(rawDataFileName,header=TRUE,sep="\t")
	fishProbe = read.table(fishProbesFileName,header=TRUE,sep="\t")
	######################
    # copy raw data fields
    chrData = rawData$chrom
    nameData = names(rawData)[5:dim(rawData)[2]]
    posData = rawData$physical.position
    cytoData = as.character(rawData$cytoband)
    rawData = as.matrix(rawData[,-(1:4)])
    rowData = dim(rawData)[1]
    colData = dim(rawData)[2]
    #######################
    # copy FISH data fields
    nameFish = names(fishProbe)[6:dim(fishProbe)[2]]
    chrFish = fishProbe$chromosome
    cytoFish = fishProbe$cytoband
    cloneFish = fishProbe$BACclone
    startFish = fishProbe$start.loc
    endFish = fishProbe$end.loc
    fishProbe = as.matrix(fishProbe[,-(1:5)])
    rowFish = dim(fishProbe)[1]
    colFish = dim(fishProbe)[2]
    #########################
    # check the files headers
    flagError = FALSE
    if(is.null(chrData)){
    	flagError = TRUE
    	cat("Please name \"chrom\" the field containing the chromosomes index in the rawData file!\n")
    }
    if(is.null(posData)){
    	flagError = TRUE
    	cat("Please name \"physical.position\" the field containing the physical position of the SNP probes in the rawData file!\n")
    }
    if(is.null(cytoData)){
    	flagError = TRUE
    	cat("Please name \"cytoband\" the field containing the cytoband location of the SNP probes in the rawData file!\n")
    }
    if(is.null(chrFish)){
    	flagError = TRUE
    	cat("Please name \"chromosome\" the field containing the chromosomes index in the FISH probes file!\n")
    }
    if(is.null(cytoFish)){
    	flagError = TRUE
    	cat("Please name \"cytoband\" the field containing the cytoband location of the FISH probes in the FISH probes file!\n")
    }
    if(is.null(cloneFish)){
    	flagError = TRUE
    	cat("Please name \"BACclone\" the field containing the clone name of the FISH probes in the FISH probes file!\n")
    }
    if(is.null(startFish)){
    	flagError = TRUE
    	cat("Please name \"start.loc\" the field containing the start location of the FISH probes in the FISH probes file!\n")
    }
    if(is.null(endFish)){
    	flagError = TRUE
    	cat("Please name \"end.loc\" the field containing the start location of the FISH probes in the FISH probes file!\n")
    }
	if(flagError) return()
    ##########################
    # define the output matrix
	normData = rawData
    normData[,5:dim(normData)[2]] = NA; 
    ###############################################
    # test if all data have a correspondent in FISH
    # and eliminate the possible additional entries
    flagHasFish = vector(mode = "numeric", length = colData)
    newFishProbes = matrix(data = NA, nrow = rowFish, ncol = colData)
    newNameFish = vector(mode = mode(nameFish), length = colData)
    for(i in 1:colData){
        j = which(nameFish == nameData[i])
        if(length(j) >= 1){
            flagHasFish[i] = 1
            newFishProbes[,i] = fishProbe[, j]
            newNameFish[i] = nameFish[j]
        } else
            newNameFish[i] = NA
    }
    if(length(flagHasFish[flagHasFish == 0]) > 0){
        cat("WARNING: There is data with no FISH associated info: Median Centering Normalization will be used\n")
    }
    fishProbe = newFishProbes
    nameFish = newNameFish
    rm(newFishProbes, newNameFish) 
    colFish = dim(fishProbe)[2] 
    #########################################
    # Associate the CN read by the FISH probe
    # to the corresponding raw SNP value     
    snpProbe = matrix(data = NA, nrow = rowFish, ncol = colData)
    for( i in 1:rowFish){
        posProbe = posData
        posProbe[chrData != chrFish[i]] = 0
        idxStartPosProbe = which(posProbe > startFish[i])
        idxEndPosProbe = which( (posProbe < endFish[i])&(posProbe != 0) )
        if( (length(idxStartPosProbe) == 0)|(length(idxEndPosProbe) == 0) ){
            cat("WARNING: Fish probe ", cloneFish[i], " has no physical correspondence with the data...\n" )
            rm(posProbe, idxStartPosProbe, idxEndPosProbe)
            next
        }
        if(idxStartPosProbe[1] > idxEndPosProbe[length(idxEndPosProbe)])
            idxPosProbe = idxEndPosProbe[length(idxEndPosProbe)]:idxStartPosProbe[1]
        else
            idxPosProbe = idxStartPosProbe[1]:idxEndPosProbe[length(idxEndPosProbe)]
        
        allSnpProbe = rawData[idxPosProbe, ]
        if(length(idxPosProbe) == 1)
        	snpProbe[i,] = allSnpProbe
        else
	        for( j in 1:colData)
    	        snpProbe[i,j] = median(allSnpProbe[,j])
        rm(posProbe, idxStartPosProbe, idxEndPosProbe, idxPosProbe)
    }
    ###################################################################
    # apply k-means on each probe to determine the normalization values
    histRawData = hist(as.vector(as.matrix(rawData)), plot = FALSE, breaks = "FD")
	breaksRawData = histRawData$breaks
    normalizingSNP = matrix(data = NA, nrow = rowFish, ncol = colData)
    for(i in 1:colData){
		if(flagHasFish[i] == 0) next
		lineData = rawData[chrData != chrData[length(chrData)], i]
		clusterLineData = FBN.kmeans(inputData = lineData, minSpan = 0.1, breaksData = breaksRawData)
		for(j in 1:rowFish){
			# find the cluster that contains the snpProbe value #
			ind = clusterLineData$cluster[lineData == snpProbe[j,i]]
			normalizingSNP[j, i] = median(lineData[clusterLineData$cluster == ind[1] ] )
		}
		rm(lineData, clusterLineData)
    }
    rm(histRawData)
    ###################################
    # Check the FISH probes consistency
    for(i in 1:colData){
    	if(flagHasFish[i] == 0) next
    	cn = fishProbe[,i]
    	snp = normalizingSNP[!is.na(cn),i]
    	probe = snpProbe[!is.na(cn),i]
    	cn = cn[!is.na(cn)]
		sortSnp = sort(snp, index.return = TRUE)
		sortCN = cn[sortSnp$ix]
		sortProbe = probe[sortSnp$ix]
		# verify if CN is in ascending order
		diffSnp = sortSnp$x[2:length(sortSnp$x)] - sortSnp$x[1:(length(sortSnp$x)-1)]
		diffCN = sortCN[2:length(sortCN)] - sortCN[1:(length(sortCN)-1)]
		diffSnp[diffSnp>0] = 1
		diffSnp[diffSnp<0] = -1
		diffCN[diffCN>0] = 1
		diffCN[diffCN<0] = -1
		flagInconsistency = FALSE
		for(j in 1:length(diffCN)){
			if(diffSnp[j] != diffCN[j])
				flagInconsistency = TRUE
		}
		if(flagInconsistency){
			flagHasFish[i] = 2
			cat("WARNING: Found inconsistent FISH info for probe: ", nameData[i], "\n")
			cat("NormalizingSNP -> fishCN -> rawSNP\n")
	    	for(j in 1:length(sortCN))
				cat(sortSnp$x[j], " -> ", sortCN[j], " -> ", sortProbe[j], "\n")
		}
	}
    
    #################################
    # Start the normalization process
    #################################
    #make three stages normalization
    # first all probes with good fish data -> flagHasFish == 1
    # second all probes with inconsistent data -> flagHasFish == 2
    # third all probes with no fish data - median centering -> flagHasFish == 0

	########################################################
    # set up the order of the CN iterations
    orderCN = c(2, 1, 3:6)
    valueCN = vector(mode = "numeric", length = length(orderCN)) #initialize the nominal values vector of the CNs
    valueCN[2] = 2
    thresholdsCN = vector(mode = "numeric", length = 6) # thresholds = {CN0-1, CN1-2, CN2-3, CN3-4, CN4-5, CN5-more}
    flagIsNormalized = vector(mode = "numeric", length = colData)
    ##############################################################################################################
    # First stage!
    # All probes with good fish data -> flagHasFish == 1
    # Normalize on the closest fishCN == 2
    for(j in 1:length(orderCN) ){
        CN = orderCN[j]
        flagNewData = FALSE
		for(i in 1:colData){
			if((flagIsNormalized[i] == 1)||(flagHasFish[i] != 1)) next
			index = which((fishProbe[,i] - CN) == 0) #look for a normalizing fish-snp pair
			if(length(index) == 0) next # current probe does not have a normalizing fish-snp pair
	    	cat("Normalize probe: ", nameData[i], "\n")
			flagNewData = TRUE
			index = index[1]
			normData[,4+i] = FBN.valueCenter(inputData = rawData[,i], normalizingValue = normalizingSNP[index,i], nominalValueCN = valueCN[CN])
			flagIsNormalized[i] = 1
			if(debugFlag){
				par(mfrow = c(2, 1))
				hist(rawData[,i], plot = TRUE, breaks = breaksRawData, main = nameData[i])
				hist(normData[,4+i], plot = TRUE, breaks = "FD")
				#dev.print(device = bmp, width=1024, height=768, paste(nameData[i],".bmp",sep="")) 
				cat("push enter to continue....")
				readline()
			}
		}
		if(!flagNewData) next
		cat("\nK-MEANS for group CN", CN, "\n")
        allNormData = as.matrix(normData[chrData != chrData[length(chrData)], 4+which(flagIsNormalized == 1)])
        allNormData = as.numeric(as.vector(allNormData))
        if(debugFlag){
        	par(mfrow = c(1, 1))
        	histAllNormData = hist(allNormData, plot = TRUE, breaks = "FD")
        }
        #set the minSpan to 0.3 for the final clustering...
        clusterAllNormData = FBN.kmeans(inputData = allNormData, minSpan = 0.3, breaksData = NULL)
		#plot(allNormData, col = clusterAllNormData$cluster)
		tempNominalCN = abs(clusterAllNormData$centers - 2)
		idxCN2 = which(tempNominalCN == min(tempNominalCN))
		# merge all clusters that are inferior to CN2 (due to possible partial hybridizations)
		clusterAllNormData$cluster[clusterAllNormData$cluster < idxCN2] = idxCN2 - 1
		indexCNs = c()
		for(i in 1:6){
			if((idxCN2-2+i) < 1) next
			if((idxCN2-2+i) > length(clusterAllNormData$centers))
				break
			indexCNs = c(indexCNs, idxCN2-2+i)
			#valueCN[i] = clusterAllNormData$centers[idxCN2-2+i]
			valueCN[i] = median(allNormData[clusterAllNormData$cluster == indexCNs[i]])
			if(i == 1)
				thresholdsCN[i] = min(allNormData[clusterAllNormData$cluster == (idxCN2-2+i)])
			else
				thresholdsCN[i] = (min(allNormData[clusterAllNormData$cluster == (idxCN2-2+i)]) + max(allNormData[clusterAllNormData$cluster == (idxCN2-2+i-1)]))/2
		}
		if(length(which(flagIsNormalized[flagHasFish == 1] == 0)) > 0)
			rm(allNormData, clusterAllNormData)
		if(debugFlag)
			rm(histAllNormData)
		cat("\nCN values: ", valueCN, "\n\n")
        if(debugFlag){
        	cat("push enter to continue....")
			readline()
		}
    }
    ##############################################################################################################
    # Second stage!
    # Normalize on each available fish data and verify all others. 
    # Keep the normalization with the highest success score
    for(i in 1:colData){
    	if((flagIsNormalized[i] == 1)||(flagHasFish[i] != 2)) next
    	score = vector(mode = "numeric", length = rowFish)
    	cat("Normalize probe: ", nameData[i], "\n")
    	#normalize on each fish info
    	for(j in 1:rowFish){
    		if(is.na(fishProbe[j, i])) next
			# need to normalize only the normalizingSNP to calculate the score. Then normalize the data for a maximum score...
			tmpNormSNP = FBN.valueCenter(inputData = normalizingSNP[,i], normalizingValue = normalizingSNP[j,i], nominalValueCN = valueCN[fishProbe[j,i]])
			#calculate score for current normalization
			tmpCN = vector(mode = "numeric", length = rowFish)
			score[j] = 0
			for(k in 1:rowFish){
				if(is.na(fishProbe[k, i])) next
				tmp = thresholdsCN - tmpNormSNP[k]
				idx1 = which(tmp <= 0)
				idx2 = which(tmp >= 0)
				if(length(idx2)>0)
					tmpCN[k] = idx2[1]-1
				else 
					tmpCN[k] = 5
				if(tmpCN[k] == fishProbe[k, i])
					score[j] = score[j] + 1
			}
		}
		idxNorm = which(score == max(score))
		cat(" with success scores", score, "\n")
		normData[,4+i] = FBN.valueCenter(inputData = rawData[,i], normalizingValue = normalizingSNP[idxNorm[1],i], nominalValueCN = valueCN[fishProbe[idxNorm[1],i]])
		flagIsNormalized[i] = 1
		if(debugFlag){
			par(mfrow = c(2, 1))
			hist(rawData[,i], plot = TRUE, breaks = breaksRawData, main = nameData[i])
			hist(normData[,4+i], plot = TRUE, breaks = "FD")
			#dev.print(device = bmp, width=1024, height=768, paste(nameData[i],".bmp",sep="")) 
			cat("push enter to continue....")
			readline()
		}
    }
    ##############################################################################################################
    # Third stage!
    # Classical Median Centering Normalization - The median value of the raw data is normalized in CN2
    for(i in 1:colData){
    	if((flagIsNormalized[i] == 1)||(flagHasFish[i] != 0)) next
    	cat("Normalize probe: ", nameData[i], "\n")
    	normData[,4+i] = FBN.valueCenter(inputData = rawData[,i], normalizingValue = median(rawData[chrData != chrData[length(chrData)],i]), nominalValueCN = 2)
		flagIsNormalized[i] = 1
    }
    ############################################################
    # Make the final clusterization and determine the thresholds
	cat("\nFinal K-MEANS for all data", "\n")
    allNormData = as.matrix(normData[chrData != chrData[length(chrData)], 4+which(flagIsNormalized == 1)])
    allNormData = as.numeric(as.vector(allNormData))
	#if(debugFlag){
    	par(mfrow = c(1, 1))
    	histAllNormData = hist(allNormData, plot = TRUE, breaks = "FD")
    #}
	clusterAllNormData = FBN.kmeans(inputData = allNormData, minSpan = 0.3, breaksData = NULL)
	#plot(allNormData, col = clusterAllNormData$cluster)
	tempNominalCN = abs(clusterAllNormData$centers - 2)
	idxCN2 = which(tempNominalCN == min(tempNominalCN))
	######################################################################################
	# merge all clusters that are inferior to CN2 (due to possible partial hybridizations)
	clusterAllNormData$cluster[clusterAllNormData$cluster < idxCN2] = idxCN2 - 1
	indexCNs = c()
	for(i in 1:6){
		if((idxCN2-2+i) < 1) next
		if((idxCN2-2+i) > length(clusterAllNormData$centers))
			break
		indexCNs = c(indexCNs, idxCN2-2+i)
		valueCN[i] = median(allNormData[clusterAllNormData$cluster == indexCNs[i]])
		if(i == 1)
			thresholdsCN[i] = min(allNormData[clusterAllNormData$cluster == (idxCN2-2+i)])
		else
			thresholdsCN[i] = (min(allNormData[clusterAllNormData$cluster == (idxCN2-2+i)]) + max(allNormData[clusterAllNormData$cluster == (idxCN2-2+i-1)]))/2
	}
	cat("\nCN values: ", valueCN, "\n\n")
    ############################################################
    # determine the thresholds
    # fit gaussian on each cluster and find their meeting points
    if(length(indexCNs) == 6)
    	indexCNs = indexCNs[1:5]
    meanCNs = vector(mode = "numeric", length = length(indexCNs))
    stdCNs = vector(mode = "numeric", length = length(indexCNs))
    for(i in 1:length(indexCNs)){
	    meanCNs[i] = mean(allNormData[clusterAllNormData$cluster == indexCNs[i]])
		stdCNs[i] = sd(allNormData[clusterAllNormData$cluster == indexCNs[i]])
    }
    #################################################################
    # thresholds = {CN0-1, CN1-2, CN2-3, CN3-4, CN4-5, CN5-6 or more}
    thresholds = vector(mode = "numeric", length = (length(indexCNs)+1))
    thresholds[1] = meanCNs[1] - 2*stdCNs[1]
    thresholds[length(indexCNs)+1] = meanCNs[length(indexCNs)] + 2*stdCNs[length(indexCNs)]
    for(i in 2:(length(indexCNs))){
    	thresholds[i] = meanCNs[i] - (meanCNs[i] - meanCNs[i-1])*stdCNs[i]/(stdCNs[i] + stdCNs[i-1])
    }
    for(i in 1:length(thresholds)){
	    cat("threshold CN", i-1, " to CN", i, " = ", thresholds[i], "\n")
	}
	write.table(normData, file = normDataFileName, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	thresholdsFileName = strsplit(normDataFileName, ".txt");
	thresholdsFileName = paste(thresholdsFileName, "_thresholds.txt", sep = "")
	b = c("CN1", "CN2", "CN3", "CN4", "CN5", "")
	b = c(b[1:length(indexCNs)], "")
	a = valueCN
	a = c(a[1:length(indexCNs)], "")

	bb = c("ThrCN0-1", "ThrCN1-2", "ThrCN2-3", "ThrCN3-4", "ThrCN4-5", "ThrCN5-6")
	bb = bb[1:(length(indexCNs)+1)]
	bb[length(bb)] = paste(bb[length(bb)], " or more", sep = "")
	aa = thresholds
	aa = aa[1:(length(indexCNs)+1)]
	outInfo = rbind(b, a, bb, aa)
	
	write.table(outInfo, file = thresholdsFileName, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
	cat("FBNormalization done!\n")
	return()
}

########################################################################################################
########################################################################################################

#############################################################################
#  FBN.valueCenter															#
#																			#
# FBN.valueCenter(inputData, normalizingValue, nominalValueCN, logScale)	#
#############################################################################
FBN.valueCenter <- function(inputData = NULL, normalizingValue = NULL, nominalValueCN = 2, logScale = FALSE){	
	if(is.null(inputData)){
        cat("WARNING: FBN.valueCenter -> Please input a valid inputData\n") 
        return(NULL)
    }
    if(is.null(normalizingValue) ){
        cat("WARNING: FBN.valueCenter -> Please input a valid normalizingValue\n") 
        return(NULL)
    }
    if(logScale)
    	normalizedData = inputData + log2(nominalValueCN) - log2(normalizingValue) 
    else
		normalizedData = inputData * nominalValueCN / normalizingValue
	return(normalizedData)
}


########################################################################################################
########################################################################################################

#########################################################################################
#		FBN k-means 																	#
#																						#
# FBN.kmeans(inputData = NULL, minSpan = 0.1, breaksData = NULL)						#
#########################################################################################
FBN.kmeans <- function(inputData = NULL, minSpan = 0.1, breaksData = NULL){
    if(is.null(inputData) ){
        cat("WARNING: FBN.kmeans -> Please input a valid inputData\n") 
        return(NULL)
    }
	maximas = FBN.histogramMaxima(inputData, minSpan = minSpan, breaksData = breaksData)
    newIteration = TRUE
    while(newIteration){
        newIteration = FALSE
        noClusters = length(maximas)
        if(noClusters == 1){
        	clusterData = c()
        	clusterData$cluster = vector(mode = "numeric", length = length(inputData)) + 1
        	clusterData$centers = median(inputData)
        	clusterData$size = length(inputData)
        	clusterData$withinss = NA
        	break
        }
        if(noClusters ==2)
            clusterData = kmeans(inputData, noClusters, iter.max = 10, nstart = 25, algorithm = "Hartigan-Wong")
        else
            clusterData = kmeans(inputData, maximas, iter.max = 100, nstart = 1, algorithm = "Hartigan-Wong")
        # put the clusters in ascending order of their centers
        ordCenters = sort(clusterData$centers, index.return = TRUE)
        difIDX = ordCenters$ix[2:length(ordCenters$ix)] - ordCenters$ix[1:(length(ordCenters$ix)-1)]
        if(length(which((ordCenters$ix[2:length(ordCenters$ix)] - ordCenters$ix[1:(length(ordCenters$ix)-1)]) != 1)) > 0){
        	#reorder
        	copyClusterData = clusterData
        	clusterData$centers = copyClusterData$centers[ordCenters$ix]
        	clusterData$size = copyClusterData$size[ordCenters$ix]
        	clusterData$withinss = copyClusterData$withinss[ordCenters$ix]
        	for(i in 1:length(ordCenters$ix))
				clusterData$cluster[copyClusterData$cluster == i] = ordCenters$ix[i]
			rm(copyClusterData)
        }
        
        # the percentage of data contained by the clusters
        percentageClusters = 100*clusterData$size/sum(clusterData$size)
        #remove clusters with less than 1%
        idxRemovableClusters = which(percentageClusters <= 1)
        idxRemainingClusters = which(percentageClusters > 1)
        if(length(idxRemovableClusters) > 0){
        	newIteration = TRUE
        	maximas = clusterData$centers[idxRemainingClusters]
        	next
        }
		#chech if the centers of the found clusters are not to close to eachother - namely 0.2
		maximas = clusterData$centers[idxRemainingClusters]
		distClusters = maximas[2:length(maximas)] - maximas[1:(length(maximas) - 1)]
		minDistance = 0.2
		idxFaultyClusters = which(distClusters < minDistance)
		if(length(idxFaultyClusters)>0){
			newIteration = TRUE
			# remove the cluster with a smaller percentage
			if(percentageClusters[idxFaultyClusters[1]] < percentageClusters[idxFaultyClusters[1]+1])
				eliminateCluster = idxFaultyClusters[1]
			else
				eliminateCluster = idxFaultyClusters[1]+1
			remainingMaximas = c()
			for(i in 1:length(maximas))
				if(i != eliminateCluster)
					remainingMaximas = c(remainingMaximas, maximas[i])
			maximas = remainingMaximas
		}
    }
	centers = vector(mode = "numeric", length = length(clusterData$centers) )
	for(i in 1:length(clusterData$centers) ){
		centers[i] = median(inputData[clusterData$cluster == i])
	}
	clusterData$centers = centers
	return(clusterData)
}

########################################################################################################
########################################################################################################

#########################################################################################
#     FBN.histogramMaxima FUNCTION   	  												#
#            	                     													#
# FBN.histogramMaxima(inputData = NULL, minSpan = 0.1, breaksData = NULL)				#
#########################################################################################
FBN.histogramMaxima <- function(inputData, minSpan = .1, breaksData = NULL){
    if(is.null(inputData) ){
        cat("WARNING: hist.profile -> Please input a valid inputData\n") 
        return(NULL)
    }
    maximas = c()
    idxMaximas = c()
    if(is.null(breaksData) ){
	    histData = hist(inputData, plot = FALSE, breaks = "FD")
	} else
	    histData = hist(inputData, plot = FALSE, breaks = breaksData)
#	if(length(histData$counts) > length(inputData)){
#		step = (max(inputData) - min(inputData))/length(inputData);
#		breaksData = min(inputData) - step
#		while(max(breaksData) <= max(inputData))
#			breaksData = c(breaksData, max(breaksData)+step)
#		histData = hist(inputData, plot = FALSE, breaks = breaksData)
#	}
    if(minSpan < 0){
	    filteredHistData = medianFilter(inputData = histData$counts, windowSize = 5)
    	filteredHistData = meanFilter(inputData = filteredHistData, windowSize = 9)
    	histData$counts = filteredHistData
    	rm(filteredHistData)
    }
    inputValues = histData$mids
    deltaHist = histData$counts[2:length(histData$counts)] - histData$counts[1:(length(histData$counts)-1)]
    deltaHist[deltaHist < 0] = -1
    deltaHist[deltaHist > 0] = 1
    i = 1
    while(i < length(deltaHist) ){
        if(deltaHist[i] != 1){ #constant or descending region
            i = i+1
            next
        }
        j = i+1
        if(deltaHist[j] == 1){ #an ascending region
            i = j
            next
        }
        #if here than it might be a clear maximum or a flat maximum region...
        while( (deltaHist[j] == 0)&&(j < length(deltaHist) ) ){
            j = j+1
        }
        if(deltaHist[j] == -1){ #clear maximum found /\
            maximas = c(maximas, inputValues[j])
            idxMaximas = c(idxMaximas, j)
            #maximas = c(maximas, inputValues[i+floor( (j-i)/2)])
            #idxMaximas = c(idxMaximas, i+floor( (j-i)/2)+2)
        }
        i = j
    }
    outlier = 1
    while((length(maximas) >= 2)&&(outlier != 0)){
    	outlier = 0
		deltaMaximas = maximas[2:length(maximas)] - maximas[1:(length(maximas)-1)]
		idx = which(deltaMaximas < abs(minSpan))
		if(length(idx) > 0){
			if(histData$counts[idxMaximas[idx[1]]] < histData$counts[idxMaximas[idx[1]+1]])
				outlier = idx[1]
			else
				outlier = idx[1]+1
			if(outlier == 1){
				maximas = maximas[2:length(maximas)]
				idxMaximas = idxMaximas[2:length(idxMaximas)]
			}
			else{
				if(outlier == length(maximas)){
					maximas = maximas[1:(length(maximas)-1)]
					idxMaximas = idxMaximas[1:(length(idxMaximas)-1)]
				}
				else{
					maximas = c(maximas[1:(outlier-1)], maximas[(outlier+1):length(maximas)])
					idxMaximas = c(idxMaximas[1:(outlier-1)], idxMaximas[(outlier+1):length(idxMaximas)])
				}
			}
		}
	}
	backMaximas = maximas	
	#Check if they are too many maximas -> apply the smoothing filter to reduce the noise...
	if((length(maximas) > 6)&&(minSpan >= 0)){
	    maximas = FBN.histogramMaxima(inputData, minSpan = -0.1, breaksData = breaksData)
	}
	if(length(maximas)==0)
		return(backMaximas)
	#for(i in 1:length(maximas)){
	#	tmp = abs(inputData - maximas[i])
	#	idx = which(tmp == min(tmp))
	#	maximas[i] = inputData[idx[1]]
	#}
    return(maximas)
}


########################################################################################################
########################################################################################################

##################################################
#     medianFilter FUNCTION                     #
#                                                #
# medianFilter(inputData, windowSize)           #
##################################################
medianFilter <- function(inputData=NULL, windowSize=3){
    if(is.null(inputData) ){
        cat("WARNING: medianFilter -> Please input a valid inputData\n") 
        return(NULL)
    }
    if(windowSize <= 1){
        cat("WARNING: medianFilter -> no filtering performed: outData = inputData\n")
        return(inputData)
    }
    outData = vector(mode = mode(inputData), length = length(inputData) )
    for(i in 1:length(inputData) ){
        if(windowSize/2 == floor(windowSize/2) ){
            startTmp = max(1, i-windowSize/2)
            endTmp = min(length(inputData), i+windowSize/2-1)
        } else{
            startTmp = max(1, i-(windowSize-1)/2)
            endTmp = min(length(inputData), i+(windowSize-1)/2)
        }
        outData[i] = median(inputData[startTmp:endTmp])
    }
    return(outData)
}

########################################################################################################
########################################################################################################


############################################
#     meanFilter FUNCTION                 #
#                                          #
# meanFilter(inputData, windowSize)       #
############################################
meanFilter <- function(inputData=NULL, windowSize=3){
    if(is.null(inputData) ){
        cat("WARNING: meanFilter -> Please input a valid inputData\n") 
        return(NULL)
    }
    if(windowSize <= 1){
        cat("WARNING: meanFilter -> no filtering performed: outData = inputData\n")
        return(inputData)
    }
    outData = vector(mode = mode(inputData), length = length(inputData) )
    for(i in 1:length(inputData) ){
        if(windowSize/2 == floor(windowSize/2) ){
            startTmp = max(1, i-windowSize/2)
            endTmp = min(length(inputData), i+windowSize/2-1)
        } else{
            startTmp = max(1, i-(windowSize-1)/2)
            endTmp = min(length(inputData), i+(windowSize-1)/2)
        }
        outData[i] = mean(inputData[startTmp:endTmp])
    }
    return(outData)
}

########################################################################################################
########################################################################################################
