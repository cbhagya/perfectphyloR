#' Extract haplotypes and return it as a list from genotype matrix by vcf file.
#'
#' @param vcfGeno  genotype matrix extract from vcf file.
#' @param snvPosns SNV positions in base paris.
#'
#' @return
#' @export internal
#'
#'
extractHaplos <- function(vcfGeno, snvPosns){
  # Input: vcfGeno = genotype matrix extract from vcf file.
  #  Ex:   HG00096 HG00097 HG00099 HG00100
  #    [1,] "0|0"   "0|0"   "0|0"   "0|0"  
  #    [2,] "0|0"   "0|0"   "0|0"   "0|0"  
  #    [3,] "0|0"   "0|0"   "0|0"   "0|0"  
  #    [4,] "0|1"   "1|1"   "0|0"   "0|0"  
  #    [5,] "0|0"   "0|0"   "0|0"   "0|0"   
  haploMat = matrix(NA, nrow = 2*ncol(vcfGeno), ncol = nrow(vcfGeno))
  
  for(i in 1:ncol(vcfGeno)){
    
    # Extract hap1 and hap2 for individual i.
    # Split "|" or "/". 
    charVec = unlist(strsplit(vcfGeno[,i], "|")) # 1 = indiv1
    # Remove "|" or "/" from the charVec.
    
    numVec = as.numeric(gsub("[^0-1.]", "", charVec))
    numVec = numVec[!is.na(numVec)]
    
    hap1 = numVec[seq(1,length(numVec), by = 2)]
    hap2 = numVec[seq(2,length(numVec), by = 2)]
    
    haploMat[2*i-1, ] = hap1 
    haploMat[2*i, ]   = hap2 
    
  }
  # Remove the columns of haploMat where all zeros and all ones, if there any.
  if(any(colSums(haploMat) == 0 | colSums(haploMat) == 2*nrow(haploMat))){
  # Extract pholymorphic SNVss.
  hapMatpoly = haploMat[, c(-which(colSums(haploMat) == 0),
                           -which(colSums(haploMat) == 2*nrow(haploMat)))]
  
  snvpos = snvPosns[c(-which(colSums(haploMat)==0),
                      -which(colSums(haploMat)==2*nrow(haploMat)))]
  
  rownames(hapMatpoly) <- rep(colnames(vcfGeno), each = 2) 
  colnames(hapMatpoly) <- paste0("SNV",1:ncol(hapMatpoly))
  
  hapData = list(haplomat = hapMatpoly, SNVposns = snvpos)
  
  }else{
    
    rownames(haploMat) <- rep(colnames(vcfGeno), each = 2) 
    colnames(haploMat) <- paste0("SNV",1:ncol(haploMat))
    
    hapData = list(haplomat = haploMat, SNVposns = snvPosns)
  }
  
  return(hapData)
}
