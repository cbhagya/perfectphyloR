
#' Create a \code{hapMat} object from variant call format (\code{vcf}) file.
#' 
#' This function creates a \code{hapMat} object from variant call format (\code{vcf}) file.   
#'
#' @param vcf_file_path File path to the \code{vcf} file.
#'
#' @return An object of class \code{hapMat}.
#' @export
#'
#' @examples
#' 
#' # Specify the file path.
#' vcf_file_path <- "C:/vcfData/vcfData.vcf.gz"
#' # Create a hapMat object from the vcf file.
#' ex_vcf_hapMat <- vcftohapMat(vcf_file_path) 
#'  
#' 
vcftohapMat <- function(vcf_file_path){
   
  if (!requireNamespace("vcfR", quietly = TRUE)) {
    stop("Package \"vcfR\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  # Import vcf file 
  vcf_data <- vcfR::read.vcfR(vcf_file_path)
  # 1. SNV positions
  SNV_pos = vcf_data@fix[,2]
  
  
  #######################
  # 2. Extract individual IDs
  
  # Remove the variable named "FORMAT" and store the resulting data as genos
  genos <- vcf_data@gt[, -1] # Genotype matrix
  
  indID = colnames(genos)
  # 3. Utility function extractHaplos
  message('Generating hapMat object...')
  haploData = extractHaplos(vcfGeno = genos, snvPosns = SNV_pos)
  
  hapMat = perfectphyloR::createHapMat(hapmat = haploData$haplomat,
                                       snvNames = colnames(haploData$haplomat),
                                       hapNames = rownames(haploData$haplomat),
                                       posns = haploData$SNVposns)
  return(hapMat)
}



