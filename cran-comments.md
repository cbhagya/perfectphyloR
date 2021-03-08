## Resubmission March 7, 2021:
  
  * In this submission, a bug in calculating dCor statistics was fixed.
  
  * Updated the version number to 0.2.1 to reflect the bug fix. 
  
  * Test Environments prior to resubmission
    * win-builder (oldrelease)

  * R CMD check results

     0 errors | 0 warnings | 0 note

## Resubmission December 8, 2019:

  * In this resubmission, the following changes were made to the package.
    
    * included a function to create an object of class hapMat
      from VCF file.
    * changed the default value of the argument, minWindow, in the function reconstructPP(), to maximum of one and 2% of the total number of the SNVs.
    * changed the function name reconsPPregion() to reconstructPPregion().
    * added an option to the function reconstructPPregion()
to reconstruct perfect phylogenies between a given range of SNV positions. 
    * added some warning messages to the function, reconstructPPregion(), to warn about the timing, and to show the current focal SNV of the hapMat object. 
    * made the packages HHG and dendextend optional to install with the package.
    * added the bioRxiv preprint of the package as the PDF vignette.
    * updated the DESCRIPTION file:
        + updated the version number to 0.2.0 to reflect the new function added, and the modifications of reconstructPP() and reconstructPPregion(). 
        + included R.rsp as VignetteBuilder.
        + made the packages HHG, dendextend, vcfR, and R.rsp as suggests.
      
  * Test Environments prior to resubmission
    * local Windows OS install, R 3.4.2
    * local Linux, R 3.5.0
    * win-builder (oldrelease)


  * R CMD check results

     0 errors | 0 warnings | 1 note

      checking CRAN incoming feasibility ... NOTE
      Maintainer: 'Charith Karunarathna <ckarunar@sfu.ca>'
      
    This submission comes actually from this maintainer, Charith Karunarathna, and not anybody else.   

## Resubmission June 2019:

### Third resubmission in June 2019

  * In this submission according to a reviewer's suggestion, I have:
  
      * fixed the 'Description' field in the DESCRIPTION start with "Reconstructs perfect ...".
      * fixed the reference in the 'Description' field of the DESCRIPTION file to be of the form: authors (year) <arXiv:...>.
      * added more small executable examples in the Rd-files to illustrate the use of the exported functions (internal functions) and enabled those examples for automatic testing.
  
  * I also updated the version number in the DESCRIPTION file to 0.1.3 to reflect the above modification.
  
### Second resubmission in June 2019

  * In this submission, vignette files were modified to comply with CRAN requirements.
  
  * Updated the version number in the DESCRIPTION file to 0.1.2 to reflect the modification in vignette folder.

### First resubmission in June 2019

  * Fixed a note related to SHLIB_OPENMP_*FLAGS in Makefiles following URL:        https://stat.ethz.ch/pipermail/r-package-devel/2019q2/003764.html.  I thank Iñaki Ucar and Dirk Eddelbuettel
for proving me the guidances.

  * Updated the version number in the DESCRIPTION file to Version: 0.1.1



## Test environments
* local Windows install, R 3.4.2
* local Linux, R 3.5.0
* win-builder (oldrelease)

## R CMD check results

0 errors | 0 warnings | 1 note

Possibly mis-spelled words in DESCRIPTION:
  Charith (8:167)
  Jinko (8:194)
  Karunarathna (8:177)
  Phylogenies (3:28)

 * The words Charith, Jinko and Karunarathna are not misspelled. These are the names of authors.
 * The word Phylogenies is not misspelled.

## Downstream Dependencies

There are currently no downstream dependencies for this package.
