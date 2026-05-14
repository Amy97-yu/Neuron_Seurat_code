hsa2mmu <- function(gl) {
  if (!require(homologene)) {
    install.packages("homologene")
    require(homologene)
  }
  genes <- homologene::homologene(gl,inTax = 9606, outTax = 10090)
  gl <- genes$`10090`
  return(gl)
}
