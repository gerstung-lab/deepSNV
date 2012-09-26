# TODO: Add comment
# 
# Author: Moritz Gerstung
###############################################################################

#' List of ambiguous base mappings
#' @noRd
AMBIGUOUS = list()
AMBIGUOUS[["R"]] = c("A","G")
AMBIGUOUS[["Y"]] = c("C","T")
AMBIGUOUS[["M"]] = c("A","C")
AMBIGUOUS[["K"]] = c("G","T")
AMBIGUOUS[["S"]] = c("G","C")
AMBIGUOUS[["W"]] = c("A","T")
AMBIGUOUS[["H"]] = c("A","C","T")
AMBIGUOUS[["B"]] = c("C","G","T")
AMBIGUOUS[["V"]] = c("A","C","G")
AMBIGUOUS[["D"]] = c("A","G","T")
AMBIGUOUS[["N"]] = c("A","C","G","T")
AMBIGUOUS[["A"]] = c("A")
AMBIGUOUS[["T"]] = c("T")
AMBIGUOUS[["C"]] = c("C")
AMBIGUOUS[["G"]] = c("G")
AMBIGUOUS[['-']] = '-'

#' Nice colors
#' 
#' Colors from \code{\link{RColorBrewer}}.
#' @noRd
nt.col <- matrix(c("#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6", "#D9D9D9",  "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A", "#525252"), ncol=2)

#' Complement nucleotide
#' @noRd
COMPLEMENT <- c("A","T","C","G","-")
names(COMPLEMENT) <- c("T","A","G","C","-")