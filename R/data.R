#' demo data for shine
#'
#' A dataset containing the name mz rt and sample intensity:
#'
#' @format A data frame with 1351 rows and 181 variables:
#' \describe{
#'   \item{name}{compound name or feature name}
#'   \item{mz}{mz}
#'   \item{rt}{rt}
#'   \item{isotope}{isotope}
#'   \item{QC1}{QC}
#'   \item{QC2}{QC}
#'   \item{QC3}{QC}
#'   \item{QC4}{QC}
#'   \item{QC5}{QC}
#'   \item{S20}{S}
#'   \item{S30}{S}
#'   \item{S39}{S}
#'   \item{S49}{S}
#'   \item{S59}{S}
#'   \item{S69}{S}
#'   \item{S78}{S}
#'   \item{S2}{S}
#'   \item{S3}{S}
#'   \item{S4}{S}
#' }
"data"
#' demo data for shine
#'
#' A dataset containing the sample.name injection.order batch,group and class:
#'
#' @format A data frame with 177 rows and 5 variables:
#' \describe{
#'   \item{sample.name}{compound name or feature name}
#'   \item{injection.order}{order}
#'   \item{class}{class}
#'   \item{batch}{batch}
#'   \item{group}{group}
#'
#' }
"sample.info"
#' pathway data for human
#'
#' A dataset containing the pathway information:
#'
#' @format A list contain pathway information:
#' \describe{
#'   \item{hsa.kegg.pathway}{ID of all compounds in hsa.kegg.pathway}
#' }
"hsa.kegg.pathway"
#' pathway data for mouse
#'
#' A dataset containing the pathway information:
#'
#' @format A list contain pathway information:
#' \describe{
#'   \item{mmu.kegg.pathway}{ID of all compounds in mmu.kegg.pathway}
#' }
"mmu.kegg.pathway"

