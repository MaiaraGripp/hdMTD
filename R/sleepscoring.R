#' Data with sleeping patterns
#'
#' The dataset contains 136,925 rows and 7 columns, representing the sleeping patterns of 151 patients
#' over the course of one night, with measurements taken at 30-second intervals. It is a collection of 151
#'  whole-night polysomno-graphic (PSG) sleep recordings (85 Male, 66 Female, mean age of 53.9 ± 15.4)
#'  collected during 2018 at the Haaglanden Medisch Centrum (HMC, The Netherlands) sleep center.
#'
#'  This dataset was acquired at Physionet.org. URL https://doi.org/10.13026/t79q-fr32.
#'
#' @format A tbl_df object with 136,925 rows representing the sleeping patterns of 151 patients.
#'
#' \describe{
#'   \item{Patient}{Identifies the patient.}
#'   \item{Date}{Date when the measurements where made.}
#'   \item{Time}{Time each measurement was made.}
#'   \item{Recorging.onset}{Time, in seconds, since the beginning of the recordings.}
#'   \item{Duration}{The duration of each lag between recordings, in seconds.}
#'   \item{SleepStage}{The annotated sleeping stage. 'W' refers to wakefulness, 'R' to REM sleep, and 'N1', 'N2', and 'N3' refer to non-REM stages 1, 2, and 3 respectively.}
#'   \item{ConscientLevel}{A measurement of the level of consciousness during different sleep stages. 0 indicates Wake, 1 represents REM sleep, and 2, 3, and 4 correspond to N1, N2, and N3 stages respectively.}
#' }
#'
#' @source {Alvarez-Estevez D, Rijsman RM (2022). Haaglanden Medisch Centrum sleep staging database (version 1.1). PhysioNet. URL https://doi.org/10.13026/t79q-fr32.}
#' @source {Alvarez-Estevez D, Rijsman RM (2021). Inter-database validation of a deep learning approach for automatic sleep scoring. PLOS ONE, 16(8), 1–27. URL https://doi.org/10.1371/journal.pone.0256111.}
#' @source {Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PC, Mark RG, Mietus JE, Moody GB, Peng CK, Stanley HE (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a New Research Resource for Complex Physiologic Signals. Circulation, 101(23), e215–e220.}
#'
#' @examples
#' data(sleepscoring)
"sleepscoring"
