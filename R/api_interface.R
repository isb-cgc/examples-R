
# interface to end points via R using httr #

# David L Gibbs
# Institute for Systems Biology
# March 4, 2016

################################################################################
#  google httr tutorial
#  Huge thanks to Hadley: https://github.com/hadley/httr/blob/master/demo/oauth2-google.r
#
# 1. Find OAuth settings for google:
#    https://developers.google.com/accounts/docs/OAuth2InstalledApp
#oauth_endpoints("google")
#
# 2. Register an application at https://cloud.google.com/console#/project
#myapp <- oauth_app("google", key = CLIENT_ID, secret = CLIENT_SECRET)
#
# 3. Get OAuth credentials
#google_token <- oauth2.0_token(oauth_endpoints("google"), myapp, scope = "https://www.googleapis.com/auth/userinfo.profile")
#
# 4. Use API
#req <- GET("https://www.googleapis.com/oauth2/v1/userinfo", config(token = google_token))
#stop_for_status(req)
#content(req)
###############################################################################


library(httr)

#' Initialize the authorization
#'
#' Creates a auth token needed for cohort related functions.
#'
#' @details Uses the httr package to create an OAuth2.0 token with Google.
#'
#' @return An Oauth token
#' @examples
#' \dontrun{
#' isb_init()
#' }
#' @export
isb_init <- function() {
	library(httr)
	# for native application - same as settings.INSTALLED_APP_CLIENT_ID
	CLIENT_ID = "907668440978-0ol0griu70qkeb6k3gnn2vipfa5mgl60.apps.googleusercontent.com"
	# NOTE: this is NOT actually a 'secret' -- we're using the 'installed
	# application' OAuth pattern here
	CLIENT_SECRET = "To_WJH7-1V-TofhNGcEqmEYi"
	EMAIL_SCOPE = "https://www.googleapis.com/auth/userinfo.email"
	myapp <- oauth_app("google", key = CLIENT_ID, secret = CLIENT_SECRET)
	isb_token <- oauth2.0_token(oauth_endpoints("google"), myapp, scope = EMAIL_SCOPE)
	isb_token
}


#' Data file filter
#'
#' Filter the results from datafilenamekeys
#'
#' The call to the datafilenamekeys API returns a list with platform types
#' in the strings. This filters down the list of files by platform.
#'
#' @param filelist The return data, a vector of characters, from sample_files
#' @param filterby A string, naming a platform
#' @return List of files
#'
#' @examples
#' \dontrun{
#'   fs <- sample_files("TCGA-A7-A6VV-10A")
#'   datafile_filter(sf, ".bam")
#'  }
#' @export
datafile_filter <- function(filelist, filterby) {
	x <- filelist$datafilenamekeys
  y <- x[str_detect(x, filterby)]
	unlist(y)
}

#' Get sample files
#'
#' Get file locations using sample barcode
#'
#' Uses the API to get a list of files for a given barcode.
#'
#' @param sample_barcode A TCGA barcode.
#' @return List of files across platforms.
#' @examples
#' \dontrun{
#'   fs <- sample_files("TCGA-A7-A6VV-10A")
#'  }
#' @export
sample_files <- function(sample_barcode) {
	b <- "https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/datafilenamekey_list_from_sample"
	content(GET(url=b, query=list(sample_barcode=sample_barcode)))
}


#' Get sample details
#'
#' Get details on a particular sample, given a sample barcode
#'
#' @param sample_barcode A TCGA sample barcode
#'
#' @returns List Details on a given sample_details
#'
#' @examples
#' \dontrun{
#'   sample_details("TCGA-A7-A6VV-10A")
#' }
#' @export
sample_details <- function(sample_barcode) {
	b <- "https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/sample_details"
	content(GET(b, query=list(sample_barcode=sample_barcode)))
}


#' Get patient details
#'
#' Get details on a particular patient, given a sample barcode
#'
#' @param sample_barcode A TCGA patient barcode like "TCGA-02-0001"
#'
#' @returns List Details on a given sample_details
#'
#' @examples
#' \dontrun{
#'   patient_details("TCGA-A7-A6VV")
#'  }
#' @export
patient_details <- function(patient_barcode) {
	# TCGA-44-7670
	b = "https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/patient_details"
	content(GET(b, query=list(patient_barcode = patient_barcode)))
}

#' List cohorts
#'
#' Using your auth token, get a list of saved cohorts
#'
#' @param a_token Auth token as returned by isb_init()
#'
#'
#' @examples
#' \dontrun{
#'   list_cohorts(mytoken)
#'  }
#' @export
list_cohorts <- function(a_token) {
	req <- GET("https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/cohorts_list", config(token = a_token))
	stop_for_status(req)
	content(req)
}


#' Save cohorts
#'
#' Using your auth token, save a new cohort
#'
#' @param a_token Your auth token
#' @param cohort_name The name of the new cohort
#' @param filter_list The list of filter parameters defining the cohort
#' @returns List
#'
#' @examples
#' \dontrun{
#'   save_cohorts(mytoken, "new_cohort", list(list(Study="BRCA")))
#' }
#' @export
save_cohorts <- function(a_token, cohort_name, filter_list) {
	b <- paste0("https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/save_cohort?name=",cohort_name)
	req <- POST(b, query=filter_list, config(token = a_token))
	stop_for_status(req)
	content(req)
}


#' Preview cohorts
#'
#' Using your auth token, create a cohort, without saving
#'
#' @param filter_list The list of filter parameters defining the cohort
#'
#' @returns List Details on the cohort
#'
#' @examples
#' \dontrun{
#'   preview_cohort(list(Study="BRCA"))
#' }
#' @export
preview_cohort <- function(filter_list) {
	# example: preview_cohort(list(Study="BRCA"))
	require(httr)
	b <- "https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/preview_cohort"
	content(POST(b, body=filter_list))
}


#' Data files from a cohort
#'
#' Using a defined cohort, get file paths to a Google bucket
#'
#' @param cohort_id The list of filter parameters defining the cohort
#' @param pipeline The pipeline used in data generation
#' @param platform The platform used for data generation
#' @param token Your auth token
#'
#' @returns List with count and datafilenamekeys
#'
#' @examples
#' \dontrun{
#'   datafiles_from_cohort("69", NULL, "IlluminaHiSeq_RNASeqV2", mytoken)
#' }
#' @export
datafiles_from_cohort <- function(cohort_id, pipeline=NULL, platform=NULL, a_token) {
	require(httr)
	b <- "https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/datafilenamekey_list_from_cohort"
	content(GET(b, query=list(cohort_id=cohort_id, pipeline=pipeline, platform=platform), config(token = mytoken)))
}

#' List of barcodes from a cohort
#'
#' Using a defined cohort, get sample and patient barcodes
#'
#' @param cohort_id The list of filter parameters defining the cohort
#' @param token Your auth token
#'
#' @returns List of barcodes
#'
#' @examples
#' \dontrun{
#'   sample_list_cohort(list("29"), mytoken)
#' }
#' @export
barcodes_from_cohort <- function(cohort_id, a_token) {
  require(httr)
  b <- "https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/cohort_patients_samples_list"
  content(GET(b, query=list(cohort_id=cohort_id), config(token = mytoken)))
}



#' Search data details list
#'
#' sample_details brings a list of data_details, this function assists in searching the list
#'
#' @param deets sample_details returned list
#' @param platform data platform used for data generation
#' @param level the TCGA level of the data
#'
#' @returns list of particulars about a given platform / level
#'
#' @examples
#' \dontrun{
#'   search_sample_details(deets)
#' }
#' @export
search_sample_details <- function(deets, platform, level) {
  n <- length(deets)

  }
