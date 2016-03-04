
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

# creates a token needed for cohort related functions.
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


# utils
datafile_filter <- function(filelist, filterby) {
	x <- filelist[[1]]$datafilenamekeys
  x[str_detect(x, filterby)]
}

# get file locations using sample barcode
sample_files <- function(sample_barcode) {
	b <- "https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/datafilenamekey_list_from_sample"
	content(GET(url=b, query=list(sample_barcode=sample_barcode)))
}


# get details on a particular sample
sample_details <- function(sample_barcode) {
	b <- "https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/sample_details"
	content(GET(b, query=list(sample_barcode=sample_barcode)))
}


# get details on a patient given a patient barcode
patient_details <- function(patient_barcode) {
	# TCGA-44-7670
	b = "https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/patient_details"
	content(GET(b, query=list(patient_barcode = patient_barcode)))
}


# list cohorts, get your token from the isb_init function.
list_cohorts <- function(a_token) {
	req <- GET("https://api-dot-mvm-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/cohorts_list", config(token = a_token))
	stop_for_status(req)
	content(req)
}


# list cohorts, get your token from the isb_init function.
save_cohorts <- function(a_token, cohort_name, filter_list) {
	b <- paste0("https://api-dot-mvm-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/save_cohort?name=",cohort_name)
	req <- POST(b, query=filter_list, config(token = a_token))
	stop_for_status(req)
	content(req)
}


# preview cohorts for the current user
preview_cohort <- function(filter_list) {
	# example: preview_cohort(list(Study="BRCA"))
	require(httr)
	b <- "https://api-dot-isb-cgc.appspot.com/_ah/api/cohort_api/v1/preview_cohort"
	content(POST(b, body=filter_list))
}
