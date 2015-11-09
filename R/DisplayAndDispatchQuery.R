#' Display and dispatch a query to BigQuery.
#'
#' This convenience method will read SQL from a local file or a url.  It will
#' then perform the text substitutions, if any, and dispatch the query to the cloud.
#'
#' @param queryUri File path or url to file containing SQL code.
#' @param project The ID of a Google Cloud Platform project of which the user is a member.
#' @param replacements A list of key/value pairs.  For each pair the key, if found in the
#'  SQL, will be replaced with the value.
#' @return The dataframe of query results.
#' @export
DisplayAndDispatchQuery <- function(queryUri, project, replacements=list()) {
  if (missing(queryUri)) {
    stop("Pass the file path or url to the file containing the query.")
  }
  if(missing(project)) {
    stop("Pass the project id of your Google Cloud Platform project.")
  }

  if (grepl("^https.*", queryUri)) {
    # Read the query from a remote location.
    querySql <- RCurl::getURL(queryUri, ssl.verifypeer=FALSE)
  } else {
    # Read the query from the local filesystem.
    querySql <- readChar(queryUri, nchars=1e6)
  }

  # If applicable, substitute values in the query template.
  for(replacement in names(replacements)) {
    querySql <- gsub(replacement, replacements[[replacement]], querySql, fixed=TRUE)
  }

  # Display the query to the terminal.
  cat(querySql)



  # Dispatch the query to BigQuery.
  bigrquery::query_exec(querySql, project)
}
