% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DisplayAndDispatchQuery.R
\name{DisplayAndDispatchQuery}
\alias{DisplayAndDispatchQuery}
\title{Display and dispatch a query to BigQuery.}
\usage{
DisplayAndDispatchQuery(queryUri, project, replacements = list())
}
\arguments{
\item{queryUri}{File path or url to file containing SQL code.}

\item{project}{The ID of a Google Cloud Platform project of which the user is a member.}

\item{replacements}{A list of key/value pairs.  For each pair the key, if found in the
SQL, will be replaced with the value.}
}
\value{
The dataframe of query results.
}
\description{
This convenience method will read SQL from a local file or a url.  It will
then perform the text substitutions, if any, and dispatch the query to the cloud.
}

