% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BNrich.R
\name{unify_path}
\alias{unify_path}
\title{Simplification networks -- applied to unifying nodes}
\usage{
unify_path(dataH, dataD, MapkG, Pathway.id)
}
\arguments{
\item{dataH}{A data frame contains (healthy) control objects data}

\item{dataD}{A data frame contains disease objects data}

\item{MapkG}{A list contains imported 187 signaling pathways}

\item{Pathway.id}{A vector contains 187 KEEG pathway IDs}
}
\value{
A list contain data_h,data_d,mapkG1 and pathway.id1
}
\description{
Unifying nodes based imported signaling pathways and GE data
}
\examples{
#All the 187 preprocessed signaling pathways can be entered in analysis by fetch_data_file().
#But here you enter a subset of those pathways to see how this package works.
files <- system.file("extdata", "test_files_to_start.RData", package = "BNrich", mustWork = TRUE)
load(files)
Data <- system.file("extdata", "Test_DATA.RData", package = "BNrich", mustWork = TRUE)
load(Data)
uni_Result <- unify_path(dataH, dataD, MapkG = sub_mapkG, Pathway.id = path.id)
}
