# ----------------------------------------------------------------------------
# testRequiredFiles
# ----------------------------------------------------------------------------
#' @export
testRequiredFiles <- function(files) {
    missing_files <- !file.exists(files)
    invisible(sapply(files[missing_files], function(file) {
        logerror("Missing required file: %s", file)
        q()
    }))
}

# ----------------------------------------------------------------------------
# testRequiredFilesWarning
# ----------------------------------------------------------------------------
#' @export
testRequiredFilesWarning <- function(files) {
    missing_files <- !file.exists(files)
    invisible(sapply(files[missing_files], function(file) {
        logwarn("Missing required file: %s", file)
    }))
}

# ----------------------------------------------------------------------------
# testRequiredArguments
# ----------------------------------------------------------------------------
#' @export
testRequiredArguments <- function(arg_list, parameters) {
    invisible(sapply(parameters, function(par) {
        if (!par %in% names(arg_list)) {
            logerror("Missing required argument: %s", par)
            q()
        }
    }))
}

# ----------------------------------------------------------------------------
# parseType
# ----------------------------------------------------------------------------
#' @export
parseType <- function(arg_list, parameters, functType) {
    invisible(sapply(parameters, function(key) {
        arg_list[key] <- functType(arg_list[key])
    }))
    return(arg_list)
}

# ----------------------------------------------------------------------------
# parseString
# ----------------------------------------------------------------------------
#' @export
parseString <- function(arg_list, parameters) {
    invisible(sapply(parameters, function(key) {
        if (key %in% names(arg_list)) {
            arg_list[key] <- str_replace_all(arg_list[key], c("\"" = "", "\'" = ""))
        }
    }))
    return(arg_list)
}

# ----------------------------------------------------------------------------
# parseList
# ----------------------------------------------------------------------------
#' @export
parseList <- function(arg_list, parameters) {
    invisible(sapply(parameters, function(key) {
        items <- str_replace_all(as.character(arg_list[key]), " ", "")
        arg_list[key] <- base::strsplit(items, split = ",")
    }))
    return(arg_list)
}

# ----------------------------------------------------------------------------
# isArgumentDefined
# ----------------------------------------------------------------------------
#' @export
isArgumentDefined <- function(arg_list, parameter) {
    return(!is.null(arg_list[parameter]) & (parameter %in% names(arg_list)) & as.character(arg_list[parameter]) != "")
}

