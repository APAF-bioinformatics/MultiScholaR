# Author(s): Ignatius Pang
# Email: ipang@cmri.org.au
# Children’s Medical Research Institute, finding cures for childhood genetic diseases

##################################################################################################################

### Function: create_id_to_attribute_hash
### Description: Create a hash function that map keys to attributes.

## Inputs:
## keys: An array of key values
## attributes: An array of attribute values

## Output:
## An environment that act as a hash to convert keys to attributes.
#' @export
createIdToAttributeHash <- function(keys, attributes) {

	keys <- as.character( as.vector(keys))
	attribute <- as.vector(attributes)

	hash <- new.env(hash = TRUE, parent = parent.frame())

	if ( length(keys) != length(attributes))  {
		warning('Length of keys != Length of attributes list.')
		return(1)
	}

	for ( i in 1:length(keys) )  {
		base::assign( keys[i], attributes[i], envir = hash)
	}

	return(hash)
}

##################################################################################################################

### Function: create_id_to_attribute_hash
### Description: Use a predefined hash dictionary to convert any Key to Attribute, return NA if key does not exists

## Inputs:
## key: A key value
## hash: The hash dictionary that maps keys to attributes

## Output:
## A value that correspond to the query key value.
#' @export
convertKeyToAttribute <- function(key, hash) {

	if ( base::exists(key, hash) ) {
		return ( base::get(key, hash))
	}
	else {
		return (NA)
	}
}

##################################################################################################################

#TODO: Move to an Rpackage RCMRI with common functions. It may be a dependency for each CMRI's project.
#' @export
createDirectoryIfNotExists <- function(file_path, mode = "0777") {

	### Create directory recursively if it doesn't exist
	if (! file.exists(file_path)){

		dir.create(file_path, showWarnings = TRUE, recursive = TRUE, mode = mode)

	}

}

#' @export
createDirIfNotExists  <- function(file_path, mode = "0777") {
  createDirectoryIfNotExists(file_path, mode = "0777")
}


##################################################################################################################

## Function to source Rmd files
# https://stackoverflow.com/questions/10966109/how-to-source-r-markdown-file-like-sourcemyfile-r
#' @export
sourceRmdFileSimple <- function(x, ...) {
	source(purl(x, output = tempfile()), ...)
}

##################################################################################################################

#' https://gist.github.com/noamross/a549ee50e8a4fd68b8b1
#' Source the R code from an knitr file, optionally skipping plots
#'
#' @param file the knitr file to source
#' @param skip_plots whether to make plots. If TRUE (default) sets a null graphics device
#'
#' @return This function is called for its side effects
#' @export
sourceRmdFile <- function(file, skip_plots = TRUE) {
	temp = tempfile(fileext=".R")
	knitr::purl(file, output=temp)

	if(skip_plots) {
		old_dev = getOption('device')
		options(device = function(...) {
			.Call("R_GD_nullDevice", PACKAGE = "grDevices")
		})
	}
	source(temp)
	if(skip_plots) {
		options(device = old_dev)
	}
}


##################################################################################################################


#TODO: create an Rpackage RCMRI with common functions. It may be a dependency for each CMRI's project.
#=====================================================================================================
#' @export
createOutputDir <- function(output_dir, no_backup) {
  if (output_dir == "") {
    logerror("output_dir is an empty string")
    q()
  }
  if (dir.exists(output_dir)) {
    if (no_backup) {
      unlink(output_dir, recursive = TRUE)
    }
    else {
      backup_name <- paste(output_dir, "_prev", sep = "")
      if (dir.exists(backup_name)) { unlink(backup_name, recursive = TRUE) }
      system(paste("mv", output_dir, backup_name)) }
  }
  dir.create(output_dir, recursive = TRUE)
}

#' @export
cmriWelcome <- function(name, autors) {
  loginfo("   ______     __    __     ______     __       ")
  loginfo('  /\\  ___\\   /\\ "-./  \\   /\\  == \\   /\\ \\      ')
  loginfo("  \\ \\ \\____  \\ \\ \\-./\\ \\  \\ \\  __<   \\ \\ \\     ")
  loginfo("   \\ \\_____\\  \\ \\_\\ \\ \\_\\  \\ \\_\\ \\_\\  \\ \\_\\    ")
  loginfo("    \\/_____/   \\/_/  \\/_/   \\/_/ /_/   \\/_/    ")
  loginfo("")
  loginfo("---- Children’s Medical Research Institute ----")
  loginfo(" Finding cures for childhood genetic diseases  ")
  loginfo("")
  loginfo(" ==============================================")
  loginfo(" %s", name)
  loginfo(" Author(s): %s", paste(autors, sep = ", "))
  loginfo(" cmri-bioinformatics@cmri.org.au")
  loginfo(" ==============================================")
  loginfo("")
}

#' @export
testRequiredFiles <- function(files) {
  missing_files <- !file.exists(files)
  for (file in files[missing_files]) {
    logerror("Missing required file: %s", file)
    q()
  }
}

#' @export
testRequiredArguments <- function(arg_list, parameters) {
  for (par in parameters) {
    if (!par %in% names(arg_list)) {
      logerror("Missing required argument: %s", par)
      q()
    }
  }
}

#' @export
parseType<-function (arg_list,parameters,functType){
  for(key in parameters){
    arg_list[key]<-functType(arg_list[key])
  }
  return (arg_list)
}

#' @export
parseString<-function (arg_list,parameters){
  for(key in parameters){
    if(key %in% names(arg_list)) {
      arg_list[key] <- str_replace_all(arg_list[key], c("\"" = "", "\'" = ""))
    }
  }
  return (arg_list)
}

#' @export
parseList<-function (arg_list,parameters){
  for(key in parameters){
    items<-str_replace_all(as.character(arg_list[key])," ","")
    arg_list[key]<- base::strsplit(items,split=",")
  }
  return (arg_list)
}

#' @export
isArgumentDefined<-function(arg_list,parameter){
  return (!is.null(arg_list[parameter]) & (parameter %in% names(arg_list)) & as.character(arg_list[parameter]) != "")
}

#' @export
cmriFormatter <- function(record) { sprintf('CMRI Bioinformatics %s [%s] | %s', record$levelname, record$timestamp, record$msg) }


#' @export
setArgsDefault <- function(args, value_name, as_func, default_val=NA ) {

  if(isArgumentDefined(args, value_name))
  {
    args<-parseType(args,
                    c(value_name)
                    ,as_func)
  }else {
    logwarn(paste0( value_name, " is undefined, default value set to ", as.character(default_val), "."))
    args[[ value_name ]] <- default_val
  }

  return(args)
}

#=====================================================================================================
