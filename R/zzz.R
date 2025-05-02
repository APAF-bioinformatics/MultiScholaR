# R/zzz.R

.onAttach <- function(libname, pkgname) {
    # Construct the message dynamically using the pkgname variable
    welcome_message <- paste0("Welcome to ", pkgname, "!\n\n")
    load_deps_message <- paste0(
        "IMPORTANT: Please run ", pkgname, "::loadDependencies() \n"
        , "           to ensure all required packages are installed and loaded.\n"
    )

    packageStartupMessage(
        "-------------------------------------------------------------------\n"
        , welcome_message
        , load_deps_message
        , "-------------------------------------------------------------------"
    )
}
