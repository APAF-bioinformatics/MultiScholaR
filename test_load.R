options(error = function() { traceback(2); stop() })
tryCatch({
    devtools::load_all(".")
    source("test_shiny.R")
}, error = function(e) {
    print(e)
    quit(status=1)
})
