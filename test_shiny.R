tryCatch({
    app <- MultiScholaR::run_app()
    print("MultiScholaR app object created successfully.")
}, error = function(e) {
    print(paste("Failed to create app:", e$message))
    quit(status=1)
})
