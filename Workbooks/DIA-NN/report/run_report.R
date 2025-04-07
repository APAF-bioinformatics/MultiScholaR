# Run report wrapper script
# This script ensures the correct Main function is called

# Source the report (this loads all functions without executing them)
rmarkdown::render("Workbooks/DIA-NN/report/DIANN_report.rmd", 
                  params = list(suffix = "Garima_firstpass"),
                  output_file = "DIANN_report_Garima_firstpass.docx")

# Or run the Main function directly with a specific suffix
# source("Workbooks/DIA-NN/report/DIANN_report.rmd", echo = FALSE)
# results <- Main("Garima_firstpass") 

# Note: If you previously had 'main' function defined in your R session,
# you may need to clear your environment first with:
# rm(list = ls()) 