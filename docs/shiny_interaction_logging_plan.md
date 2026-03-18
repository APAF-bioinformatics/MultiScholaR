# Shiny UI Interaction Logging Plan

## 1. Current State Assessment
After reviewing `R/utils_shiny_logging.R` and `R/app_server.R`, it is clear that MultiScholaR currently has a custom logging setup (`appender_shiny` via the `logger` package). However, this setup **only logs backend events and system messages** triggered explicitly by `logger::log_info()`, etc. 

It **does not** currently track user UI interactions (button clicks, dropdown selections, text inputs, etc.). We need to add telemetry for the frontend to enable AI-powered test generation.

## 2. Proposed Implementation Strategy

To capture *every* click and interaction automatically without manually binding observers to hundreds of UI elements, we will integrate the `r_shinylogs` (or `shinylogs` package from CRAN) framework.

### Step 2.1: Dependency Addition
Add the `shinylogs` package to the `Imports:` section of the `DESCRIPTION` file.

### Step 2.2: Server Integration (`R/app_server.R`)
At the very top of `app_server.R` (inside the `app_server` function), we will inject the tracking function. 

```R
app_server <- function(input, output, session) {
  
  # Ensure log directory exists
  log_dir <- file.path(getwd(), "logs", "ui_interactions")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  
  # Track all UI interactions to JSON files
  shinylogs::track_usage(
    storage_mode = shinylogs::store_json(path = log_dir),
    session = session
  )
  
  # (Optional) Track to stdout simultaneously via a custom storage mode
  shinylogs::track_usage(
    storage_mode = shinylogs::store_custom(function(data) {
      if (length(data$inputs) > 0) {
        cat(sprintf("[UI INTERACTION] %s modified at %s\\n", 
            names(data$inputs), Sys.time()))
      }
    }),
    session = session
  )
  
  # ... Existing MultiScholaR setup code ...
}
```

## 3. Log File Location
Based on this plan, the interaction logs will be physically located at:
**`[Project_Root]/logs/ui_interactions/`**

Each user session will automatically generate a highly structured `.json` file in this folder (e.g., `shinylogs_20260319_103025.json`) containing exact timestamps, input IDs, and the resulting values changed. 

## 4. How this Benefits `shinytest2` / `chromote`
Because `shinylogs` captures the exact `input_id` and the `value` payload on every specific click, Gemini AI can simply parse the generated JSON file from `logs/ui_interactions/` and translate it 1:1 into `app$set_inputs()` and `app$click()` commands for your `shinytest2` test suite.
