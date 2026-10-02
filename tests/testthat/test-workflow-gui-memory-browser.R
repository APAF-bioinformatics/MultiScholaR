# DIA-025 standalone browser attribution. No full GUI journey or application change.
test_that("bounded Glimma replay isolates table DOM allocation without a GUI journey", {
  receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
  skip_if(!nzchar(receipt), "guarded diagnostic only")
  library <- Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY")
  expect_identical(normalizePath(getNamespaceInfo(asNamespace("MultiScholaR"), "path")),
    normalizePath(file.path(library, "MultiScholaR")))
  n <- 1024L
  annotation <- as.data.frame(stats::setNames(rep(list(paste0("value_", seq_len(n))), 48L),
    paste0("annotation_", seq_len(48L))))
  # Glimma supplies its reserved `gene` key from the counts row names.
  annotation$Protein_Ids <- paste0("P", seq_len(n))
  counts <- matrix(seq_len(n * 12L), n, 12L,
    dimnames = list(annotation$Protein_Ids, paste0("sample_", seq_len(12L))))
  widget <- Glimma::glimmaXY(x = seq_len(n) / n, y = seq_len(n) / n,
    anno = annotation, counts = counts,
    groups = factor(rep(c("A", "B"), each = 6L)), transform.counts = "none")
  expect_identical(anyDuplicated(widget$x$data$cols), 0L)
  output <- file.path(dirname(receipt), "glimma-replay")
  stopifnot(dir.create(output))
  observation <- list(status = "running", application_helper = "deferProtGlimmaTableRows",
    meaning = "isolated Glimma DOM ownership, not exact large-run RSS attribution")
  on.exit(jsonlite::write_json(observation,
    file.path(dirname(receipt), "repair-browser-table.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null"), add = TRUE)
  browser <- chromote::ChromoteSession$new()
  on.exit(browser$close(), add = TRUE)
  browser$Page$addScriptToEvaluateOnNewDocument(source = "
    window.diaErrors=[];
    addEventListener('error', function(e) {
      if(window.diaErrors.length<16) window.diaErrors.push(String(e.message ||
        (e.target && e.target.src) || 'resource error').slice(0,1000));
    }, true);
    addEventListener('unhandledrejection', function(e) {
      if(window.diaErrors.length<16) window.diaErrors.push(String(e.reason).slice(0,1000));
    });
    window.alert=function(message) {
      if(window.diaErrors.length<16) window.diaErrors.push('alert: '+String(message).slice(0,1000));
    };")
  observe <- function(deferred) {
    page <- file.path(output, paste0("glimma-", deferred, ".html"))
    trial <- widget
    if (deferred) trial <- deferProtGlimmaTableRows(trial)
    expect_identical(trial$x, widget$x)
    htmlwidgets::saveWidget(trial, page, selfcontained = FALSE)
    navigation <- browser$Page$navigate(paste0("file://", page))
    if (!is.null(navigation$errorText)) stop(navigation$errorText)
    query <- "JSON.stringify((function(){var result={ready:document.readyState,
      errors:window.diaErrors || [],jquery:!!window.jQuery,tableCount:0,
      datatables:!!(window.jQuery && jQuery.fn.dataTable),
      bindings:window.HTMLWidgets ? HTMLWidgets.widgets.map(function(b){return b.name;}) : [],
      widgetCount:document.querySelectorAll('.html-widget').length,
      boundCount:document.querySelectorAll('.html-widget-static-bound').length,
      plotContainers:document.querySelectorAll('.plotContainer').length,
      tablesInDom:document.querySelectorAll('table').length};
      if(!window.jQuery || !jQuery.fn.dataTable) return result;
      var tables=jQuery.fn.dataTable.tables();result.tableCount=tables.length;
      if(!tables.length) return result;
      var table=jQuery(tables[0]).DataTable();result.records=table.rows().count();
      result.rowNodes=table.rows().nodes().toArray().filter(Boolean).length;result.columns=table.columns().count();
      result.deferred=table.settings()[0].oFeatures.bDeferRender;
      result.globalDeferred=jQuery.fn.dataTable.defaults.bDeferRender;return result;})())"
    deadline <- Sys.time() + 30
    requested_render <- FALSE
    inspected_render <- FALSE
    repeat {
      evaluated <- browser$Runtime$evaluate(query, returnByValue = TRUE)
      value <- evaluated$result$value
      data <- if (is.character(value)) jsonlite::fromJSON(value) else NULL
      observation$latest <<- list(deferred = deferred, data = data,
        exception = evaluated$exceptionDetails, requested_render = requested_render)
      if (!is.null(data) && isTRUE(data$records == n)) break
      if (!inspected_render && identical(data$ready, "complete") &&
          isTRUE(data$boundCount == 1L) && isTRUE(data$tableCount == 0L) &&
          Sys.time() > deadline) {
        inspected_render <- TRUE
        direct <- browser$Runtime$evaluate("(function(){
          var el=document.querySelector('.glimmaXY');
          var instance=HTMLWidgets.getInstance(el);
          var script=document.querySelector('script[type=\"application/json\"][data-for=\"'+el.id+'\"]');
          var result={instanceKeys:Object.keys(instance || {}),hasScript:!!script};
          try { var payload=JSON.parse(script.textContent);result.payloadKeys=Object.keys(payload);
            instance.renderValue(payload.x);result.invoked=true;
          } catch(e) {result.error=String(e);result.stack=e.stack;}
          return JSON.stringify(result);})()", returnByValue = TRUE)
        observation$direct_render <<- list(
          result = if (is.character(direct$result$value)) jsonlite::fromJSON(direct$result$value) else NULL,
          exception = direct$exceptionDetails)
      }
      if (!requested_render && identical(data$ready, "complete") &&
          isTRUE(data$boundCount == 0L) && "glimmaXY" %in% data$bindings) {
        requested_render <- TRUE
        rendered <- browser$Runtime$evaluate("HTMLWidgets.staticRender()")
        observation$render_exception <<- rendered$exceptionDetails
      }
      if (Sys.time() > deadline) stop("bounded Glimma table did not initialize; see retained page/errors")
      Sys.sleep(0.1)
    }
    data$dom <- browser$Memory$getDOMCounters()
    result <- browser$Runtime$evaluate("(function(){var t=jQuery(jQuery.fn.dataTable.tables()[0]).DataTable();
      t.search('P1024').draw();return t.rows({search:'applied'}).count();})()",
      returnByValue = TRUE)$result$value
    expect_equal(result, 1)
    selected <- browser$Runtime$evaluate("(function(){var t=jQuery(jQuery.fn.dataTable.tables()[0]).DataTable();
      jQuery(t.rows({search:'applied'}).nodes()).trigger('click');
      return document.querySelector('.geneDisplay').textContent;})()",
      returnByValue = TRUE)$result$value
    expect_match(selected, "P1024", fixed = TRUE)
    cleared <- browser$Runtime$evaluate("(function(){jQuery('.clearSubset').trigger('click');
      return document.querySelector('.geneDisplay').textContent;})()",
      returnByValue = TRUE)$result$value
    expect_identical(cleared, "")
    data
  }
  tryCatch({
    observation$native <- observe(FALSE)
    observation$deferred_trial <- observe(TRUE)
    expect_identical(observation$native$tableCount, 1L)
    expect_identical(observation$deferred_trial$tableCount, 1L)
    expect_length(observation$native$errors, 0L)
    expect_length(observation$deferred_trial$errors, 0L)
    expect_null(observation$direct_render)
    expect_false(observation$native$globalDeferred)
    expect_false(observation$deferred_trial$globalDeferred)
    expect_identical(observation$native$records, observation$deferred_trial$records)
    expect_identical(observation$native$columns, observation$deferred_trial$columns)
    expect_lt(observation$deferred_trial$rowNodes, observation$native$rowNodes)
    observation$status <- "observed"
  }, error = function(error) {
    observation$status <<- "failed"
    observation$error <<- conditionMessage(error)
    stop(error)
  })
})
