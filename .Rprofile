options(
  languageserver.server_capabilities = list(
    definitionProvider = TRUE
  ),
  languageserver.snippet_support = FALSE,
  r.lsp.rich_documentation = FALSE,
  languageserver.formatting_style = function(options) {
    styler::tidyverse_style(scope = "indention", indent_by = options$tabSize)
  }
)

setHook(
  packageEvent("languageserver", "onLoad"),
  function(...) {
    options(languageserver.default_linters = lintr::with_defaults(
      line_length_linter = NULL,
      object_name_linter = NULL,
      object_usage_linter = NULL
    ))
  }
)
