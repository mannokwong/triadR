trmGlobal <- new.env(parent=globalenv())
trmGlobal$suffixes <- c(".j", ".a", ".p",
                        ".ja", ".jp", ".ap", ".aj", ".pj", ".pa",
                        ".jap", ".jpa", ".apj", ".ajp", ".pja", ".paj")
trmGlobal$style <- "perception"


# set options for printing results etc.
#' @export
trm.style <- function(style="perception", suffixes = NA) {
  trmGlobal$style <- style <- match.arg(style, c("perception", "cube"))

  if (any(is.na(suffixes))) {
    if (style=="perception") {
      trmGlobal$suffixes <- c(".j", ".a", ".p",
                              ".ja", ".jp", ".ap", ".aj", ".pj", ".pa",
                              ".jap", ".jpa", ".apj", ".ajp", ".pja", ".paj")
    } else
      if (style=="cube") {
        trmGlobal$suffixes <- c(".l", ".r", ".c",
                                ".lr", ".lc", ".rc", ".rl", ".cl", ".cr",
                                ".lrc", ".lcr", ".rcl", ".rlc", ".clr", ".crl")
      }

  } else {
    trmGlobal$suffixes <- suffixes
  }
}

