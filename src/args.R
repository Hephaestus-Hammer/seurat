# Small argument helpers shared by the entrypoints.

# ob passes every parameter as text: "true"/"false" (also TRUE/FALSE, 1/0, yes/no) -> logical.
as_flag <- function(x, name) {
  v <- tolower(trimws(as.character(x)))
  if (v %in% c("true", "t", "1", "yes")) return(TRUE)
  if (v %in% c("false", "f", "0", "no")) return(FALSE)
  stop(name, " must be true or false, not ", x)
}
