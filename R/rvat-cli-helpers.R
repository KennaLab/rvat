#' @include rvat_cli_options.R

# collect & parse command-line args
collect_args <- function() {
  parsed_opts <- optparse::OptionParser(option_list = rvat_cli_options, add_help_option = FALSE)
  args <- optparse::parse_args(parsed_opts, print_help_and_exit = FALSE)
  args_raw <- unlist(lapply(strsplit(commandArgs(TRUE), split = "=|\\s"), function(x) x[[1]]))
  args_raw <- args_raw[grepl("^-", args_raw)]
  args_raw <- gsub("^-{1,2}", "", args_raw)
  list(args = args, args_raw = args_raw)
}

# print specified arguments
print_args <- function(args, args_raw, flags) {
  message("The following options are specified:\n\n")
  args_raw <- args_raw[args_raw != "help"]
  for(i in args_raw) {
    if(i %in% flags) {
      message(sprintf("--%s\n", i))
    } else {
      message(sprintf("--%s=%s\n", i, args[[i]]))
    }
  }
  cat("\n")
}

# check expected arguments
check_expected_args <- function(args, expected, func_name) {
  if(!all(args %in% expected)) {
    stop(sprintf("Unexpected argument(s) for %s: %s",
                 func_name, paste(paste0("--", args[!args %in% expected]), collapse=",")),
         call. = FALSE
    )
  }
}

# check required arguments
check_required_args <- function(args, required,required_one_of = NULL, func_name) {
  if(!all(required %in% args)) {
    stop(sprintf("The following required arguments for %s are missing: %s",
                 func_name,
                 paste(paste0("--", required[!required %in% args]), collapse=",")),
         call. = FALSE
    )
  }
  if(!is.null(required_one_of)) {
    if(sum(required_one_of %in% args) == 0) {
      stop(sprintf("One of the following arguments should be specified: %s",
                   paste(paste0("--", required_one_of), collapse=",")),
           call. = FALSE
      )
    }
  }
}

# check arguments
check_args <- function(func_name, args, help, required, expected, required_one_of = NULL) {
  check_required_args(args = args, required = required, required_one_of = required_one_of, func_name = func_name)
  check_expected_args(args = args, expected = expected, func_name = func_name)
}

# stop cli quietly
stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
} 

# print help pages
check_help <- function (args_raw) {
  
  if ("help" %in% args_raw && sum(rvat_cli_methods_flags %in% args_raw) == 1) {
      message(rvat_cli_help[[rvat_cli_methods_flags[rvat_cli_methods_flags %in% args_raw]]])
      stopQuietly()
    } else if ( "help" %in% args_raw && sum(rvat_cli_methods_flags %in% args_raw) != 1 ) {
      message("")
      stopQuietly()
    }
  
  if ( (length(args_raw) == 1) && args_raw %in% rvat_cli_methods_flags ) {
    message(rvat_cli_help[[args_raw]])
    stopQuietly()
  }
  invisible(NULL)
}

.retrieve_formals_dispatch <- function(generic, class) {
  method <- getMethod(generic, class)
  formals(body(method@.Data)[[2]][[3]])
}