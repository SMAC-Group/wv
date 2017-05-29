#' @importFrom utils packageDescription
.onAttach = function(...) {
  
  local_version = utils::packageDescription('wv')
  
  packageStartupMessage('Version: ', local_version$Version, ' built on ', local_version$Date)
  packageStartupMessage('To see the user guides use: browseVignettes("wv").')
}
