library(git2r)
library(credentials)

setup_git_auth <- function() {
  # Check if PAT exists
  if (Sys.getenv("GITHUB_PAT") == "") {
    message("No GITHUB_PAT found. Let's set it up...")
    set_github_pat()
  }
  
  # Test authentication
  tryCatch({
    cred <- cred_token()
    # Try to list repositories (test authentication)
    repos <- git2r::repository(".")
    message("Git authentication successful!")
    return(TRUE)
  }, error = function(e) {
    message("Authentication failed: ", e$message)
    return(FALSE)
  })
}

# Run setup
setup_git_auth()