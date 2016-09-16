.onLoad <- function(libname = find.package("bigmelon"), pkgname = "bigmelon") {
  options(runLast = TRUE)

  .GlobalEnv$.Last <- function(x = ls(.GlobalEnv)){
    for(i in x){
      if(inherits(get(i), "gds.class")){
#        cat("Closing any unclosed gds files! \n")
        a <- try(closefn.gds(get(i)), silent = TRUE)
#        if(!inherits(a, "try-error")) cat('Closed:', get(i)$filename, '\n')
      }
    }
  }
}


