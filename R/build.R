## Build the per genus-data set.
suppressMessages(library(dplyr))
build.woodiness.genus <- function(extreme=FALSE) {
  dat <- load.woodiness()
  lookup <- load.genus.order.lookup()[c("genus", "family", "order")]
  tpl <- load.theplantlist()

  if (!identical(extreme, FALSE))
    dat$woodiness <-
      summarise.count(parse.count(dat$woodiness.count), extreme)

  res <-
    tpl %>% left_join(dat, "gs") %>%
    group_by(genus, family) %>%
    summarise(W=sum(woodiness == "W", na.rm=TRUE),
              V=sum(woodiness == "V", na.rm=TRUE),
              H=sum(woodiness == "H", na.rm=TRUE),
              N=length(woodiness)) %>%
    mutate(K = W + H) %>%
    mutate(p = W / K) %>%
    left_join(lookup, c("genus", "family"))

  cols <- c("genus", "family", "order", "W", "V", "H", "N", "K", "p")
  
  message(sprintf("Final set: %d genera, %d with data, %d species known",
                  nrow(res), sum(res$K > 0), sum(res$K)))

  res[cols]
}

## Check the classification by pulling apart the count.
parse.count <- function(x) {
  res <- t(sapply(strsplit(x, ";", fixed=TRUE), as.integer))
  colnames(res) <- c("H", "V", "W")
  drop(res)
}
summarise.count <- function(x, extreme=FALSE) {
  if (identical(extreme, FALSE)) {
    ans <- ifelse(x[,"W"] > x[,"H"], "W", "H")
    ans[(x[,"W"] == 0 & x[,"H"] == 0 & x[,"V"] > 0) |
        x[,"W"] == x[,"H"]] <- "V"
  } else if (identical(extreme, "woody")) {
    ans <- ifelse(x[,"W"] > 0 | x[,"V"] > 0, "W", "H")
  } else if (identical(extreme, "herbaceous")) {
    ans <- ifelse(x[,"H"] > 0 | x[,"V"] > 0, "H", "W")
  } else {
    stop("Invalid argument for 'extreme'")
  }
  ans
}
