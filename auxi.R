
k2 <- function (df){
    kable( df , format = "html", booktabs = T, caption = "cancer studies") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"))
}