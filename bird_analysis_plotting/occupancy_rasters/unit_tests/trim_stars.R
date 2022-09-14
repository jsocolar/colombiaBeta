library(stars)

s1 <- st_as_stars(matrix(c(NA, 2, 3, NA, 5:6, rep(NA, 3)), nrow=3))
s2 <- st_as_stars(matrix(c(1:6, rep(NA, 3)), nrow=3))

s3 <- do.call(c, list(s1, s1, s1)) %>%
    merge

s4 <- do.call(c, list(s2, s2, s2)) %>%
    merge

s5 <- c(s3, s4)

y <- trim_stars(s5)

dim(y) == c(3, 2, 3)