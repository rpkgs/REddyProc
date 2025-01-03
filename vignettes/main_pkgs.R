pacman::p_load(
  data.table, dplyr, lubridate
)

make_date_yj <- function(year, doy, tz = "UTC") {
  make_date(year, tz = tz) + ddays(doy - 1) 
}

make_datetime_yjh <- function(year, doy, hour = 0, tz = "UTC") {
  make_datetime(year, tz = tz) + ddays(doy - 1) + dhours(hour)
}

#' create_season
#' @param time a vector of POSIXct
#' @param info season info, with the column of `year` and `doy`
#' @export 
create_season <- function(time, info) {
  year_beg <- info$year %>% min()
  year_end <- info$year %>% max()
  date_beg <- make_datetime(year_beg)
  date_end <- make_datetime(year_end, 12, 31, 23, 59, 59)

  if (info$doy[1] != 1) {
    info %<>% rbind(data.frame(doy = 1, year = year_beg), .)
  }
  info %<>% mutate(flag = sprintf("%d%03d", year, doy))

  brks <- make_date_yday(info$year, info$doy) %>%
    as.POSIXct() %>%
    c(date_end)
  inds <- cut(time, breaks = brks, labels = FALSE)
  as.factor(info$flag)[inds]
}
