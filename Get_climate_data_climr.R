# Get_climate_data_climr.R
# Goal: attach monthly climate to herbarium records using the climr package.
# - Gets 1961–1990 climate normals
# - Gets annual climate for each record’s collection year
# - Returns monthly Tave_01–12 and PPT_01–12 plus their normal means
# Sean Michaletz, August 2025 (sean.michaletz@ubc.ca)

# 1. Load packages
library(climr)
library(dplyr)

# Set working directory (edit this path for your machine)
setwd("[path]")

# 2. Read and validate
library(readr)   # for parse_double
library(dplyr)

# Load the input CSV (edit file_in to match your filename)
file_in <- "[filename].csv"  # Update with your file name here
df <- read.csv(file_in, stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN", "N/A"))

# Ensure required columns exist before proceeding
req <- c("Latitude","Longitude","Elevation_m","Year")
if (! all(req %in% names(df))) {
  stop("CSV must have columns: ", paste(req, collapse=", "))
}

# Define function to convert elevation text (character) to meters (numeric)
# - Accepts values like "1234", "1234 m", "4000 ft", or ranges ("1200–1300 m")
# - Converts feet to meters and returns the midpoint for ranges
elev_to_m <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  x0 <- gsub(",", "", x0)                 # remove thousands separators
  x0 <- gsub("\u00A0", " ", x0)           # non-breaking space -> space
  x0 <- gsub("[–—]", "-", x0)             # dashes -> hyphen
  is_ft <- grepl("\\b(ft|feet)\\b", x0)   # any flavor of feet
  # extract all numbers (handles single values and ranges)
  nums  <- regmatches(x0, gregexpr("\\d+\\.?\\d*", x0))
  vals  <- vapply(nums, function(v) {
    if (length(v) == 0) return(NA_real_)
    mean(as.numeric(v))                   # midpoint if a range
  }, numeric(1))
  # convert feet to meters; assume meters otherwise (including no unit)
  vals[is_ft & !is.na(vals)] <- vals[is_ft & !is.na(vals)] * 0.3048
  vals
}

# Keep only valid decimal lat/lon in real-world ranges
lat_num <- suppressWarnings(parse_double(as.character(df$Latitude)))
lon_num <- suppressWarnings(parse_double(as.character(df$Longitude)))
keep_ll <- !is.na(lat_num) & !is.na(lon_num) &
  lat_num >= -90 & lat_num <= 90 &
  lon_num >= -180 & lon_num <= 180
if (sum(!keep_ll)) message("Dropping ", sum(!keep_ll), " rows with non-decimal or out-of-range lat/lon.")

# Parse elevation to meters and year to numeric
elev_m  <- elev_to_m(df$Elevation_m)
year_num <- suppressWarnings(parse_double(as.character(df$Year)))

# Keep rows with good lat/lon and year; overwrite columns with cleaned numerics
keep_all <- keep_ll & !is.na(year_num)

df <- df[keep_all, ] %>%
  mutate(
    Latitude    = lat_num[keep_all],
    Longitude   = lon_num[keep_all],
    Elevation_m = elev_m[keep_all],   # now numeric, meters
    Year        = year_num[keep_all]
  )

# 3. Build coordinates table (adds id; forces western longitudes negative)
df2 <- df %>%
  filter(!is.na(Latitude),
         !is.na(Longitude),
         !is.na(Elevation_m),
         !is.na(Year)) %>%
  mutate(
    lon  = ifelse(Longitude > 0, -Longitude, Longitude),  # assume Western Hemisphere
    lat  = Latitude,
    elev = Elevation_m,
    year = Year,
    id   = row_number()
  )

coords <- df2 %>% select(id, lon, lat, elev)


# 4. List the 24 monthly variables we want (12 temps + 12 precipitation)
vars24 <- c(
  paste0("Tave_", sprintf("%02d", 1:12)),
  paste0("PPT_",  sprintf("%02d", 1:12))
)

# 5. Download 1961–1990 normals for each record (timed)
message("Starting normals downscale at ", Sys.time())
t0_normals <- Sys.time()

normals_raw <- downscale(
  xyz              = coords,
  return_refperiod = TRUE,
  vars             = vars24,
  out_spatial      = FALSE
)

t1_normals <- Sys.time()
message("Finished normals at ", Sys.time(),
        " (", round(difftime(t1_normals, t0_normals, units = "mins"), 2),
        " minutes)")

# Keep id + monthly normals; rename with "Mean_" prefix
normals_df <- as.data.frame(normals_raw) %>%
  select(id, all_of(vars24)) %>%
  rename_with(~ paste0("Mean_", .), .cols = all_of(vars24))

# 6. Download annual climate for the collection year of each record (timed)
message("Starting annual downscale at ", Sys.time())
t0_annual <- Sys.time()

annual_raw <- downscale(
  xyz              = coords,
  obs_years        = df2$year,
  obs_ts_dataset   = "climatena",
  return_refperiod = FALSE,
  vars             = vars24,
  out_spatial      = FALSE
)

t1_annual <- Sys.time()
message("Finished annual at ", Sys.time(),
        " (", round(difftime(t1_annual, t0_annual, units = "mins"), 2),
        " minutes)")

# 7. Select the row that matches each record’s year and keep monthly columns
annual_df2 <- as.data.frame(annual_raw) %>%
  filter(DATASET == "climatena") %>%
  left_join(df2 %>% select(id, year = Year), by = "id") %>%
  filter(PERIOD == as.character(year)) %>%
  select(id, all_of(vars24))

# 8. (Re)build normals table with "Mean_" prefixes (same content; kept for clarity)
normals_df2 <- as.data.frame(normals_raw) %>%
  select(id, all_of(vars24)) %>%
  rename_with(~ paste0("Mean_", .), .cols = all_of(vars24))

# 9. Combine annual climate with normals (one row per specimen)
climate_df <- annual_df2 %>%
  left_join(normals_df2, by = "id")

#   Resulting columns:
#   id, Tave_01..Tave_12, PPT_01..PPT_12,
#   Mean_Tave_01..Mean_Tave_12, Mean_PPT_01..Mean_PPT_12

# 10. Join climate back to records and compute simple “spring” metrics
#     (Uses a 3-month window based on the dataset’s average DOY)
df_final <- df2 %>%
  left_join(climate_df, by = "id") %>%
  { 
    m <- mean(.$DOY, na.rm = TRUE)
    start <- floor((m - 90)/30) + 1
    months <- pmax(1, pmin(start:(start+2), 12))
    t_ann <- paste0("Tave_", sprintf("%02d", months))
    t_nrm <- paste0("Mean_Tave_", sprintf("%02d", months))
    p_ann <- paste0("PPT_", sprintf("%02d", months))
    p_nrm <- paste0("Mean_PPT_", sprintf("%02d", months))
    mutate(.,
           Annual_Spring_T   = rowMeans(select(., all_of(t_ann)),  na.rm=TRUE),
           Normal_Spring_T   = rowMeans(select(., all_of(t_nrm)), na.rm=TRUE),
           Spring_T_Anomaly  = Annual_Spring_T - Normal_Spring_T,
           Annual_Spring_PPT = rowSums(select(., all_of(p_ann)),  na.rm=TRUE),
           Normal_Spring_PPT = rowSums(select(., all_of(p_nrm)), na.rm=TRUE)
    )
  }

# 11. Quick quality checks
#  a) one output row per input record
stopifnot(nrow(df_final) == nrow(df2))
#  b) no duplicate ids
stopifnot(!any(duplicated(df_final$id)))
#  c) climate columns have no missing values
all_monthly <- grep("^(Tave|PPT)_[0-9]{2}$", names(df_final), value=TRUE)
stopifnot(!any(is.na(dplyr::select(df_final, dplyr::all_of(all_monthly))))) # An error here indicates that for some specimen(s), climr didn't return a valid monthly value. If this is the case, check your lat/long values as there is likely an error.
#  d) normals must be the 1961–1990 reference period
if (any(normals_raw$PERIOD != "1961_1990")) {
  stop("Some normals are not 1961-1990!")
}

# 12. Write the final CSV next to the input (adds “_with_climate” to the name)
# e.g., "Data_1901-2023_stm_with_climate.csv"
out_file <- paste0(tools::file_path_sans_ext(file_in), "_with_climate.csv")
write.csv(df_final, out_file, row.names = FALSE)
message("Done: ", out_file, " saved.")
