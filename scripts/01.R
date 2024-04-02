
library(tidyverse)
library(stars)
library(units)
library(furrr)

options(future.fork.enable = T)
plan(multicore)
options(future.globals.maxSize = 1000*1024^2)

source("https://raw.github.com/carlosdobler/spatial-routines/master/tile.R")
source("https://raw.github.com/carlosdobler/spatial-routines/master/tile_land.R")


# cloud directory of raw input data
dir_gs <- "/mnt/bucket_mine/cmip6/nex/daily"

# local directory where temp files will be stored
dir_disk <- "/mnt/pers_disk_"

# cloud directory where final results will be uploaded
dir_res <- "gs://clim_data_reg_useast1/cmip6/nex"


models <- c("GFDL-ESM4",
            "MPI-ESM1-2-HR",
            "MRI-ESM2-0",
            "UKESM1-0-LL",
            "IPSL-CM6A-LR")


var_combos <- list(c("maximum_temperature", "minimum_temperature") %>%
                     set_names("tasmax", "tasmin"),
                   c("average_temperature", "precipitation"),
                   c("vapor_pressure_deficit"))


# MODEL LOOP ------------------------------------------------------------------

model <- models[1]



# DOWNLOAD DATA --------------------------------------------------------------

var_combo <- var_combos[[1]]
  


print(str_glue("Downloading raw data"))

# directory to download raw input data
dir_data <- str_glue("{dir_disk}/data")
fs::dir_create(dir_data)

# get directories of variables to use
dir_vars <-
  dir_gs %>%
  fs::dir_ls(regexp = model) %>%
  fs::dir_ls(regexp = var_combo %>% str_flatten("|"))

# download all files in directories (parallel)
dir_vars %>%
  map(fs::dir_ls) %>%
  unlist(use.names = F) %>%
  str_replace("/mnt/bucket_mine", "gs://clim_data_reg_useast1") %>%
  
  future_walk(function(f) {
    
    system(str_glue("gsutil cp {f} {d}", d = dir_data),
           ignore.stdout = T, ignore.stderr = T)
    
  })



## GENERATE TILES -------------------------------------------------------------


tiles <- 
  dir_data %>% 
  fs::dir_ls() %>% 
  first() %>% 
  read_ncdf(proxy = T) %>% 
  fn_tile(25)

land <- 
  dir_data %>% 
  fs::dir_ls() %>% 
  first() %>% 
  read_ncdf(proxy = F,
            ncsub = cbind(start = c(1,1,1),
                          count = c(NA,NA,1))) %>% 
  adrop() %>% 
  drop_units() %>% 
  setNames("v") %>% 
  mutate(v = if_else(is.na(v), NA, 1))
  
tiles_land <- 
  tiles %>% 
  fn_tile_land(land)



## PREP TIME ------------------------------------------------------------------

# All dates
time_vector <- 
  dir_data %>% 
  fs::dir_ls(regexp = names(var_combo)[1]) %>% 
  future_map(read_ncdf, proxy = T) %>% 
  suppressMessages() %>% 
  map(st_get_dimension_values, "time") %>% 
  map(as.character) %>% 
  unlist() %>% 
  str_sub(end = 10)

# Obtain calendar type
max_feb <- 
  time_vector[str_sub(time_vector, 6,7) == "02"] %>% # filter feb months
  str_sub(9,10) %>% # extract days
  as.numeric() %>% 
  max()

model_cal <- 
  case_when(max_feb == 30 ~ "360_day",
            max_feb == 29 ~ "gregorian",
            max_feb == 28 ~ "noleap")

print(str_glue("   Calendar: {model_cal} --- Max feb: {max_feb}"))

# update time_vector
time_vector <- PCICt::as.PCICt(time_vector, cal = model_cal)



# table for seasonal and dekadal splits
tb_time <- 
  tibble(date = time_vector,
         yr = str_sub(date, end = 4) %>% as.numeric(),
         mon = str_sub(date, 6,7) %>% as.numeric(),
         dy = str_sub(date, 9,10) %>% as.numeric(),
         seas = case_when(mon %in% c(12, 1, 2) ~ "DJF",
                          mon %in% c(3:5) ~ "MAM",
                          mon %in% c(6:8) ~ "JJA",
                          mon %in% c(9:11) ~ "SON"),
         yr_shift = if_else(mon %in% c(1, 2), yr - 1, yr),
         dek = case_when(dy <= 10 ~ 1,
                         dy <= 20 ~ 2,
                         TRUE ~ 3))


seas <- 
  str_glue("{tb_time$yr_shift}_{tb_time$seas}") %>% 
  as.vector()

time_seas <- 
  tb_time %>% 
  group_by(yr_shift, seas) %>% 
  mutate(s = str_glue("{first(yr_shift)}-{first(mon)}-{first(dy)}")) %>% 
  pull(s) %>% 
  as.vector() %>% 
  unique() %>% 
  PCICt::as.PCICt(cal = "gregorian")




## TILES LOOP -----------------------------------------------------------------

# Create dir where all resulting tiles will be saved
dir_tiles <- str_glue("{dir_disk}/tiles")
fs::dir_create(dir_tiles)

pwalk(tiles_land %>% filter(cover == T), function(tile_id, ...) {
  
  
  # print(str_glue(" "))
  # print(str_glue("PROCESSING TILE {tile_i} / {nrow(tiles_wcover)}"))
  
  
  nc <- 
    cbind(start = c(tiles$start_x[tile_id],
                    tiles$start_y[tile_id],
                    1),
          count = c(tiles$end_x[tile_id] - tiles$start_x[tile_id] + 1,
                    tiles$end_y[tile_id] - tiles$start_y[tile_id] + 1,
                    NA))
  
  
  l_s <- 
    names(var_combo) %>%
    set_names() %>% 
    map(function(variable){
      
      s <- 
        dir_data %>% 
        fs::dir_ls(regexp = variable) %>% 
        future_map(read_ncdf, ncsub = nc) %>%                                   # if not at upper level
        suppressMessages() %>% 
        do.call(c, .) %>% 
        setNames("v")
      
      if (str_detect(variable, "tas")) {
        
        s %>% 
          mutate(v = 
                   v %>% 
                   set_units(degC) %>% 
                   set_units(NULL)) %>% 
          setNames(variable)
        
      } else if (variable == "pr") {
        
        s %>% 
          mutate(v =
                   v %>% 
                   set_units(kg/m^2/d) %>% 
                   set_units(NULL) %>% 
                   if_else(. < 0, 0, .)) %>% 
          setNames(variable)
     
      } 
    })
    
  
  
  
  
  # CALCULATE INDICATORS
  
  
  
  if (var_combo[1] == "maximum")) {
  
  
  
  # Maximum num of consecutive frost days 
  
  s <- 
    l_s %>% 
    pluck("tasmin") %>% 
    st_apply(c(1,2), function(x) {
      
      if (any(is.na(x))) {
        rep(NA, length(time_seas))
      } else {
        
        aggregate(x, by = list(seas), function(xx) {
          
          r <- rle(xx < 0)
          l <- r$lengths[r$values == T]
          
          if (length(l) < 1) {
            0
          } else {
            max(l)
          }
          
        })$x
        
      }
      
    },
    FUTURE = T,
    .fname = "time") %>% 
    aperm(c(2,3,1)) %>% 
    st_set_dimensions(3, values = time_seas)
  
  saveRDS(s, str_glue("{dir_tiles}/cfd_{tile_id}.rds"))
  
  
  # Cold-spell duration index
  bl_i <- which(year(time_vector) == 1970) %>% first()
  bl_f <- which(year(time_vector) == 2000) %>% last()
  zoo::rollmean(x, 5, na.pad = T, align = "center")[bl_i:bl_f] %>% 
    quantile(c(0.1))
  
  
  
    
  }
  
  
  
  
  
})

