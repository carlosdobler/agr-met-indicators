
library(tidyverse)
library(stars)
library(units)
library(furrr)
library(climdex.pcic)

options(future.fork.enable = T)
plan(multicore)
options(future.globals.maxSize = 1000*1024^2)

source("https://raw.github.com/carlosdobler/spatial-routines/master/tile.R")
source("https://raw.github.com/carlosdobler/spatial-routines/master/tile_land.R")
source("https://raw.github.com/carlosdobler/spatial-routines/master/write_nc.R")


# cloud directory of raw input data
dir_gs <- "/mnt/bucket_mine/cmip6/nex/daily"

# local disk where temp files will be stored
dir_disk <- "/mnt/pers_disk"

# local directory to download raw input data
dir_data <- str_glue("{dir_disk}/data")

# Create dir where all resulting tiles will be saved
dir_tiles <- str_glue("{dir_disk}/tiles")

# cloud directory where final results will be uploaded
dir_res <- "gs://clim_data_reg_useast1/cmip6/nex"


models <- c("GFDL-ESM4",
            "MPI-ESM1-2-HR",
            "MRI-ESM2-0",
            "UKESM1-0-LL",
            "IPSL-CM6A-LR")


var_combos <- list(c("maximum_temperature", "minimum_temperature") %>%
                     set_names("tasmax", "tasmin"),
                   c("average_temperature", "precipitation") %>% 
                     set_names("tas", "pr"),
                   c("vapor_pressure_deficit"))


tb_vars <- 
  read_csv("/mnt/bucket_mine/misc_data/tb_vars.csv")

tb_vars$var_short[9] <- "bedd-10"



fn_calc_indicator <- function(s, fn, v) {
  
  s %>% 
    st_apply(c(1,2), function(x) {
      
      if (any(is.na(x))) {
        
        rep(NA, length(time_seas))
        
      } else {
        
        fn(x)
        
      }
      
    },
    .fname = "time") %>% 
    aperm(c(2,3,1)) %>% 
    st_set_dimensions(3, values = time_seas) %>% 
    
    write_nc(time_seas, str_glue("{dir_tiles}/{v}_{tile_id}.nc"))
  
}



# MODEL LOOP ------------------------------------------------------------------

for (model in models[5]) {
  
  print(str_glue("PROCESSING MODEL {model}"))
  
  
  
  # GENERATE TILES -------------------------------------------------------------
  
  print(str_glue("    Tiling..."))
  
  # dir_bucket_gs <- "gs://clim_data_reg_useast1" # hide
  # dir_nex_gs <- str_glue("{dir_bucket_gs}/cmip6/nex")
  
  
  s_proxy <- 
    
    # str_glue("gcloud storage ls {dir_nex_gs}/daily/{model}/precipitation/*.nc") %>% 
    # system(intern = T) %>% 
    # first() %>% 
    # str_replace("gs:/", "/vsigs") %>% 
    # read_mdim(count = c(NA,NA,1)) %>% 
    # adrop()
    
    dir_gs %>% 
    fs::dir_ls(regexp = model) %>%
    fs::dir_ls(regexp = "precipitation") %>% # any variable
    fs::dir_ls() %>% 
    first() %>% 
    read_ncdf(proxy = F,
              ncsub = cbind(start = c(1,1,1),
                            count = c(NA,NA,1))) %>% 
    suppressMessages() %>% 
    adrop()
  
  tiles <- 
    s_proxy %>% 
    fn_tile(25)
  
  land <- 
    s_proxy %>% 
    drop_units() %>% 
    setNames("v") %>% 
    mutate(v = if_else(is.na(v), NA, 1))
  
  tiles_land <- 
    tiles %>% 
    fn_tile_land(land)
  
  
  
  # VAR COMBO LOOP --------------------------------------------------------------
  
  for (vc in seq_along(var_combos[1:2])) {
    
    var_combo <- var_combos[[vc]]
    
    
    print(str_glue("  PROCESSING INDs W VARS {str_flatten(var_combo, ' + ')}"))
    
    tb_vars_sub <- 
      tb_vars %>% 
      filter(under_var %in% var_combo | str_detect(under_var, var_combo[1]))
    
    
    # DOWNLOAD DATA --------------------------------------------------------------
    
    
    print(str_glue("    Downloading raw data..."))
    
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
    
    
    
    
    # PREP TIME ------------------------------------------------------------------
    
    print(str_glue("    Prep time vectors..."))
    
    # All dates
    time_vector <- 
      dir_data %>% 
      fs::dir_ls(regexp = str_glue("{names(var_combo)[1]}_")) %>%
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
    
    
    
    # table for seasonal splits
    tb_time <- 
      tibble(date = time_vector,
             yr = str_sub(date, end = 4) %>% as.numeric(),
             mon = str_sub(date, 6,7) %>% as.numeric(),
             dy = str_sub(date, 9,10) %>% as.numeric(),
             seas = case_when(mon %in% c(12, 1, 2) ~ "1", # DJF
                              mon %in% c(3:5) ~ "2", # MAM
                              mon %in% c(6:8) ~ "3", # JJA
                              mon %in% c(9:11) ~ "4"), # SON
             yr_shift = if_else(mon %in% c(1, 2), yr - 1, yr)) %>% 
      group_by(yr) %>% 
      mutate(jdy = row_number()) %>% 
      ungroup()
    
    
    seas <- 
      str_glue("{tb_time$yr_shift}_{tb_time$seas}") %>% 
      as.vector() %>% 
      as.factor()
    
    time_seas <- 
      tb_time %>% 
      group_by(yr_shift, seas) %>% 
      mutate(s = str_glue("{first(yr_shift)}-{first(mon)}-{first(dy)}")) %>% 
      pull(s) %>% 
      as.vector() %>% 
      unique() %>% 
      PCICt::as.PCICt(cal = "gregorian")
    
    first_day_sh <- which(month(time_vector) == 7) %>% first() %>% {.-1}
    
    time_vector_sh <- head(time_vector, -first_day_sh)
    
    
    time_yrs <- 
      time_vector %>% 
      year() %>% 
      unique() %>% 
      {str_glue("{.}-01-01")} %>% 
      PCICt::as.PCICt(cal = "gregorian")
    
    
    
    
    # TILES LOOP -----------------------------------------------------------------
    
    fs::dir_create(dir_tiles)
    
    tiles_wland <- tiles_land %>% filter(cover == T)
    
    
    for (tilerow in 86:nrow(tiles_wland)) {                                 # ******************
      
      print(str_glue(" "))
      print(str_glue("    Calculating tile {tilerow} / {nrow(tiles_wland)}"))
      
      tile_id <- tiles_wland$tile_id[tilerow]
      
      
      # Load data
      
      print(str_glue("      Loading data..."))
      tictoc::tic("      done")
      
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
            future_map(read_ncdf, ncsub = nc) %>%
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
      
      tictoc::toc()
      
      
      
      # CALCULATE INDICATORS -------------------------------------------------------
      
      print(str_glue("      Calculating indices..."))
      
      future_walk(tb_vars_sub$var_short, function(var_sh) {
        
        # print(var_sh)
        
        # *************
        
        if (var_sh == "cfd") {
          
          fn_calc_indicator(l_s$tasmin,
                            function(x) {
                              spell.length.max(x,
                                               date.factor = seas,
                                               threshold = 0,
                                               op = "<",
                                               spells.can.span.years = F)
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "csdi") {
          
          fn_calc_indicator(l_s$tasmin,
                            function(x) {
                              
                              q <- get.outofbase.quantiles(tmin = x,
                                                           tmin.dates = time_vector,
                                                           base.range = c(1971,2000),
                                                           temp.qtiles = 0.1,
                                                           n = 5) %>%
                                unlist() %>%
                                unname()
                              
                              threshold.exceedance.duration.index(daily.temp = x,
                                                                  date.factor = seas,
                                                                  jdays = tb_time$jdy,
                                                                  thresholds = q,
                                                                  op = "<",
                                                                  spells.can.span.years = F,
                                                                  max.missing.days = 1000)
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "wsdi") {
          
          fn_calc_indicator(l_s$tasmax,
                            function(x) {
                              
                              q <- get.outofbase.quantiles(tmax = x,
                                                           tmax.dates = time_vector,
                                                           base.range = c(1971,2000),
                                                           temp.qtiles = 0.9,
                                                           n = 5) %>%
                                unlist() %>%
                                unname()
                              
                              threshold.exceedance.duration.index(daily.temp = x,
                                                                  date.factor = seas,
                                                                  jdays = tb_time$jdy,
                                                                  thresholds = q,
                                                                  op = ">",
                                                                  spells.can.span.years = F,
                                                                  max.missing.days = 1000)
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "csu") {
          
          fn_calc_indicator(l_s$tasmax,
                            function(x) {
                              spell.length.max(x,
                                               date.factor = seas,
                                               threshold = 25,
                                               op = ">",
                                               spells.can.span.years = F)
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "dtr") {
          
          fn_calc_indicator(l_s$tasmax - l_s$tasmin,
                            function(x) {
                              aggregate(x, by = list(seas), mean)$x
                            },
                            var_sh)
        }
        
        
        # *************
        
        if (var_sh == "fd") {
          
          fn_calc_indicator(l_s$tasmin,
                            function(x) {
                              aggregate(x < 0, by = list(seas), sum)$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "id") {
          
          fn_calc_indicator(l_s$tasmax,
                            function(x) {
                              aggregate(x < 0, by = list(seas), sum)$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "su") {
          
          fn_calc_indicator(l_s$tasmax,
                            function(x) {
                              aggregate(x > 25, by = list(seas), sum)$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "tn") {
          
          fn_calc_indicator(l_s$tasmin,
                            function(x) {
                              aggregate(x, by = list(seas), mean)$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "tnn") {
          
          fn_calc_indicator(l_s$tasmin,
                            function(x) {
                              aggregate(x, by = list(seas), min)$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "tnx") {
          
          fn_calc_indicator(l_s$tasmin,
                            function(x) {
                              aggregate(x, by = list(seas), max)$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "tr") {
          
          fn_calc_indicator(l_s$tasmin,
                            function(x) {
                              aggregate(x > 20, by = list(seas), sum)$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "tx") {
          
          fn_calc_indicator(l_s$tasmax,
                            function(x) {
                              aggregate(x, by = list(seas), mean)$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "txn") {
          
          fn_calc_indicator(l_s$tasmax,
                            function(x) {
                              aggregate(x, by = list(seas), min)$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "txx") {
          
          fn_calc_indicator(l_s$tasmax,
                            function(x) {
                              aggregate(x, by = list(seas), max)$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "cdd") {
          
          fn_calc_indicator(l_s$pr,
                            function(x) {
                              spell.length.max(x,
                                               date.factor = seas,
                                               threshold = 1,
                                               op = "<",
                                               spells.can.span.years = F)
                            },
                            var_sh)
        } 
        
        # *************
        
        if (var_sh == "cwd") {
          
          fn_calc_indicator(l_s$pr,
                            function(x) {
                              spell.length.max(x,
                                               date.factor = seas,
                                               threshold = 1,
                                               op = ">=",
                                               spells.can.span.years = F)
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "ww") {
          
          fn_calc_indicator(c(l_s$tas, l_s$pr, along = "v"),
                            function(x) {
                              
                              q_tas <- 
                                get.outofbase.quantiles(tmin = x[,1],
                                                        tmin.dates = time_vector,
                                                        base.range = c(1971,2000),
                                                        temp.qtiles = 0.75,
                                                        n = 5) %>%
                                unlist() %>%
                                unname()
                              
                              
                              t_index <- year(time_vector) >= 1971 & year(time_vector) <= 2000
                              q_pr <-
                                aggregate(x[,2][t_index], 
                                          by = list(yday(time_vector[t_index])), 
                                          function(x) quantile(x, 0.75))$x
                              
                              
                              tibble(tas = x[,1],
                                     pr = x[,2],
                                     time = time_vector,
                                     jdy = yday(time_vector)) %>% 
                                left_join(tibble(q_tas = q_tas) %>% mutate(jdy = row_number())) %>% 
                                left_join(tibble(q_pr = q_pr) %>% mutate(jdy = row_number())) %>% 
                                suppressMessages() %>% 
                                mutate(ww = if_else(tas > q_tas & pr >= 1 & pr > q_pr, 1, 0)) %>% 
                                pull(ww) %>% 
                                
                                aggregate(by = list(seas), sum) %>% 
                                .$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "bedd") {
          
          fn_calc_indicator(l_s$tas,
                            function(x) {
                              aggregate(x > 10 & x < 30, by = list(seas), sum)$x
                            },
                            "bedd-10")
        }
        
        # *************
        
        if (var_sh == "gsl") {
          
          nh <- st_get_dimension_values(l_s$tas, "lat")[1] >= 0
          
          l_s$tas %>% 
            st_apply(c(1,2), function(x) {
              
              if (any(is.na(x))) {
                
                rep(NA, length(unique(year(time_vector))))
                
              } else {
                
                if (nh == TRUE) {
                  
                  growing.season.length(x, 
                                        factor(year(time_vector)), 
                                        time_vector,
                                        northern.hemisphere = T)
                  
                } else {
                  
                  growing.season.length(tail(x, -first_day_sh), 
                                        factor(year(time_vector_sh)), 
                                        time_vector_sh,
                                        northern.hemisphere = T)
                  
                }
              }
            }, 
            .fname = "time") %>% 
            aperm(c(2,3,1)) %>% 
            st_set_dimensions(3, values = time_yrs) %>% 
            write_nc(time_yrs, str_glue("{dir_tiles}/{var_sh}_{tile_id}.nc"))
        }
        
        # *************
        
        if (var_sh == "r10mm") {
          
          fn_calc_indicator(l_s$pr,
                            function(x) {
                              aggregate(x >= 10, by = list(seas), sum)$x
                            },
                            var_sh)
        }   
        
        # *************
        
        if (var_sh == "r20mm") {
          
          fn_calc_indicator(l_s$pr,
                            function(x) {
                              aggregate(x >= 20, by = list(seas), sum)$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "rr") {
          
          fn_calc_indicator(l_s$pr,
                            function(x) {
                              aggregate(x, by = list(seas), sum)$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "rr1") {
          
          fn_calc_indicator(l_s$pr,
                            function(x) {
                              aggregate(x >= 1, by = list(seas), sum)$x
                            },
                            var_sh)
        }
        
        # *************
        
        if (var_sh == "sdii") {
          
          fn_calc_indicator(l_s$pr,
                            function(x) {
                              aggregate(x, 
                                        by = list(seas), 
                                        function(x) {
                                          if (all(x < 1)) {
                                            0
                                          } else {
                                            mean(x[x >= 1])
                                          }
                                          
                                        })$x
                            },
                            var_sh)
          
          # simple.precipitation.intensity.index(x, seas)
        }
        
        # *************
        
        if (var_sh == "tg") {
          
          fn_calc_indicator(l_s$tas,
                            function(x) {
                              aggregate(x, by = list(seas), mean)$x
                            },
                            var_sh)
        }
        
      })
      
      
    } # End of tiles loop
    
    
    
    
    # MOSAIC ----------------------------------------------------------------------
    
    
    vars <- 
      dir_tiles %>% 
      fs::dir_ls() %>% 
      fs::path_ext_remove() %>% 
      fs::path_file() %>% 
      str_split("_", simplify = T) %>% 
      .[,1] %>% 
      unique()
    
    
    for (var in vars) {
      
      print(str_glue(" "))
      print(str_glue("    Mosaicking var: {var}"))
      
      
      if (var == "gsl") {
        # annual (as opposed to seasonal)
        yrs_i <- cut(year(time_yrs), seq(1960, 2100, by = 10), right = F, labels = F)
      } else {
        yrs_i <- cut(year(time_seas), c(1959, seq(1960, 2100, by = 10)), right = F, labels = F)
      }
      
      
      
      # loop through blocks of yrs
      for (yr_i in unique(yrs_i)) {
        
        print(str_glue("      year block: {yr_i}"))
        
        pos <- 
          which(yrs_i == yr_i)
        
        if (var == "gsl") {
          time_ <- time_yrs[pos]
        } else {
          time_ <- time_seas[pos]
        }
        
        
        cols_1var <-
          
          # loop through cols
          future_map(unique(tiles$start_x), function(col_i) {
            
            r <- 
              tiles %>% 
              filter(start_x == col_i) %>% 
              left_join(st_drop_geometry(tiles_land), by = "tile_id") %>% 
              pmap(function(tile_id, start_x, end_x, start_y, end_y, cover, ...) {
                
                if (cover == FALSE) {
                  
                  s <- s_proxy[, start_x:end_x, start_y:end_y]
                  
                  empty_array <- 
                    array(NA, dim = c(dim(s)[1], 
                                      dim(s)[2],
                                      last(pos) - first(pos) +1))
                  
                  names(dim(empty_array)) <- c("lon", "lat", "time")
                  
                  # format as stars
                  empty_stars <- 
                    st_as_stars(empty_array)
                  
                  st_dimensions(empty_stars)[1] <- st_dimensions(s)[1]
                  st_dimensions(empty_stars)[2] <- st_dimensions(s)[2]
                  empty_stars <- st_set_dimensions(empty_stars, 3, values = time_)
                  
                  return(empty_stars)
                  
                  
                } else {
                  
                  dir_tiles %>% 
                    fs::dir_ls(regexp = str_glue("/{var}_{tile_id}.nc")) %>%
                    read_ncdf(ncsub = cbind(start = c(1, 1, first(pos)),
                                            count = c(NA,NA,last(pos) - first(pos) +1))) %>% 
                    suppressMessages()
                }
                
                
                
              })
            
            r <- 
              do.call(c, c(r, along = 2))
            
            st_dimensions(r)[2] <- st_dimensions(s_proxy)[2]
            
            return(r)
            
          })
        
        
        # merge all columns
        mos <- 
          do.call(c, c(cols_1var, along = 1))
        
        # fix lon dimension
        st_dimensions(mos)[1] <- st_dimensions(s_proxy)[1]
        
        
        # SAVING ------------------
        
        if (var == "gsl") {
          per <- "yr" %>% set_names("annual")
        } else {
          per <- "seas" %>% set_names("seasonal")
        }
        
        
        filename <-
          str_glue("{tempdir()}/{var}_{per}_{model}_{year(first(time_))}-{year(last(time_))}.nc")
        
        un <- tb_vars %>% filter(var_short == var) %>% pull(unit)
        var_long <- tb_vars %>% filter(var_short == var) %>% pull(var_long)
        
        write_nc(mos %>% 
                   setNames("v") %>% 
                   mutate(v = units::set_units(v, !!un)) %>% 
                   setNames(var),
                 time_seas[pos],
                 filename)
        
        filename_gs <- str_glue("{dir_res}/{names(per)}/{model}/{var_long}/{fs::path_file(filename)}")
        
        str_glue("gsutil cp {filename} {filename_gs}") %>% 
          system(ignore.stdout = T, ignore.stderr = T)
        
        fs::file_delete(filename)
        
        
      } # End of block of years loop
      
      # remove tiles
      # dir_tiles %>%
      #   fs::dir_ls() %>% #regexp = str_glue("{var}_")) %>%
      #   walk(fs::file_delete)
      
    } # End of vars loop
    
    fs::dir_delete(dir_data)
    fs::dir_delete(dir_tiles)
    
  }
  
}









