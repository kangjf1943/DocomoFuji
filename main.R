# Statement ----
# Calculate consumer surplus of trailheads of Fujisan. Compare results of different resolutions. 

# Preparation ----
pacman::p_load(haven, dplyr, sf, tmap, jpmesh, purrr, lubridate, showtext, ggplot2)
showtext_auto()

# Data ----
# Trailhead points. 
trailhead <- matrix(
  c("yoshida", 35.39433, 138.73259,
    "subashiri", 35.36593, 138.77845,
    "gotemba", 35.3358, 138.79465,
    "fujinomiya", 35.3366, 138.73346
  ),
  ncol = 3, byrow = TRUE
) %>%
  as.data.frame() %>%
  rename_with(~ c("head", "lat", "long")) %>%
  mutate(
    lat = as.numeric(lat), 
    long = as.numeric(long), 
    mesh = as.character(coords_to_mesh(longitude = long, latitude = lat))
  ) %>% 
  st_as_sf(coords = c("long", "lat")) %>%
  st_set_crs(4326) %>%
  st_transform(6668)

# Prefecture names and codes. 
prefcode <- read.csv("data_raw/prefcode_citycode_master_UTF-8.csv") %>%
  tibble() %>% 
  select(prefcode, prefname) %>%
  mutate(prefcode = as.character(prefcode)) %>% 
  distinct()

# Distance from prefectures to trailheads. Inherited from Uryu's calculation for "22_AgoopFuji" project.
osrm_res <-
  readRDS("data_raw/osrm/prefecture_office_to_4route_entrance.rds") %>%
  st_drop_geometry() %>%
  mutate(head = case_when(
    route == "吉田" ~ "yoshida",
    route == "須走" ~ "subashiri",
    route == "御殿場" ~ "gotemba",
    route == "富士宮" ~ "fujinomiya"
  )) %>% 
  mutate(trv_cost = distance / 21.7 * 126)

# Population of prefectures. 
prefpop <- read.csv(
  "data_raw/FEI_PREF_230406150517.csv", skip = 8
) %>% 
  tibble() %>% 
  select(
    "YEAR", "AREA.Code", "AREA", "A1101_Total.population..Both.sexes..person."
  ) %>% 
  rename_with(~ c("year", "pref_code", "pref", "pop")) %>% 
  mutate(pref_code = pref_code / 1000) %>% 
  mutate(pref_code = sprintf("%02d", pref_code)) %>% 
  filter(year == 2016) %>% 
  mutate(
    pref_code = as.character(as.numeric(pref_code)), 
    pop = as.numeric(gsub(",", "", pop))
  )

# Docomo raw data. 
docomo_all <- read_dta("data_raw/Docomo_Fuji/02_KyotoUniv_pref.dta") %>% 
  rename(mesh = area, prefcode = residence, vis_pop = population) %>% 
  mutate(mesh = as.character(mesh), prefcode = as.character(prefcode)) 

# The meshes covered by DOCOMO data. 
docomo_all_mesh <- docomo_all %>% 
  select(mesh) %>% 
  distinct() %>% 
  mutate(rowid = row_number()) %>% 
  group_by(rowid) %>% 
  mutate(mesh = as_meshcode(as.character(mesh))) %>% 
  meshcode_sf(., mesh) %>% 
  st_union()

# Docomo data of trailheads. 
docomo_head <- docomo_all %>% 
  left_join(trailhead, by = "mesh") %>% 
  filter(!is.na(head)) %>% 
  mutate(
    date = as_date(as.character(date)), 
    year = year(date), month = month(date), day = day(date)
  )

# Function to get raw data for travel cost model regression. 
get_reg_raw <- function(x, res, vis_pop_method) {
  if (vis_pop_method == "sum") {
    x_grp <- x %>% 
      group_by(date, get(res), head, mesh, prefcode) %>% 
      summarise(vis_pop = sum(vis_pop), .groups = "drop")
  } else if (vis_pop_method == "max") {
    x_grp <- x %>% 
      group_by(date, get(res), head, mesh, prefcode) %>% 
      summarise(vis_pop = max(vis_pop), .groups = "drop")
  }
  x_grp %>% 
    group_by(`get(res)`, head, mesh, prefcode) %>% 
    summarise(vis_pop = sum(vis_pop), .groups = "drop") %>% 
    # Add distance. 
    left_join(prefcode, by = "prefcode") %>% 
    left_join(osrm_res, by = c("head", "prefname" = "src_loc")) %>% 
    left_join(prefpop, by = c("prefcode" = "pref_code")) %>% 
    mutate(vis_rate = vis_pop / pop) %>% 
    rename("res" = "get(res)") %>% 
    return()
}

# Function to apply travel cost model regression. 
calc_val <- function(x) {
  group_by(x, head) %>% 
    summarise(
      demand_model = list(summary(lm(log(vis_rate) ~ trv_cost))), 
      all_coef = lapply(demand_model, coef), 
      tc_coef = lapply(all_coef, function(x) x[2, 1]) %>% do.call("c", .), 
      p = lapply(all_coef, function(x) x[2, 4]) %>% do.call("c", .), 
      avg_cs = -1 / tc_coef, 
      # Bug: Am I right about that? Since we used visiting rate (unit: trips per capita) in demand curve, so the unit for the consumer surplus is JPY per trip per capita. According to Kubo et al's beach research, the total consumer surplus can be calculated by multiplying this consumer surplus by visitor number of each visitor center.
      vis_pop = sum(vis_pop), 
      tot_cs = avg_cs * vis_pop, 
      sample_size = n()
    ) %>% 
    return()
}

# Analysis ----
# DOCOMO mesh does not cover all trailheads. 
tm_shape(docomo_all_mesh) + 
  tm_polygons() + 
  tm_shape(trailhead) + 
  tm_dots(size = 0.5, col = "red4") + 
  tm_text(text = "head", just = "right")

# Data for travel cost model regression. 
reg_raw <- rbind(
  # Results by "max" method. 
  rbind(
    get_reg_raw(docomo_head, "year", "max") %>% 
      mutate(res = "year", size = 3), 
    get_reg_raw(docomo_head, "month", "max") %>% 
      mutate(res = "month", size = 1.5), 
    get_reg_raw(docomo_head, "date", "max") %>% 
      mutate(res = "day", size = 1)
  ) %>% 
    mutate(vis_pop_method = "max"), 
  # Results by "sum" method. 
  rbind(
    get_reg_raw(docomo_head, "year", "sum") %>% 
      mutate(res = "year", size = 3), 
    get_reg_raw(docomo_head, "month", "sum") %>% 
      mutate(res = "month", size = 1.5), 
    get_reg_raw(docomo_head, "date", "sum") %>% 
      mutate(res = "day", size = 1)
  ) %>% 
    mutate(vis_pop_method = "sum") 
)

# Compare relationship between visiting rate and travel cost under different time scales and visiting population counting method. 
ggplot(reg_raw, aes(trv_cost, log(vis_rate))) + 
  geom_point(aes(col = res, size = size), alpha = 0.9) + 
  geom_smooth(aes(col = res), method = "lm", se = FALSE, alpha = 0.5) + 
  facet_grid(vis_pop_method ~ head) + 
  theme_bw() + 
  labs(x = "Travel cost", y = "Ln(visiting rate)")

# Apply travel cost model. 
head_val <- map2(
  rep(c("year", "month", "day"), 2), rep(c("max", "sum"), each = 3), 
  function(x, y) {
    filter(reg_raw, res == x, vis_pop_method == y) %>% 
      calc_val() %>% 
      mutate(res = x, method = y)
  }
) %>% 
  bind_rows()

# Compare values and population under different time scales and visiting population counting methods. 
ggplot(head_val) + 
  geom_point(aes(head, avg_cs, col = res)) + 
  facet_grid(.~ method) + 
  theme_bw() + 
  labs(x = "", y = "Per trip CS")
ggplot(head_val) + 
  geom_point(aes(head, tot_cs, col = res)) + 
  facet_grid(.~ method) + 
  theme_bw() + 
  labs(x = "", y = "Total CS")
ggplot(head_val) + 
  geom_point(aes(head, vis_pop, col = res)) + 
  facet_grid(.~ method) + 
  theme_bw() + 
  labs(x = "", y = "Visiting population")

# Box plot for daily travel cost. 
reg_raw %>% 
  filter(res == "day") %>% 
  mutate(tot_trv_cost = trv_cost * vis_pop) %>% 
  ggplot() + 
  geom_point()

# Box plot for daily total CS. 
# Use day*max per trip CS. 
docomo_head %>% 
  group_by(head, date) %>% 
  summarise(vis_pop = max(vis_pop), .groups = "drop") %>% 
  left_join(
    filter(head_val, res == "day", method == "max") %>% 
      select(head, avg_cs), 
    by = "head"
  ) %>% 
  mutate(tot_cs = avg_cs * vis_pop) %>% 
  ggplot(aes(head, tot_cs)) + 
  geom_boxplot() + 
  geom_jitter()

# Use day*sum per trip CS. 
docomo_head %>% 
  group_by(head, date) %>% 
  summarise(vis_pop = sum(vis_pop), .groups = "drop") %>% 
  left_join(
    filter(head_val, res == "day", method == "sum") %>% 
      select(head, avg_cs), 
    by = "head"
  ) %>% 
  mutate(tot_cs = avg_cs * vis_pop) %>% 
  ggplot(aes(head, tot_cs)) + 
  geom_boxplot() + 
  geom_jitter()

# Get daily travel cost ----
# Holiday data. 
holiday <- read.csv("data_raw/Japan_holiday_2016.csv") %>%
  mutate(
    year = 2016, month = substr(.$月日, 1, 2), day  = substr(.$月日, 4, 5)
  ) %>% 
  rename(holiday_name = 名称) %>%
  mutate(date = as.Date(paste(year, month, day, sep = "-"))) %>%
  dplyr::select(date, holiday_name)

# Function to get boxplot raw data. 
get_box_raw <- function(x, vis_pop_method) {
  if (vis_pop_method == "sum") {
    x_grp <- x %>% 
      group_by(date, month, head, mesh, prefcode) %>% 
      summarise(vis_pop = sum(vis_pop), .groups = "drop")
  } else if (vis_pop_method == "max") {
    x_grp <- x %>% 
      group_by(date, month, head, mesh, prefcode) %>% 
      summarise(vis_pop = max(vis_pop), .groups = "drop")
  }
  x_grp %>% 
    group_by(date, month, head, mesh, prefcode) %>% 
    # Add distance. 
    left_join(prefcode, by = "prefcode") %>% 
    left_join(osrm_res, by = c("head", "prefname" = "src_loc")) %>% 
    left_join(prefpop, by = c("prefcode" = "pref_code")) %>% 
    left_join(holiday, by = "date") %>% 
    group_by(date, month, head, holiday_name) %>% 
    summarise(trv_cost = sum(trv_cost)) %>% 
    return()
}
# Method by "sum" and by "max", color by month or weekday. 
library(patchwork)

(
  (
    get_box_raw(docomo_head, "sum") %>% 
      ggplot(aes(head, trv_cost)) + 
      geom_boxplot() + 
      geom_jitter(
        aes(col = as_factor(month)), width = 0.2, height = 100, alpha = 0.5
      ) + 
      labs(col = "Month", x = "", y = "Travel cost", title = "Sum-Month")
  ) | (
    get_box_raw(docomo_head, "sum") %>% 
      mutate(weekday = weekdays(date)) %>% 
      mutate(weekday = case_when(
        weekday == "Saturday" | weekday == "Sunday" | !is.na(holiday_name) ~ 
          "weekend/holiday", 
        TRUE ~ "weekday"
      )) %>% 
      ggplot(aes(head, trv_cost)) + 
      geom_boxplot() + 
      geom_jitter(
        aes(col = weekday), width = 0.2, height = 100, alpha = 0.5
      ) + 
      labs(col = "Weekday", x = "", y = "Travel cost", title = "Sum-Weekday")
  )
) / (
  (
    get_box_raw(docomo_head, "max") %>% 
      ggplot(aes(head, trv_cost)) + 
      geom_boxplot() + 
      geom_jitter(
        aes(col = as_factor(month)), width = 0.2, height = 100, alpha = 0.5
      ) + 
      labs(col = "Month", x = "", y = "Travel cost", title = "Max-Month")
  ) | (
    get_box_raw(docomo_head, "max") %>% 
      mutate(weekday = weekdays(date)) %>% 
      mutate(weekday = case_when(
        weekday == "Saturday" | weekday == "Sunday" | !is.na(holiday_name) ~ 
          "weekend/holiday", 
        TRUE ~ "weekday"
      )) %>% 
      ggplot(aes(head, trv_cost)) + 
      geom_boxplot() + 
      geom_jitter(
        aes(col = weekday), width = 0.2, height = 100, alpha = 0.5
      ) + 
      labs(col = "Weekday", x = "", y = "Travel cost", title = "Max-Weekday")
  )
)

