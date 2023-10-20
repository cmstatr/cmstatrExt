library(tidyverse)

fff_shear <- read_csv("data-raw/fff_shear.csv") %>%
  mutate(Specimen = case_when(
    Specimen == 1 ~ "A",
    Specimen == 2 ~ "B",
    Specimen == 3 ~ "C"
  ))

fff_shear %>%
  ggplot(aes(x = Strain, y = Stress, color = Specimen)) +
  geom_point()

usethis::use_data(fff_shear, overwrite = TRUE)
