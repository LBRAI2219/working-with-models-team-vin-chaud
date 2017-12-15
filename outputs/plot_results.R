install.packages("tidyverse")
library("tidyverse")

# Remove the "space" in the tab files

# Load the data
path = "/home/brieuc/Documents/GitHub/working-with-models-team-vin-chaud/outputs/"
setwd(path) # Tell R where work

rs = read_tsv("nitrate = 100/tabled_output.tab")
rs$simulation = "100" # Add a column to differentiate the different simulations

rs2 = read_tsv("nitrate_x0.2/tabled_output.tab")
rs2$simulation = "x0.2" # Add a column to differentiate the different simulations

rs3 = read_tsv("nitratex2/tabled_output.tab")
rs3$simulation = "x2" # Add a column to differentiate the different simulations

# Combine the data
rs = rbind(rs, rs2, rs3)

unique(rs$name)


# Plot the data
rs %>%
  filter(name == "plantNutrientUptake") %>%
  filter(path == "//plants/bean/nitrate") %>%
  #filter(time >= 18) %>%
  ggplot(aes(time, value, colour=simulation)) + 
    geom_point() + 
    geom_line()

rs %>%
  filter(name == "concentration") %>%
  filter(path == "//environment/soil/nitrate") %>%
  # filter(time >= 18) %>%
  ggplot(aes(time, value, colour=simulation)) + 
  geom_point() + 
  geom_line()



rs %>%
  filter(name == "rootLength") %>%
  #filter(path == "//plants/bean/nitrate") %>%
  #filter(time >= 18) %>%
  ggplot(aes(time, value, colour=simulation)) + 
  geom_point() + 
  geom_line()
