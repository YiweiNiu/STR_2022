library(viridis)
library(ggsci)
library(RColorBrewer)

# label 2 color
label2Col <- function(labels) {
  if (!is.factor(labels)) labels <- factor(labels)
  num <- length(levels(labels))
  colors <-
    if (num <= 9) {
      RColorBrewer::brewer.pal(num, "Set1")
    } else {
      gg_color_hue(num)
    }
  colors[labels]
}

# more colors (ggplot)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  grDevices::hcl(h=hues, l=65, c=100)[1:n]
}

period_color_map = c(
  "Di-" = "#79AF97FF",
  "Tri-" = "#EFC000FF",
  "Tet-" = "#7AA6DCFF",
  "Pen-" = "#f0b98d",
  "Hex-" = "#CD534CFF"
)

pop_color_maps = c(
  "Unassigned" = "#d6d6d6",
  "1KGP" = "#4A6990FF",
  "NyuWa" = "#003C67FF",
  "AFR" = "#0073C2FF",
  "ACB" = "#023fa5", "ASW" = "#7d87b9", "ESN" = "#bec1d4", "GWD" = "#d6bcc0",
  "LWK" = "#bb7784", "MSL" = "#8e063b", "YRI" = "#4a6fe3",
  "AMR" = "#EFC000FF",
  "CLM" = "#8595e1", "MXL" = "#b5bbe3", "PEL" = "#e6afb9", "PUR" = "#e07b91",
  "EAS" = "#7AA6DCFF",
  "CHN.NyuWa" = "#d33f6a", "CHS.NyuWa" = "#11c638",
  "CHB.1KGP" = "#8dd593", "CHS.1KGP" = "#c6dec7", "CDX" = "#ead3c6",
  "JPT" = "#f0b98d", "KHV" = "#ef9708",
  "EUR" = "#868686FF",
  "CEU" = "#0fcfc0", "FIN" = "#9cded6", "GBR" = "#d5eae7",
  "IBS" = "#f3e1eb", "TSI" = "#f6c4e1",
  "SAS" = "#CD534CFF",
  "BEB" = "#f79cd4", "GIH" = "#7f7f7f", "ITU" = "#c7c7c7",
  "PJL" = "#1CE6FF", "STU" = "#336600"
)

spop_shape_map = c(
  "AFR" = 7,
  "AMR" = 19,
  "EUR" = 15,
  "EAS" = 11,
  "SAS" = 17
)

region_color_maps = c(
  "Central China" = "#d62728",
  "East China" = "#aa40fc",
  "North China" = "#8c564b",
  "Northeast China" = "#e377c2",
  "Northwest China" = "#7f7f7f",
  "South China" = "#b5bd61",
  "Southwest China" = "#17becf"
)

# three colors for heatmap
c("#476fa9", "#ffffff", "#ca3226")
# or
c("#3658d3", "#ffffff", "#900504")

# negative and positive
colors_2 = c("#476fa9", "#ca3226")

# Un, neg, pos
colors_3 = c("#d6d6d6", "#476fa9", "#ca3226")

# two color (mild)
c("#7da1d5", "#bdde85")



