scCOL =  function (n)
{
  colorSpace <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                  "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4",
                  "#222F75", "#1B9E77", "#B2DF8A", "#E3BE00", "#FB9A99",
                  "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261",
                  "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69",
                  "#8DD3C7", "#999999")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}
