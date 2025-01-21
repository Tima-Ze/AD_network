# let's comment this function

wTO_consensus <- function(x) {
  cns <- x %>%
    igraph::simplify(edge.attr.comb = "random") # To make sure there is no duplicate or multiplex edges
  cns <- igraph::as_data_frame(cns, what = "edges") %>%
    rename(Node.1 = from, Node.2 = to)

  # calculate the weight for each consensus edge
  data_x <- cns[, -c(1, 2)]
  sum_x <- apply(data_x, 1, sum) # calculate the sum of each row (all weights of that edge)
  div <- (data_x / sum_x) * data_x # divide each weight (in a row) by the sum of the row and multiply by the original value
  wTO_cons <- apply(div, 1, sum) # Sum of the div values per row is assigned as weight of the that consensus edge

  # discard edges with mismatch signs
  cons_wto <- cns %>%
    mutate(CONS = wTO_cons) %>%
    mutate(CN = case_when(
      if_all(starts_with("wTO_"), ~ . > 0) ~ CONS,
      if_all(starts_with("wTO_"), ~ . < 0) ~ CONS,
      TRUE ~ 0
    )) %>%
    mutate(CN = round(CN, digits = 2)) %>%
    dplyr::select(Node.1, Node.2, CN) %>%
    filter(abs(CN) > 0)

  return(cons_wto)
}
