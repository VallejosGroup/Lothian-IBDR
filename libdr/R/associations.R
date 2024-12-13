#' Percentage bar plots for categorical data
#' @param dat Data frame holding categorical and cluster assignment data.
#' @param var Character. The name of the categorical variable of interest.
#' @param class Character. The name of the class assignment variable. Assumed
#'   to be \code{"class_combined"} if not manually specified
#' @returns A \code{\link[patchwork]{patchwork}} object
#' @export
plotCat <- function(dat, var, class = "class_combined") {
  cluster <- NULL

  fill.vec <- c("#C2F8CB", "#EFC7E5", "#F6BD60", "blue")
  col.vec <- c("#82C68F", "#C795BB", "#C89216", "blue")

  if (!(var %in% colnames(dat))) stop("Column for categorical variable not found.")
  if (!(class %in% colnames(dat))) stop("Class column not found.")

  labels <- sort(unique(dat[, var]))

  # Create df to hold percentages
  perc.table <- data.frame(
    label = character(),
    cluster = character(),
    perc = numeric()
  )

  # Calculate percentage per cluster for each category
  for (g in sort(unique(dat[, class]))) {
    temp.1 <- dat[dat[, class] == g, ]
    for (label in labels) {
      perc <- nrow(subset(temp.1, eval(parse(text = var)) == label)) / nrow(temp.1)
      perc.table <- rbind(
        perc.table,
        data.frame(
          label = label,
          cluster = g,
          perc = perc
        )
      )
    }
  }

  # Set order of clusters in percentage data frame
  perc.table$cluster <- factor(perc.table$cluster,
    levels = levels(dat[, class])
  )
  # Per sub plot
  for (i in 1:length(labels)) {
    totalPerc <- nrow(subset(dat, eval(parse(text = var)) == labels[i])) /
      nrow(dat)

    if (i == 1) {
      p <- perc.table %>%
        filter(label == labels[i]) %>%
        ggplot(aes(x = cluster, y = perc)) +
        geom_bar(stat = "identity", fill = fill.vec[i], color = col.vec[i]) +
        geom_hline(yintercept = totalPerc, linetype = "dashed", color = "#4D4730") +
        theme_minimal() +
        labs(x = "Cluster", y = "Percentage") +
        ggtitle(labels[i])
    } else {
      p <- p +
        perc.table %>%
        filter(label == labels[i]) %>%
        ggplot(aes(x = cluster, y = perc)) +
        geom_bar(stat = "identity", fill = fill.vec[i], color = col.vec[i]) +
        geom_hline(yintercept = totalPerc, linetype = "dashed", color = "#4D4730") +
        theme_minimal() +
        labs(x = "Cluster", y = "Percentage") +
        ggtitle(labels[i])
    }
  }
  p <- p +
    patchwork::plot_layout(nrow = length(labels), ncol = 1, guides = "collect") &
    # patchwork::plot_annotation(tag_levels = "A") &
    scale_y_continuous(labels = scales::percent, limits = c(0, 1))
  return(p)
}

.getCI <- function(m) {
  Var1 <- Var2 <- var <- SE <- Estimate <- z <- Upper <- Lower <- NULL

  CI_factor <- qnorm(0.975)

  # A chatGPT edited version of the original code :-)
  m_summary <- summary(m)
  m_coeff <- as.data.frame(m_summary$coefficients) %>%
    tibble::rownames_to_column("Var1") %>%
    tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "Estimate") %>%
    dplyr::mutate(var = paste(Var1, Var2))
  m_se <- as.data.frame(m_summary$standard.errors) %>%
    tibble::rownames_to_column("Var1") %>%
    tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "SE") %>%
    dplyr::mutate(var = paste(Var1, Var2)) %>%
    dplyr::select(c(var, SE))
  m_coeff <- merge(m_coeff, m_se, by = "var") %>%
    dplyr::mutate(
      z = Estimate / SE,
      p.val = (1 - pnorm(abs(z), 0, 1)) * 2,
      Lower = Estimate - CI_factor * SE,
      Upper = Estimate + CI_factor * SE,
      Sig = ifelse(Upper < 0 | Lower > 0, TRUE, FALSE)
    ) %>%
    dplyr::arrange(Var2, Var1) %>%
    subset(Var2 != "(Intercept)")

  return(m_coeff)
}

.plotCI <- function(tab, variable) {
  Var2 <- Var1 <- Estimate <- Lower <- Upper <- Sig <- Model <- NULL

  if (is.null(tab$Model)) {
    p <- tab %>%
      dplyr::filter(Var2 == variable) %>%
      ggplot(aes(
        x = Var1,
        y = Estimate,
        ymin = Lower,
        ymax = Upper,
        color = ifelse(Sig == TRUE, "red", "black")
      )) +
      geom_errorbar() +
      geom_point(size = 3.5) +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() + # flip coordinates (puts labels on y axis)
      xlab("") +
      ylab("Estimate (95% CI)") +
      theme_minimal() +
      scale_color_manual(values = c("black", "#FF007F")) +
      theme(legend.position = "none") +
      ggtitle(variable)
  } else {
    # Filter for the specified variable and plot
    p <- tab %>%
      filter(Var2 == variable) %>%
      ggplot(aes(
        x = Var1,
        y = Estimate,
        ymin = Lower,
        ymax = Upper,
        color = Model
      )) + # Use model as color
      geom_errorbar(position = position_dodge(width = 0.7)) + # Dodge for separation
      geom_point(position = position_dodge(width = 0.7), size = 2.5) +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() + # Flip coordinates (puts labels on y axis)
      xlab("") +
      ylab("Estimate (95% CI)") +
      theme_minimal() +
      ggtitle(variable)
  }

  return(p)
}

#' Forest plot of multinomial logistic regression model
#' @param dat  Data frame holding covariate and cluster assignment data.
#' @param var Character. The name of the covariate(s) of interest. If given as
#'   a character vector then a multivariate model is fitted using the specified
#'   covariates.
#' @param class Character. The name of the class assignment variable. Assumed
#'   to be \code{"class_combined"} if not manually specified.
#' @param prob Character. The name of the variable which gives posterior
#'   probabilities for cluster membership. Assumed to be \code{"probmax"} if not
#'   manually specified.
#' @param minprob Numeric. The minimum posterior probability for cluster
#'   membership required for a subject to be included in the analysis. Assumed
#'   to be 0.5 if not manually specified
#' @param extern Optional data frame containing multinomial logistic
#'   regression results obtained from an external dataset. If not given then it
#'   is assumed that no external data is required.
#' @export
mlrPlot <- function(dat, var, class = "class_combined", prob = "probmax", minprob = 0.5, extern = NULL) {

  #print("Multivariate analysis")
  # Multivariate analysis
  ## Fit the model
  mlr_multi <- nnet::multinom(formula = reformulate(var, class), data = dat, trace = FALSE)
  ## Extract coeffs and CIs
  tab_multi <- .getCI(mlr_multi)
  ## Plot
  p_multi <- list()
  for (variable in unique(tab_multi$Var2)) {
    p_multi[[variable]] <- .plotCI(tab_multi, variable)
  }

  #print("Multivariate analysis - excluding those with prob < minprob")
  # Multivariate analysis - excluding those with prob < minprob
  ## Fit the model
  dat.minprob <- dat %>% filter(!!sym(prob) >= minprob)
  mlr_multi_minprob <- nnet::multinom(
    formula = reformulate(var, class),
    data = dat.minprob, trace = FALSE
  )
  ## Extract coeffs and CIs
  tab_multi_minprob <- .getCI(mlr_multi_minprob)
  ## Plot
  p_multi_minprob <- list()
  for (variable in unique(tab_multi$Var2)) {
    p_multi_minprob[[variable]] <- .plotCI(tab_multi_minprob, variable)
  }

  if (length(var) > 1) {
    #print("Everyone")
    ### Everyone
    ## Fit the model + Extract coeffs and CIs

    if (is.null(extern)) {
    # Univariate analysis

      tab_uni <- NULL
      for (variable in var) {
        mlr_uni <- nnet::multinom(
          formula = reformulate(variable, class),
          data = dat, trace = FALSE
        )
        tab_uni <- rbind(tab_uni, .getCI(mlr_uni))
      }
      ## Plot
      p_uni <- list()
      for (variable in unique(tab_multi$Var2)) {
        p_uni[[variable]] <- .plotCI(tab_uni, variable)
      }


      tab_uni$Model <- "Univariate"
      tab_multi$Model <- "Multivariate"
      tab_both <- bind_rows(tab_uni, tab_multi)
      tab_both$Model <- as.factor(tab_both$Model)
      p_both <- list()
      for (variable in unique(tab_multi$Var2)) {
        p_both[[variable]] <- .plotCI(tab = tab_both, variable)
      }

      #print("Excluding those with prob < minprob")
      ### Excluding those with prob < minprob
      ## Fit the model + Extract coeffs and CIs
      tab_uni_minprob <- NULL
      for (variable in var) {
        mlr_uni_minprob <- nnet::multinom(
          formula = reformulate(variable, class),
          data = dat.minprob, trace = FALSE
        )
        tab_uni_minprob <- rbind(tab_uni_minprob, .getCI(mlr_uni_minprob))
      }
      ## Plot
      p_uni_minprob <- list()
      for (variable in unique(tab_multi$Var2)) {
        p_uni_minprob[[variable]] <- .plotCI(tab_uni_minprob, variable)
      }
      ## Combined plot
      # Add a model identifier to each dataset and combine the two tables
      tab_uni_minprob$Model <- paste0("Univariate (prob >=", minprob, ")")
      tab_multi_minprob$Model <- paste0("Multivariate (prob >=", minprob, ")")
      tab_both_minprob <- bind_rows(tab_uni_minprob, tab_multi_minprob)
      tab_both_minprob$Model <- as.factor(tab_both_minprob$Model)
      p_both_minprob <- list()
      for (variable in unique(tab_uni$Var2)) {
        p_both_minprob[[variable]] <- .plotCI(tab = tab_both_minprob, variable)
      }

      ## Full-combined plot
      tab_everything <- bind_rows(
        tab_uni, tab_multi,
        tab_uni_minprob, tab_multi_minprob
      )
      tab_everything$Model <- as.factor(tab_everything$Model)
      p_everything <- list()
      for (variable in unique(tab_multi$Var2)) {
        p_everything[[variable]] <- .plotCI(tab = tab_everything, variable)
      }
    # Output
    if (length(var) == 1) {
      out <- list(
        tab_multi = tab_multi,
        plot_multi = p_multi,
        tab_multi_minprob = tab_multi_minprob,
        plot_multi_minprob = p_multi_minprob
      )
    } else {
      out <- list(
        tab_multi = tab_multi,
        plot_multi = p_multi,
        tab_uni = tab_uni,
        plot_uni = p_uni,
        plot_both = p_both,
        tab_everything = tab_everything,
        plot_everything = p_everything
      )
    }
    return(out)
  } else {

    tab_multi$Model <- "LIBDR Multivariate"
    tab_both <- bind_rows(tab_multi, extern)
    tab_both$Model <- as.factor(tab_both$Model)
    p_both <- list()
    for (variable in unique(tab_multi$Var2)) {
      p_both[[variable]] <- .plotCI(tab = tab_both, variable)
    }

    tab_multi_minprob$Model <- paste0("LIBDR Multivariate (prob >=", minprob, ")")

    extern <- subset(extern, Model %in% c("Multivariate", "Multivariate (prob >=0.5)"))
    extern$Model <- paste("Danish", extern$Model)


    ## Full-combined plot
    tab_everything <- bind_rows(
      tab_multi,
      tab_multi_minprob,
      extern
    )
    tab_everything$Model <- factor(tab_everything$Model,
                                   levels = c("Danish Multivariate (prob >=0.5)",
                                              "Danish Multivariate",
                                              "LIBDR Multivariate (prob >=0.5)",
                                              "LIBDR Multivariate"))

    p_everything <- list()
    for (variable in unique(tab_multi$Var2)) {
      p_everything[[variable]] <- .plotCI(tab = tab_everything, variable) +
        guides(color = guide_legend(reverse = TRUE)) +
        scale_color_manual(values = c("#7D0A1D", "#DC1435", "#00488C", "#0071DA"))
    }

  # Output
  if (length(var) == 1) {
    out <- list(
      tab_multi = tab_multi,
      plot_multi = p_multi,
      tab_multi_minprob = tab_multi_minprob,
      plot_multi_minprob = p_multi_minprob
    )
  } else {
    out <- list(
      tab_multi = tab_multi,
      plot_multi = p_multi,
      plot_both = p_both,
      tab_everything = tab_everything,
      plot_everything = p_everything
    )
  }
  return(out)
  }
  }
}

