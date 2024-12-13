# Legacy code deleted from the associations.qmd file

##### Montreal classification of age for Crohn's disease

Unlike ulcerative colitis, age at diagnosis is part of the Montreal
classification for Crohn's disease.  As such, we will explore potential
associations between Montreal age classification and cluster membership for
subjects diagnosed with Crohn's disease. The Montreal age categories are as
follows:

  * A1: $\leq$ 16 years
* A2: 17—40 years
* A3: $>$ 40 years

Montreal classification of age for Crohn's disease was not found to be
significantly associated with cluster membership.

```{R}
ageLabels <- c("\u226416", "17—40", "> 40")

myDF.fc.cd <- subset(myDF.fc, diagnosis == "Crohn's Disease")

myDF.fc.cd$ageCat <- cut(myDF.fc.cd$age,
  breaks = c(0, 16, 40, Inf),
  labels = ageLabels,
  include.lowest = TRUE,
  right = FALSE
)

myDF.fc.cd %>%
  filter(probmax_combined > 0.5) %>%
  with(fisher.test(class_combined,
                   ageCat,
                   workspace = 2000000,
                   simulate.p.value = TRUE)) %>%
  pander()
```

```{R}
myDF.fc.cd %>%
  filter(probmax_combined > 0.5) %>%
  mlrPlot(var = "ageCat")
```

#### Time to first line AT usage

```{r}
myDF.fc <- myDF.fc %>%
  mutate(AT_line_1_censor = if_else(AT == 0, 7, as.numeric(AT_line_1))) %>%
  mutate(AT_line_1_censor = if_else(AT_line_1_censor > 7, 7, AT_line_1_censor),
         AT_censor = if_else(AT_line_1_censor == 7, 0, AT))

p_tte <- survfit(Surv(AT_line_1_censor, AT_censor) ~ class_order, data = myDF.fc) %>%
  ggsurvplot(pval = TRUE) +
  ggtitle("CD + UC + IBDU")

ATcum<- cmprsk::cuminc(
  ftime = myDF.fc$AT_line_1_censor,
  fstatus = myDF.fc$AT_censor,
  group = myDF.fc$class_order)

ATcum_cd <- cmprsk::cuminc(
  ftime = myDF.fc$AT_line_1_censor,
  fstatus = myDF.fc$AT_censor,
  group = myDF.fc$class_order,
  subset = ifelse(myDF.fc$diagnosis == "Crohn's Disease", TRUE, FALSE))

ATcum_uc <- cmprsk::cuminc(
  ftime = myDF.fc$AT_line_1_censor,
  fstatus = myDF.fc$AT_censor,
  group = myDF.fc$class_order,
  subset = ifelse(myDF.fc$diagnosis == "Ulcerative Colitis", TRUE, FALSE))

  ggcompetingrisks(fit) + facet_wrap(ncol = 2)
  ggtitle("CD + UC + IBDU")

  ggcompetingrisks(ATcum_uc)
```


#### All lines

```{R}
fc.cd.at <-
fc.cd.at$Line <- as.factor(str_replace(fc.cd.at$Line, "AT_line_", ""))

myDF.fc.cd <- subset(myDF.fc, diagnosis == "Crohn's Disease")



scaling <- table(fc.cd.at$class)/table(myDF.fc.cd$class)


plots <- list()




for (g in 1:8) {
  prediction <- pred[, c(g, g + 8, g + 16, 25)]
  names(prediction) <- c("mean", "lower", "upper", "time")

  scaler <- scaling[[g]]

  temp <- fc.cd.at %>%
    filter(class == paste0(g)) %>%
    filter(value < 7)

  dens <- density(temp$value, kernel = "gaussian", from = 0, to = 7)
  dens$y <- dens$y / max(dens$y) # Scale density to 1
  dens <- data.frame(x = dens$x, y = dens$y)
  dens$y <- dens$y * scaler * 4

  plots[[g]] <- dens %>%
    ggplot() +
    geom_area(aes(x = x,
                  y = y),
                 fill = "#AB87FF",
color =  "#7903F0",
linewidth = 1.2,
lty = 3,
alpha = 0.7) +
  geom_line(aes(x = time, y = mean),
            prediction,
            color = "red",
            linewidth = 1.2) +
  geom_hline(yintercept = log(250), lty = 2, color =  "#007add") +
  scale_y_continuous(
    name = "Log (FC (\u03BCg/g))",
    limits = c(0, ylimit),
    sec.axis = sec_axis(transform = ~./ylimit,
                        name = "Advanced therapy density")) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))  +
  xlab("Time (years)") +
  ggtitle(paste0("FC", title.mapping[g]))
}

mapping <- c(5, 7, 8, 3, 6, 2, 1, 4)

for (i in 1:8) {
  if (i == 1) {
    p <- plots[[mapping[i]]]
  } else {
    p <- p + plots[[mapping[i]]]
  }
}

p <- p + plot_layout(ncol = 2)
ggsave("plots/fc-cd-at.png", p, width = 10, height = 14, units = "in")
p
```


```{R}
fc.cd.at <- reshape2::melt(dict.fc.cd,
                           id.vars = "ids",
                           measure.vars = paste0("AT_line_", seq(1, 8)),
                           variable.name = "Line") %>%
  drop_na(value)
fc.cd.at$Line <- as.factor(str_replace(fc.cd.at$Line, "AT_line_", ""))

myDF.fc.cd <- subset(myDF.fc, diagnosis == "Crohn's Disease")


fc.cd.at <- merge(fc.cd.at,
                  myDF.fc[, c("ids", "class")],
                  by = "ids",
                  all.x = TRUE,
                  all.y = FALSE)


time.pred <- seq(0, 7, by = 0.01)

pred.fc.df <- data.frame(
  calpro_time = c(time.pred, time.pred),
  diagnosis = c(
    rep("Crohn's Disease", length(time.pred)),
    rep("Ulcerative Colitis", length(time.pred))
  )
)
pred.fc.df.update <- lcmm::predictY(model.fc,
                                    pred.fc.df,
                                    var.time = "calpro_time",
                                    draws = TRUE
)$pred

pred <- predictY(model.fc, pred.fc.df, var.time = "calpro_time", draws = TRUE)$pred

pred <- as.data.frame(pred[1:length(time.pred), ])
pred$time <- time.pred

plots <- list()

ylimit <- log(2500)
title.mapping <- c(7, 6, 4, 8, 1, 5, 2, 3)


for (g in 1:8) {
  prediction <- pred[, c(g, g + 8, g + 16, 25)]
  names(prediction) <- c("mean", "lower", "upper", "time")


  temp <- fc.cd.at %>%
    filter(class == paste0(g)) %>%
    filter(value < 7)


  temp %>%
    summarise()

  plots[[g]] <- temp %>%
    ggplot() +
    geom_histogram(aes(x = value, y = after_stat(count)),
                   fill = "#AB87FF",
                   color =  "#7903F0") +
    geom_line(aes(x = time, y = mean),
              prediction,
              color = "red",
              linewidth = 1.2) +
    geom_hline(yintercept = log(250), lty = 2, color =  "#007add") +
    scale_y_continuous(
      name = "Log (FC (\u03BCg/g))",
      limits = c(0, ylimit),
      sec.axis = sec_axis(transform = ~./ylimit,
                          name = "Advanced therapy density")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))  +
    xlab("Time (years)") +
    ggtitle(paste0("FC", title.mapping[g]))
}

mapping <- c(5, 7, 8, 3, 6, 2, 1, 4)

for (i in 1:8) {
  if (i == 1) {
    p <- plots[[mapping[i]]]
  } else {
    p <- p + plots[[mapping[i]]]
  }
}

p <- p + plot_layout(ncol = 2)
ggsave("plots/fc-cd-at-hist.png", p, width = 10, height = 14, units = "in")
p
```

#### First-line only

```{R}
myDF.fc.cd <- subset(myDF.fc, diagnosis == "Crohn's Disease")
fc.cd.at <- fc.cd.at %>%
  filter(Line == 1)

scaling <- table(fc.cd.at$class)/table(myDF.fc.cd$class)

time.pred <- seq(0, 7, by = 0.01)

pred.fc.df <- data.frame(
  calpro_time = c(time.pred, time.pred),
  diagnosis = c(
    rep("Crohn's Disease", length(time.pred)),
    rep("Ulcerative Colitis", length(time.pred))
  )
)
pred.fc.df.update <- lcmm::predictY(model.fc,
                                    pred.fc.df,
                                    var.time = "calpro_time",
                                    draws = TRUE
)$pred

pred <- predictY(model.fc, pred.fc.df, var.time = "calpro_time", draws = TRUE)$pred

pred <- as.data.frame(pred[1:length(time.pred), ])
pred$time <- time.pred

plots <- list()

ylimit <- log(2500)
title.mapping <- c(7, 6, 4, 8, 1, 5, 2, 3)


for (g in 1:8) {
  prediction <- pred[, c(g, g + 8, g + 16, 25)]
  names(prediction) <- c("mean", "lower", "upper", "time")

  scaler <- scaling[[g]]

  temp <- fc.cd.at %>%
    filter(class == paste0(g)) %>%
    filter(value < 7)

  dens <- density(temp$value, kernel = "gaussian", from = 0, to = 7)
  dens$y <- dens$y / max(dens$y) # Scale density to 1
  dens <- data.frame(x = dens$x, y = dens$y)
  dens$y <- dens$y * scaler * 4

  plots[[g]] <- dens %>%
    ggplot() +
    geom_area(aes(x = x,
                  y = y),
              fill = "#AB87FF",
              color =  "#7903F0",
              linewidth = 1.2,
              lty = 3,
              alpha = 0.7) +
    geom_line(aes(x = time, y = mean),
              prediction,
              color = "red",
              linewidth = 1.2) +
    geom_hline(yintercept = log(250), lty = 2, color =  "#007add") +
    scale_y_continuous(
      name = "Log (FC (\u03BCg/g))",
      limits = c(0, ylimit),
      sec.axis = sec_axis(transform = ~./ylimit,
                          name = "Advanced therapy density")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))  +
    xlab("Time (years)") +
    ggtitle(paste0("FC", title.mapping[g]))
}

mapping <- c(5, 7, 8, 3, 6, 2, 1, 4)

for (i in 1:8) {
  if (i == 1) {
    p <- plots[[mapping[i]]]
  } else {
    p <- p + plots[[mapping[i]]]
  }
}

p <- p + plot_layout(ncol = 2)
ggsave("plots/fc-cd-at-firstline.png", p, width = 10, height = 14, units = "in")
p
```


```{R}
dict.fc.cd <- subset(myDF.fc, diagnosis == "Crohn's Disease" & ids %in% fcal$ids)

#Apply censoring
dict.fc.cd <- dict.fc.cd %>%
  mutate(AT_line_1  = if_else(AT == 0, 7, as.numeric(AT_line_1))) %>%
  mutate(AT = if_else(AT_line_1 > 7, 0, 1)) %>%
  mutate(AT_line_1 = if_else(AT_line_1  > 7, 7, AT_line_1),
         AT = if_else(AT_line_1 == 7, 0, AT))
survfit(Surv(AT_line_1, AT) ~ class_order, data = dict.fc.cd) %>%
  ggsurvplot(pval = TRUE)
```

```{R}
#| echo: false
png("plots/at-survival-fc-ci.png",
    width = 7,
    height = 5,
    units = "in",
    res = 300
)
survfit(Surv(AT_line_1, AT) ~ class_order, data = dict.fc.cd) %>%
  ggsurvplot(pval = TRUE, conf.int = TRUE)
invisible(dev.off())
```

### Advanced therapy in UC

#### All lines

```{R}
dict.fc.uc <- dict %>%
  filter(ids %in% fcal$ids) %>%
  filter(diagnosis == "Ulcerative Colitis")
```

```{R}
fc.uc.at <- reshape2::melt(dict.fc.uc,
                           id.vars = "ids",
                           measure.vars = paste0("AT_line_", seq(1, 8)),
                           variable.name = "Line") %>%
  drop_na(value) %>%
  mutate(value = as.numeric(value))
fc.uc.at$Line <- as.factor(str_replace(fc.uc.at$Line, "AT_line_", ""))

myDF.fc.uc <- subset(myDF.fc, diagnosis == "Ulcerative Colitis")


fc.uc.at <- merge(fc.uc.at,
                  myDF.fc[, c("ids", "class")],
                  by = "ids",
                  all.x = TRUE,
                  all.y = FALSE)

scaling <- table(fc.uc.at$class)/table(myDF.fc.uc$class)

time.pred <- seq(0, 7, by = 0.01)

pred.fc.df <- data.frame(
  calpro_time = c(time.pred, time.pred),
  diagnosis = c(
    rep("Crohn's Disease", length(time.pred)),
    rep("Ulcerative Colitis", length(time.pred))
  )
)
pred.fc.df.update <- lcmm::predictY(model.fc,
                                    pred.fc.df,
                                    var.time = "calpro_time",
                                    draws = TRUE
)$pred

pred <- predictY(model.fc, pred.fc.df, var.time = "calpro_time", draws = TRUE)$pred

pred <- as.data.frame(pred[1:length(time.pred), ])
pred$time <- time.pred

plots <- list()

ylimit <- log(2500)
title.mapping <- c(7, 6, 4, 8, 1, 5, 2, 3)


for (g in 1:8) {
  prediction <- pred[, c(g, g + 8, g + 16, 25)]
  names(prediction) <- c("mean", "lower", "upper", "time")

  scaler <- scaling[[g]]

  temp <- fc.uc.at %>%
    filter(class == paste0(g)) %>%
    filter(value < 7)

  dens <- density(temp$value, kernel = "gaussian", from = 0, to = 7)
  dens$y <- dens$y / max(dens$y) # Scale density to 1
  dens <- data.frame(x = dens$x, y = dens$y)
  dens$y <- dens$y * scaler * 4

  plots[[g]] <- dens %>%
    ggplot() +
    geom_area(aes(x = x,
                  y = y),
              fill = "#AB87FF",
              color =  "#7903F0",
              linewidth = 1.2,
              lty = 3,
              alpha = 0.7) +
    geom_line(aes(x = time, y = mean),
              prediction,
              color = "red",
              linewidth = 1.2) +
    geom_hline(yintercept = log(250), lty = 2, color =  "#007add") +
    scale_y_continuous(
      name = "Log (FC (\u03BCg/g))",
      limits = c(0, ylimit),
      sec.axis = sec_axis(transform = ~./ylimit,
                          name = "Advanced therapy density")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))  +
    xlab("Time (years)") +
    ggtitle(paste0("FC", title.mapping[g]))
}


mapping <- c(5, 7, 8, 3, 6, 2, 1, 4)

for (i in 1:8) {
  if (i == 1) {
    p <- plots[[mapping[i]]]
  } else {
    p <- p + plots[[mapping[i]]]
  }
}

p <- p + plot_layout(ncol = 2)
ggsave("plots/fc-uc-at.png", p, width = 10, height = 14, units = "in")
p
```


#### First-line only

```{R}
myDF.fc.uc <- subset(myDF.fc, diagnosis == "Ulcerative Colitis")
fc.uc.at <- fc.uc.at %>%
  filter(Line == 1)

scaling <- table(fc.uc.at$class)/table(myDF.fc.uc$class)

time.pred <- seq(0, 7, by = 0.01)

pred.fc.df <- data.frame(
  calpro_time = c(time.pred, time.pred),
  diagnosis = c(
    rep("Crohn's Disease", length(time.pred)),
    rep("Ulcerative Colitis", length(time.pred))
  )
)
pred.fc.df.update <- lcmm::predictY(model.fc,
                                    pred.fc.df,
                                    var.time = "calpro_time",
                                    draws = TRUE
)$pred

pred <- predictY(model.fc, pred.fc.df, var.time = "calpro_time", draws = TRUE)$pred

pred <- as.data.frame(pred[1:length(time.pred), ])
pred$time <- time.pred

plots <- list()

ylimit <- log(2500)
title.mapping <- c(7, 6, 4, 8, 1, 5, 2, 3)


for (g in 1:8) {
  prediction <- pred[, c(g, g + 8, g + 16, 25)]
  names(prediction) <- c("mean", "lower", "upper", "time")

  scaler <- scaling[[g]]

  temp <- fc.uc.at %>%
    filter(class == paste0(g)) %>%
    filter(value < 7)

  dens <- density(temp$value, kernel = "gaussian", from = 0, to = 7)
  dens$y <- dens$y / max(dens$y) # Scale density to 1
  dens <- data.frame(x = dens$x, y = dens$y)
  dens$y <- dens$y * scaler * 4

  plots[[g]] <- dens %>%
    ggplot() +
    geom_area(aes(x = x,
                  y = y),
              fill = "#AB87FF",
              color =  "#7903F0",
              linewidth = 1.2,
              lty = 3,
              alpha = 0.7) +
    geom_line(aes(x = time, y = mean),
              prediction,
              color = "red",
              linewidth = 1.2) +
    geom_hline(yintercept = log(250), lty = 2, color =  "#007add") +
    scale_y_continuous(
      name = "Log (FC (\u03BCg/g))",
      limits = c(0, ylimit),
      sec.axis = sec_axis(transform = ~./ylimit,
                          name = "Advanced therapy density")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))  +
    xlab("Time (years)") +
    ggtitle(paste0("FC", title.mapping[g]))
}


mapping <- c(5, 7, 8, 3, 6, 2, 1, 4)

for (i in 1:8) {
  if (i == 1) {
    p <- plots[[mapping[i]]]
  } else {
    p <- p + plots[[mapping[i]]]
  }
}

p <- p + plot_layout(ncol = 2)
ggsave("plots/fc-uc-at-firstline.png", p, width = 10, height = 14, units = "in")
p
```

```{R}
dict.fc.uc <- subset(myDF.fc, diagnosis == "Ulcerative Colitis" & ids %in% fcal$ids)

#Apply censoring
dict.fc.uc <- dict.fc.uc %>%
  mutate(AT_line_1  = if_else(AT == 0, 7, as.numeric(AT_line_1 ))) %>%
  mutate(AT = if_else(AT_line_1 > 7, 0, 1)) %>%
  mutate(AT_line_1 = if_else(AT_line_1  > 7, 7, AT_line_1),
         AT = if_else(AT_line_1 == 7, 0, AT))
survfit(Surv(AT_line_1, AT) ~ class_order, data = dict.fc.uc) %>%
  ggsurvplot(pval = TRUE)
```

### Advanced therapy for CD and UC

#### All lines

```{R}

```

```{R}
dict.fc <- dict %>%
  filter(ids %in% fcal$ids) %>%
  filter(diagnosis == "Crohn's Disease" | diagnosis == "Ulcerative Colitis")
```

```{R}
fc.at <- reshape2::melt(dict.fc,
                        id.vars = "ids",
                        measure.vars = paste0("AT_line_", seq(1, 8)),
                        variable.name = "Line") %>%
  drop_na(value) %>%
  mutate(value = as.numeric(value))
fc.at$Line <- as.factor(str_replace(fc.at$Line, "AT_line_", ""))

myDF.fc.temp <- subset(myDF.fc, diagnosis == "Crohn's Disease" |
                         diagnosis == "Ulcerative Colitis")


fc.at <- merge(fc.at,
               myDF.fc[, c("ids", "class")],
               by = "ids",
               all.x = TRUE,
               all.y = FALSE)

scaling <- table(fc.at$class)/table(myDF.fc.temp$class)

time.pred <- seq(0, 7, by = 0.01)

pred.fc.df <- data.frame(
  calpro_time = c(time.pred, time.pred),
  diagnosis = c(
    rep("Crohn's Disease", length(time.pred)),
    rep("Ulcerative Colitis", length(time.pred))
  )
)
pred.fc.df.update <- lcmm::predictY(model.fc,
                                    pred.fc.df,
                                    var.time = "calpro_time",
                                    draws = TRUE
)$pred

pred <- predictY(model.fc, pred.fc.df, var.time = "calpro_time", draws = TRUE)$pred

pred <- as.data.frame(pred[1:length(time.pred), ])
pred$time <- time.pred

plots <- list()

ylimit <- log(2500)
title.mapping <- c(7, 6, 4, 8, 1, 5, 2, 3)


for (g in 1:8) {
  prediction <- pred[, c(g, g + 8, g + 16, 25)]
  names(prediction) <- c("mean", "lower", "upper", "time")

  scaler <- scaling[[g]]

  temp <- fc.at %>%
    filter(class == paste0(g)) %>%
    filter(value < 7)

  dens <- density(temp$value, kernel = "gaussian", from = 0, to = 7)
  dens$y <- dens$y / max(dens$y) # Scale density to 1
  dens <- data.frame(x = dens$x, y = dens$y)
  dens$y <- dens$y * scaler * 4

  plots[[g]] <- dens %>%
    ggplot() +
    geom_area(aes(x = x,
                  y = y),
              fill = "#AB87FF",
              color =  "#7903F0",
              linewidth = 1.2,
              lty = 3,
              alpha = 0.7) +
    geom_line(aes(x = time, y = mean),
              prediction,
              color = "red",
              linewidth = 1.2) +
    geom_hline(yintercept = log(250), lty = 2, color =  "#007add") +
    scale_y_continuous(
      name = "Log (FC (\u03BCg/g))",
      limits = c(0, ylimit),
      sec.axis = sec_axis(transform = ~./ylimit,
                          name = "Advanced therapy density")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))  +
    xlab("Time (years)") +
    ggtitle(paste0("FC", title.mapping[g], "CD: "))
}


mapping <- c(5, 7, 8, 3, 6, 2, 1, 4)

for (i in 1:8) {
  if (i == 1) {
    p <- plots[[mapping[i]]]
  } else {
    p <- p + plots[[mapping[i]]]
  }
}

p <- p + plot_layout(ncol = 2)
ggsave("plots/fc-at.png", p, width = 10, height = 14, units = "in")
p
```






```{R}

dat <- data.frame(Cluster = character(), CD = numeric(), UC = numeric())


for (g in 1:8) {
  total <- myDF.fc %>%
    mutate(AT_line_1 = if_else(AT == 0, 7, as.numeric(AT_line_1))) %>%
    mutate(AT = if_else(AT_line_1 >= 7, 0, AT)) %>%
    mutate(AT_line_1  = if_else(AT_line_1  > 7, 7, AT_line_1)) %>%
    mutate(AT_line_1 = if_else(AT_line_1 < 0, 0, AT_line_1)) %>%
    filter(class_order == paste0("FC", g)) %>%
    summarise(CD = sum(diagnosis == "Crohn's Disease"),
              UC = sum(diagnosis == "Ulcerative Colitis"))

  therapy <- myDF.fc %>%
    mutate(AT_line_1 = if_else(AT == 0, 7, as.numeric(AT_line_1))) %>%
    mutate(AT = if_else(AT_line_1 >= 7, 0, AT)) %>%
    mutate(AT_line_1  = if_else(AT_line_1  > 7, 7, AT_line_1)) %>%
    mutate(AT_line_1 = if_else(AT_line_1 < 0, 0, AT_line_1)) %>%
    filter(class_order == paste0("FC", g)) %>%
    filter(AT == 1) %>%
    summarise(CD = sum(diagnosis == "Crohn's Disease"),
              UC = sum(diagnosis == "Ulcerative Colitis"))

  dat <- rbind(dat,
               data.frame(Cluster = paste0("FC", g),
                          CD = (therapy$CD/total$CD) * 100,
                          UC = (therapy$UC/total$UC) * 100))
}
```

## CRP

```{R}
model.crp$pprob$class <- plyr::mapvalues(
  model.crp$pprob$class,
  from = c(8, 9),
  to = c(7, 8)
)
```

As Cluster CRP2 and CRP4 are characterised by very similar longitudinal
profiles, we have elected to merge these clusters. All subjects in cluster
4 have been migrated to cluster 2 with clusters 5-9 renumbered accordingly.

Similar to FC, we will explore associations only for subjects with a
sufficiently high posterior probability of being in their assigned cluster. As
cluster assignment is exclusive, we are able to add the posterior probabilities
for the merged clusters together to determine the posterior probability for
CRP2 after the merge.

### Sex

```{R}
dict.crp <- subset(dict, ids %in% model.crp$pprob$ids)
dict.crp <- merge(dict.crp, model.crp$pprob[, c("ids", "class")])


crp.median <-
  crp.median %>%
  select(-sex) %>%
  merge(dict.crp[, c("ids", "class", "sex")],
        all.x = TRUE,
        all.y = FALSE
  )

chisq.test(dict.crp$class, dict.crp$sex) %>%
  pander()
```

### Age

Age has been calculated by subtracting year of diagnosis from year of birth.

```{R}
dict.crp$ageCat <- cut(dict.crp$age,
                       breaks = c(0, 18, 30, 50, 70, Inf),
                       labels = c("<18", "18-29", "30-49", "50-69", "\u226570"),
                       include.lowest = TRUE,
                       right = FALSE
)

crp.median <- merge(crp.median,
                    dict.crp[, c("ids", "ageCat")],
                    by = "ids",
                    all.x = TRUE,
                    all.y = FALSE
)

aov(class ~ age, data = dict.crp) %>%
  summary() %>%
  pander()
```


```{R}
#| label: fig-crp-age-dist
#| fig-cap: "Distribution of age at diagnosis by CRP cluster assignment."
#| fig-height: 4
#| fig-width: 9.142857
#| warning: false
p2 <- dict.crp %>%
  mutate(class = factor(class, labels = paste("CRP", seq(1, 8)))) %>%
  ggplot(aes(x = class, y = age)) +
  geom_violin(fill = "#D4CDF4", color = "#7A6F9B") +
  theme_minimal() +
  xlab("Cluster") +
  ylab("Age at diagnosis")

p <- p1 / p2 + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave("plots/age-dist.png", p, width = 16 * 4 / 7, height = 8)
ggsave("plots/age-dist.pdf", p, width = 16 * 4 / 7, height = 8)
print(p2)
```


### IBD type

As with the FC  model, IBD type has been included in the cluster assignment sub
model as a covariate. Crohn's disease is the reference IBD type with cluster
CRP 8 being the reference cluster.


```{R}
#| label: fig-crp-ibd
#| fig-cap: "Forest plot of IBD type for the chosen CRP model."
n.clust <- 8
id <- 1:length(model.crp$best)

indice <- rep(id * (id + 1) / 2)

se <- sqrt(model.crp$V[indice])
coefficient <- coef(model.crp)
var.names <- names(coefficient)

df <- data.frame(coefficient = coefficient, se = se)
df <- df[1:(3 * (n.clust - 1)), ]

df$lower <- df$coefficient - (qnorm(0.975) * df$se)
df$upper <- df$coefficient + (qnorm(0.975) * df$se)
df$labels <- names(coefficient)[1:(3 * (n.clust - 1))]

df$labels <- factor(df$labels, levels = rev(df$labels))

df$sig <- FALSE

df[df$upper < 0, "sig"] <- TRUE
df[df$lower > 0, "sig"] <- TRUE

df <- df[n.clust:(3 * (n.clust - 1)), ] # Remove intercept

df %>%
  ggplot(aes(
    x = labels,
    y = coefficient,
    ymin = lower,
    ymax = upper,
    color = sig
  )) +
  geom_errorbar() +
  geom_point(size = 3.5) +
  geom_hline(yintercept = 0, lty = 2) +
  coord_flip() + # flip coordinates (puts labels on y axis)
  xlab("") +
  ylab("Estimate (95% CI)") +
  theme_bw() +
  ylim(-4, 4) +
  scale_color_manual(values = c("black", "#FF007F")) +
  theme(legend.position = "none")
```

```{R}
#| label: tbl-crp-ibd
#| fig-cap: "Table of β estimates and associated p-values obtained via univariate Wald tests."
wald <- c()
p.vals <- c()

for (index in seq(n.clust, 3 * (n.clust - 1))) {
  temp <- WaldMult(model.crp, pos = index)
  wald <- c(wald, temp[, 1])
  p.vals <- c(p.vals, temp[, 2])
}

df$wald <- wald
df$p.vals <- p.vals

rownames(df) <- df$labels

df %>%
  select(coefficient, se, wald, p.vals) %>%
  knitr::kable(digits = 3, col.names = c("Estimate", "SE", "Wald", "p-value"))
```

```{R}
dict.crp$diagnosis <- factor(dict.crp$diagnosis,
  levels = c(
    "Crohn's Disease",
    "Ulcerative Colitis",
    "IBDU"
  )
)

dict.crp$class <- factor(dict.crp$class,
  levels = seq(1, 8),
  labels = paste0("CRP ", seq(1, 8))
)

p2 <- dict.crp %>%
  ggplot() +
  geom_mosaic(aes(x = product(class), fill = diagnosis), show.legend = FALSE) +
  scale_fill_manual(values = c("#F8766D", "#619CFF", "#00BA38")) +
theme_minimal() +
  labs(
    x = "Cluster",
    y = "IBD type"
  ) +
  ggtitle(
    "IBD type",
    "CRP clustering"
  )

ggsave("plots/mosaic/crp-by-diag.png",
       p2,
       width = 7.5,
       height = 5,
       units = "in",
       dpi = 300,
       create.dir = TRUE
)
ggsave("plots/mosaic/crp-by-diag.pdf", p2, width = 7.5, height = 5)

p6 <- dict.crp %>%
  ggplot(aes(x = class, fill = diagnosis, color = diagnosis)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("#F8766D", "#619CFF", "#00BA38")) +
  scale_color_manual(values = c("#CF544B", "#417BD2", "#01932A")) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Cluster",
    y = "",
    fill = "IBD type",
    color = "IBD type"
  )
p6
```

```{R}
#| label: fig-crp-ibd-facet
#| fig-cap: "Percentage of (A) Crohn's disease, (B) ulcerative colitis, and (C) IBDU within each CRP cluster. Dotted lines indicate the overall percentage of each IBD type."
perc.table <- data.frame(
  diagnosis = character(),
  cluster = character(),
  perc = numeric()
)
for (g in 1:8) {
  temp.1 <- subset(dict.crp, class == paste0("CRP ", g))
  for (diag in c("Crohn's Disease", "Ulcerative Colitis", "IBDU")) {
    perc <- nrow(subset(temp.1, diagnosis == diag)) / nrow(temp.1)
    perc.table <- rbind(
      perc.table,
      data.frame(
        diagnosis = diag,
        cluster = paste0("CRP ", g),
        perc = perc
      )
    )
  }
}

totalCdPerc <- nrow(subset(dict.crp, diagnosis == "Crohn's Disease")) /
  nrow(dict.crp)

p1 <- perc.table %>%
  filter(diagnosis == "Crohn's Disease") %>%
  ggplot(aes(x = cluster, y = perc)) +
  geom_bar(stat = "identity", fill = "#F8766D", color = "#CF544B") +
  geom_hline(yintercept = totalCdPerc, linetype = "dashed", color = "#4D4730") +
  theme_minimal() +
  labs(x = "Cluster", y = "Percentage")

totalCdPerc <- nrow(subset(dict.crp, diagnosis == "Ulcerative Colitis")) /
  nrow(dict.crp)

p2 <- perc.table %>%
  filter(diagnosis == "Ulcerative Colitis") %>%
  ggplot(aes(x = cluster, y = perc)) +
  geom_bar(stat = "identity", fill = "#619CFF", color = "#417BD2") +
  geom_hline(yintercept = totalCdPerc, linetype = "dashed", color = "#4D4730") +
  theme_minimal() +
  labs(x = "Cluster", y = "Percentage")


totalCdPerc <- nrow(subset(dict.crp, diagnosis == "IBDU")) /
  nrow(dict.crp)

p3 <- perc.table %>%
  filter(diagnosis == "IBDU") %>%
  ggplot(aes(x = cluster, y = perc)) +
  geom_bar(stat = "identity", fill = "#00BA38", color = "#01932A") +
  geom_hline(yintercept = totalCdPerc, linetype = "dashed", color = "#4D4730") +
  theme_minimal() +
  labs(x = "Cluster", y = "Percentage")

p <- p1 / p2 / p3 + patchwork::plot_annotation(tag_levels = "A")  &
  scale_y_continuous(labels = scales::percent, limits = c(0, 1))
ggsave("plots/ibd-type-crp.png", p, width = 7.5, height = 7.5)
ggsave("plots/ibd-type-crp.pdf", p, width = 7.5, height = 7.5)
print(p)
```

### Advanced therapy in CD

```{R}
dict.crp.cd <- readRDS(paste0(prefix, "processed/dict-crp-cd.RDS"))

dict.crp.cd <- merge(dict.crp.cd,
                     dict.crp[, c("ids", "class")],
                     by = "ids",
                     all.x = TRUE,
                     all.y = FALSE
)
```



```{R}
pred.crp.df <- data.frame(
  crp_time = c(seq(0, 6.25, 0.01), seq(0, 6.25, 0.01)),
  diagnosis = c(
    rep("Crohn's Disease", 626),
    rep("Ulcerative Colitis", 626)
  )
)
pred.crp.df.update <- predictY(model.crp, pred.crp.df,
                               var.time = "crp_time",
                               draws = TRUE
)$pred

temp <- data.frame(Time = NULL, Cluster = NULL, Value = NULL)

for (g in 1:8) {
  temp <- rbind(
    temp,
    data.frame(
      Time = seq(0, 6.25, 0.01),
      Cluster = g,
      Value = pred.crp.df.update[, g]
    )
  )
}


dict.crp.cd$timeToTherapy <- with(dict.crp.cd, startDate - date.of.diag) / 365.25

for (g in 1:8) {
  p1 <- temp %>%
    filter(Cluster == g) %>%
    ggplot(aes(x = Time, y = Value)) +
    geom_line(color = "#558C8C") +
    theme_minimal() +
    ylim(0, 6) +
    xlab("") +
    ylab("Log (CRP (\u03BCg/mL))")
  p2 <- dict.crp.cd %>%
    filter(class == paste0("CRP ", g)) %>%
    filter(timeToTherapy > 0) %>%
    ggplot(aes(x = timeToTherapy)) +
    geom_density(color = "#5D002F", fill = "#82204A") +
    theme_minimal() +
    xlim(0, 7) +
    ylim(0, NA) +
    xlab("Time") +
    ylab("Density of time to first advanced therapy")
  p <- p1 / p2 + plot_annotation(
    title = paste("CRP ", g),
    subtitle = "Crohn's disease",
    tag_levels = "A"
  )
  ggsave(paste0("plots/at-density-crp/g=", g, ".png"),
         p,
         width = 7,
         height = 7,
         units = "in",
         create.dir = TRUE
  )
  print(p)
}
```

```{R}
temp <- dict.crp.cd %>%
  mutate(timeToTherapy = if_else(advancedT == 1, 7, as.numeric(timeToTherapy))) %>%
  mutate(advancedT = if_else(advancedT == 1, 0, 1)) %>%
  mutate(advancedT = if_else(timeToTherapy > 7, 0, 1)) %>%
  mutate(timeToTherapy = if_else(timeToTherapy > 7, 7, timeToTherapy))
survfit(Surv(timeToTherapy, advancedT) ~ class, data = temp) %>%
  ggsurvplot(pval = TRUE)
```

```{R}
#| echo: false
png("plots/at-survival-crp-ci.png",
    width = 7,
    height = 5,
    units = "in",
    res = 300
)
survfit(Surv(timeToTherapy, advancedT) ~ class, data = temp) %>%
  ggsurvplot(pval = TRUE, conf.int = TRUE)
invisible(dev.off())
```


```{R}
p <- dict.crp.cd %>%
  ggplot() +
  geom_mosaic(aes(x = product(class), fill = firstY), show.legend = FALSE) +
  theme_minimal() +
  labs(
    x = "Cluster",
    y = "Advanced therapy within first year of diagnosis"
  ) +
  ggtitle(
    "Crohn's disease subjects",
    "CRP clustering"
  )
ggsave("plots/mosaic/crp-cd.png",
       p,
       width = 7.5,
       height = 5,
       units = "in",
       dpi = 300,
       create.dir = TRUE
)
ggsave("plots/mosaic/crp-cd.pdf", p, width = 7.5, height = 5)
p
```

```{R}
pred.crp.df <- data.frame(
  crp_time = c(seq(0, 7, 0.01), seq(0, 7, 0.01)),
  diagnosis = c(
    rep("Crohn's Disease", 701),
    rep("Ulcerative Colitis", 701)
  )
)
pred.crp.df.update <- predictY(model.crp, pred.crp.df,
                               var.time = "crp_time",
                               draws = TRUE
)$pred

temp <- data.frame(Time = NULL, Cluster = NULL, Value = NULL)

for (g in 1:8) {
  temp <- rbind(
    temp,
    data.frame(
      Time = seq(0, 7, 0.01),
      Cluster = g,
      Value = pred.crp.df.update[, g]
    )
  )
}

traj.crp.outlines <- list()

for (g in 1:8) {
  traj.crp.outlines[[g]] <- temp %>%
    filter(Cluster == g) %>%
    ggplot(aes(x = Time, y = Value)) +
    geom_line(color = "#757575") +
    theme_classic() +
    theme(
      text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(color = "#dfdfdf")
    ) +
    ylim(0, 6)
}


p1 <- dict.crp.cd %>%
  ggplot(aes(x = class, fill = firstY, color = firstY)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("#F8766D", "#619CFF")) +
  scale_color_manual(values = c("#CF544B", "#417BD2")) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Cluster",
    y = "",
    fill = desc,
    color = desc
  )


p <- (p1) /
  (traj.crp.outlines[[1]] +
     traj.crp.outlines[[2]] +
     traj.crp.outlines[[3]] +
     traj.crp.outlines[[4]] +
     traj.crp.outlines[[5]] +
     traj.crp.outlines[[6]] +
     traj.crp.outlines[[7]] +
     traj.crp.outlines[[8]] +
     plot_layout(ncol = 8, guides = "collect", widths = unit(2.8, "cm"))) +
  plot_layout(heights = c(4, 1), nrow = 2)

ggsave("plots/mosaic/crp-by-AT-with-traj.pdf",
       p,
       width = 7 * 1.65,
       height = 5 * 1.65
)
ggsave("plots/mosaic/crp-by-AT-with-traj.png",
       p,
       width = 7 * 1.65,
       height = 5 * 1.65
)
```

### Mutlivariate modelling

```{R}
mlr <- multinom(class ~ age + sex + diagnosis, data = dict.crp, trace = FALSE)
z <- summary(mlr)$coefficients / summary(mlr)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2
knitr::kable(p)
```

#### Overall cluster-specific trajectories

Here, we extract overall cluster-specific trajectories as these will be used
for visualisation purposes in order to better understand patterns of AT use.
Note that model outputs do not match the reordered clusters (based on cumulative
                                                             inflammation) used throughout this report. As such, we use `title.mapping` to
re-order the trajectories when these are plotted.

```{r cluster specific trajectories}
time.pred <- seq(0, 7, by = 0.01)

pred.fc.df <- data.frame(
  calpro_time = c(time.pred, time.pred),
  diagnosis = c(
    rep("Crohn's Disease", length(time.pred)),
    rep("Ulcerative Colitis", length(time.pred))
  )
)
pred.fc.df.update <- lcmm::predictY(model.fc,
                                    pred.fc.df,
                                    var.time = "calpro_time",
                                    draws = TRUE
)$pred

pred <- predictY(model.fc, pred.fc.df, var.time = "calpro_time", draws = TRUE)$pred

pred <- as.data.frame(pred[1:length(time.pred), ])
pred$time <- time.pred

ylimit <- log(2500)
title.mapping <- c(7, 6, 4, 8, 1, 5, 2, 3)
```

This is based on the first AT prescription for each patient.

W/ Density

```{R}
#| column: page
#| fig-width: 12
#| fig-height: 12
dens.cd.big <- data.frame(x = numeric(),
                          y = numeric(),
                          CumulativeY = numeric(),
                          class = character())

dens.uc.big <- data.frame(x = numeric(),
                          y = numeric(),
                          CumulativeY = numeric(),
                          class = character())

km.df <- data.frame(time = numeric(),
                    cumhaz = numeric(),
                    class = character(),
                    diag = character())

for (g in 1:8) {

  # Calculate cumulative patterns
  temp.cd <- myDF.fc %>%
    filter(class_order == paste0("FC", g)) %>%
    filter(diagnosis == "Crohn's Disease")

  km <- survfit(Surv(AT_line_1_cens, AT_7Y) ~ 1, data = temp.cd)
  km.df <- rbind(km.df,
                 data.frame(time = km$time,
                            cumhaz =  1 - km$surv,
                            class = paste0("FC", g),
                            diag = "CD"))

  temp.uc <- myDF.fc %>%
    filter(class_order == paste0("FC", g)) %>%
    filter(diagnosis == "Ulcerative Colitis")

  km <- survfit(Surv(AT_line_1_cens, AT_7Y) ~ 1, data = temp.uc)
  km.df <- rbind(km.df,
                 data.frame(time = km$time,
                            cumhaz =  1 - km$surv,
                            class = paste0("FC", g),
                            diag = "UC"))

  # Calculate densities
  temp.cd <- temp.cd %>%
    filter(AT_7Y == 1)
  dens.cd <- density(temp.cd$AT_line_1, kernel = "gaussian", from = 0, to = 7, na.rm = TRUE)
  dens.cd <- data.frame(x = dens.cd$x, y = dens.cd$y)
  dens.cd$CumulativeY <- cumsum(dens.cd$y)
  dens.cd$class <- paste0("FC", g)
  dens.cd.big <- rbind(dens.cd.big,
                       dens.cd)

  temp.uc <- temp.uc %>%
    filter(AT_7Y == 1)
  dens.uc <- density(temp.uc$AT_line_1, kernel = "gaussian", from = 0, to = 7, na.rm = TRUE)
  dens.uc <- data.frame(x = dens.uc$x, y = dens.uc$y)
  #  dens.uc$y <- dens.uc$y/nrow(temp.uc)
  dens.uc$CumulativeY <- cumsum(dens.uc$y)
  dens.uc$class <- paste0("FC", g)
  dens.uc.big <- rbind(dens.uc.big,
                       dens.uc)
}

dens.cd.big$diag <- "CD"
dens.uc.big$diag <- "UC"

dens <- rbind(dens.cd.big, dens.uc.big)


p1 <- ggplot(dens, aes(x = x, y = y, color = diag, fill = diag)) +
  geom_area(alpha = .5, position = 'identity') +
  geom_line(aes(y = cumhaz, x = time), data = km.df, lty = 2) +
  facet_wrap(~class, ncol = 2) +
  theme_minimal()
p1
```

This is based on the first AT prescription for each patient.

W/ Density

```{R}
#| column: page
#| fig-width: 12
#| fig-height: 12
dens.cd.big <- data.frame(x = numeric(),
                          y = numeric(),
                          CumulativeY = numeric(),
                          class = character())

dens.uc.big <- data.frame(x = numeric(),
                          y = numeric(),
                          CumulativeY = numeric(),
                          class = character())

km.df <- data.frame(time = numeric(),
                    cumhaz = numeric(),
                    class = character(),
                    diag = character())

for (g in 1:8) {

  # Calculate cumulative patterns
  temp.cd <- myDF.crp %>%
    filter(class_order == paste0("CRP", g)) %>%
    filter(diagnosis == "Crohn's Disease")

  km <- survfit(Surv(AT_line_1_cens, AT_7Y) ~ 1, data = temp.cd)
  km.df <- rbind(km.df,
                 data.frame(time = km$time,
                            cumhaz =  1 - km$surv,
                            class = paste0("CRP", g),
                            diag = "CD"))

  temp.uc <- myDF.crp %>%
    filter(class_order == paste0("CRP", g)) %>%
    filter(diagnosis == "Ulcerative Colitis")

  km <- survfit(Surv(AT_line_1_cens, AT_7Y) ~ 1, data = temp.uc)
  km.df <- rbind(km.df,
                 data.frame(time = km$time,
                            cumhaz =  1 - km$surv,
                            class = paste0("CRP", g),
                            diag = "UC"))

  # Calculate densities
  temp.cd <- temp.cd %>%
    filter(AT_7Y == 1)
  dens.cd <- density(temp.cd$AT_line_1, kernel = "gaussian", from = 0, to = 7, na.rm = TRUE)
  dens.cd <- data.frame(x = dens.cd$x, y = dens.cd$y)
  dens.cd$CumulativeY <- cumsum(dens.cd$y)
  dens.cd$class <- paste0("CRP", g)
  dens.cd.big <- rbind(dens.cd.big,
                       dens.cd)

  temp.uc <- temp.uc %>%
    filter(AT_7Y == 1)
  dens.uc <- density(temp.uc$AT_line_1, kernel = "gaussian", from = 0, to = 7, na.rm = TRUE)
  dens.uc <- data.frame(x = dens.uc$x, y = dens.uc$y)
  #  dens.uc$y <- dens.uc$y/nrow(temp.uc)
  dens.uc$CumulativeY <- cumsum(dens.uc$y)
  dens.uc$class <- paste0("CRP", g)
  dens.uc.big <- rbind(dens.uc.big,
                       dens.uc)
}

dens.cd.big$diag <- "CD"
dens.uc.big$diag <- "UC"

dens <- rbind(dens.cd.big, dens.uc.big)


p1 <- ggplot(dens, aes(x = x, y = y, color = diag, fill = diag)) +
  geom_area(alpha = .5, position = 'identity') +
  geom_line(aes(y = cumhaz, x = time), data = km.df, lty = 2) +
  facet_wrap(~class, ncol = 2) +
  theme_minimal()
p1
```

w/ trajectories

```{R}
#| column: page
#| fig-width: 12
#| fig-height: 12
km.df <- data.frame(time = numeric(),
                    cumhaz = numeric(),
                    class = character(),
                    diag = character())

traj.df <- data.frame(prediction = numeric(),
                      time = numeric(),
                      class = character())


for (g in 1:8) {
  traj <- pred[, c(g, ncol(pred)), ]
  names(traj) <- c("prediction", "time")
  traj$class <- paste0("FC", g)
  traj.df <- rbind(traj.df, traj)

  temp.all <- myDF.fc %>%
    filter(class_order == paste0("FC", g))

  # Calculate cumulative patterns
  temp.cd <- myDF.fc %>%
    filter(class_order == paste0("FC", g)) %>%
    filter(diagnosis == "Crohn's Disease")

  km <- survfit(Surv(AT_line_1_cens, AT_7Y) ~ 1, data = temp.cd)
  km.df <- rbind(km.df,
                 data.frame(time = km$time,
                            cumhaz =  1 - km$surv,
                            class = paste0("FC", g),
                            diag = "CD"))

  temp.uc <- myDF.fc %>%
    filter(class_order == paste0("FC", g)) %>%
    filter(diagnosis == "Ulcerative Colitis")

  km <- survfit(Surv(AT_line_1_cens, AT_7Y) ~ 1, data = temp.uc)
  km.df <- rbind(km.df,
                 data.frame(time = km$time,
                            cumhaz =  1 - km$surv,
                            class = paste0("FC", g),
                            diag = "UC"))

  temp.all <- myDF.fc %>%
    filter(class_order == paste0("FC", g))

  km <- survfit(Surv(AT_line_1_cens, AT_7Y) ~ 1, data = temp.all)
  km.df <- rbind(km.df,
                 data.frame(time = km$time,
                            cumhaz =  1 - km$surv,
                            class = paste0("FC", g),
                            diag = "All"))
}

traj.df$class <- plyr::mapvalues(traj.df$class,
                                 from = paste0("FC", seq(1, 8)),
                                 to = paste0("FC", c(7, 6, 4, 8, 1, 5, 2, 3)))

p1 <- km.df %>%
  ggplot(aes(x = time, y = cumhaz)) +
  geom_line(aes(color = diag), lty = 1) +
  geom_line(aes(x = time, y = prediction/log(2500) ), data = traj.df, color = "black") +
  facet_wrap(~class, ncol = 2) +
  theme_minimal() +
  scale_y_continuous(name = "Cumulative",

                     sec.axis = sec_axis(transform = ~.*log(2500),
                                         name = "Log (FC (\u03BCg/g))"))
p1

p2 <- km.df %>%
  subset(diag != "All") %>%
  ggplot(aes(x = time, y = cumhaz)) +
  geom_line(aes(color = diag), lty = 1) +
  geom_line(aes(x = time, y = prediction/log(2500) ), data = traj.df, color = "black") +
  facet_wrap(~class, ncol = 2) +
  theme_minimal() +
  scale_y_continuous(name = "Cumulative",

                     sec.axis = sec_axis(transform = ~.*log(2500),
                                         name = "Log (FC (\u03BCg/g))"))
p2

p3 <- km.df %>%
  subset(diag != "All") %>%
  ggplot(aes(x = time, y = cumhaz)) +
  geom_line(aes(color = diag), lty = 1, lwd = 1.2) +
  facet_wrap(~class, ncol = 2) +
  theme_minimal() +
  scale_y_continuous(name = "Cumulative")
p3
```
