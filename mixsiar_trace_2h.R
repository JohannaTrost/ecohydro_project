# Brian Stock
# March 8, 2016
# Script file to run wolves example without GUI

# Load MixSIAR package
library(MixSIAR)
library(readxl)
library(readr)
library(dplyr)
library(TAF)
library(readr)
library(tidyr)
library(ggplot2)

get_summary <- function (jags.1, mix, source, model_dir) {

  combined <- combine_sources(jags.1, mix, source, alpha.prior=1,
                            groups=list(Label=c("Label"),
                                        Precipitation=c("Precipitation"),
                                        Soil = c("Soil")
                            ))

  summary_output_string <- paste(capture.output(summary_stat(combined)), collapse = "\n")

  # Remove head
  summary_output_string <- strsplit(summary_output_string, "#\n\n")[[1]][[2]]

  # Split the input string into lines
  lines <- strsplit(summary_output_string, "\n")[[1]]

  # Remove the header line
  lines <- lines[-1]

  # Parse each line
  parsed_data <- lapply(lines, function(line) {
    parts <- strsplit(line, "\\s+")[[1]]
    parts <- parts[parts != ""]

    variable <- parts[1]
    values <- as.numeric(parts[2:length(parts)])

    data.frame(
      variable = variable,
      Mean = values[1],
      SD = values[2],
      `2.5%` = values[3],
      `25%` = values[4],
      `50%` = values[5],
      `75%` = values[6],
      `97.5%` = values[7],
      stringsAsFactors = FALSE
    )
  })

  # Combine all parsed lines into a single data frame
  df <- do.call(rbind, parsed_data)

  # Set column names
  colnames(df) <- c("variable", "Mean", "SD", "2.5%", "25%", "50%", "75%", "97.5%")

  # Parse variable names into treeid, level, and source
  df <- df %>%
    separate(variable, into = c("p", "source", "level", "treeid"), sep = "\\.", extra = "drop") %>%
    filter(!is.na(treeid)) %>%  # Remove rows where parsing failed (e.g., Epsilon rows)
    select(-p) %>%
    tibble(.)

  write_csv(df, file.path(model_dir, "summary_statistics.csv"))

  return(df)
}


################################################################################
# Load mixture data, i.e. your:
#    Consumer isotope values (trophic ecology / diet)
#    Mixed sediment/water tracer values (sediment/hydrology fingerprinting)

# 'filename' - name of the CSV file with mix/consumer data
# 'iso_names' - column headings of the tracers/isotopes you'd like to use
# 'random_effects' - column headings of any random effects
# 'cont_effects' - column headings of any continuous effects

# 'iso_names', 'random_effects', and 'cont_effects' can be a subset of your columns
#   i.e. have 3 isotopes in file but only want MixSIAR to use 2,
#   or have data by Region and Pack but only want MixSIAR to use Region

indir <- "/Users/johanna/Uni/ws_25_26/ecohydrology/ecohydro_project/data/Water stable isotope analysis"
outdir <- "/Users/johanna/Uni/ws_25_26/ecohydrology/ecohydro_project/output"

# Read data
filename <- file.path(indir, "20260309_water_stable_isotope_data_corrected.csv")
data <- read_csv(filename)

data <- data %>%
  rename("iso2H" = "d.D_H.corr.drift", "iso18O" = "d.18_16.corr.drift") %>%
  mutate(Type = case_when(
      substr(Identifier_1, 1, 2) == "SD" ~ "Soil",
      substr(Identifier_1, 1, 2) == "SK" ~ "Soil",
      substr(Identifier_1, 1, 1) == "D"  ~ "Tree",
      substr(Identifier_1, 1, 1) == "K"  ~ "Tree",
      substr(Identifier_1, 1, 1) == "P"  ~ "Precipitation", # Puddle
      substr(Identifier_1, 1, 1) == "C"  ~ "Label", # Control label
      TRUE ~ NA_character_
    ),
     treeid = case_when(
      substr(Identifier_1, 1, 1) == "D"  ~ substr(Identifier_1, 1, 2),
      substr(Identifier_1, 1, 1) == "K"  ~ substr(Identifier_1, 1, 2),
      substr(Identifier_1, 1, 2) == "SD"  ~ substr(Identifier_1, 2, 3),
      substr(Identifier_1, 1, 2) == "SK"  ~ substr(Identifier_1, 2, 3),
     ))

# Generate and save mixing data df
mix_df <- data %>%
  filter(Type == "Tree", Date != "2026-02-25")
mix_df[substr(mix_df$Identifier_1, 3, 4) == "_C", "level"] <- "mid"

for (day in unique(mix_df$Date)) {
    for (spec in unique(mix_df$species)) {
        mix.filename <- file.path(outdir, paste0("mix_data_", spec, "_", as.character(as.Date(day)), ".csv"))
        mix_df_curr <- mix_df %>%
          filter(Date == day, species == spec) %>%
          select(iso2H, iso18O, level, treeid)
        write_csv(mix_df_curr, mix.filename)
    }
}

# Generate source data and save
source_df <- data %>%
  filter(Type %in% c("Precipitation", "Soil", "Label")) %>%
  select(Type, species, iso2H, iso18O, treeid)

for (spec in unique(mix_df$species)) {
    source.filename <- file.path(outdir, paste0("source_data_", spec, ".csv"))
    src_df_curr <- source_df %>%
      filter(ifelse(Type == "Soil", species == spec, TRUE))
    write_csv(src_df_curr, source.filename)
}
write_csv(source_df, source.filename)

# Generate discrimination data
discr.filename <- file.path(outdir, "discr_data.csv")
discr_df <- data.frame(
  Type = c("Precipitation", "Soil", "Label"),
  Meaniso2H = 0,
  SDiso2H = 0,
  Meaniso18O = 0,
  SDiso18O = 0
)
write_csv(discr_df, discr.filename)

# --------------------------Run model for each day and species -------------------------- #

for (spec in c(unique(mix_df$species)[2])) {
    source.filename <- file.path(outdir, paste0("source_data_", spec, ".csv"))

    for (day in unique(mix_df$Date)) {
        model_id <- paste(spec, as.character(as.Date(day)), sep = "_")
        mix.filename <- file.path(
          outdir, paste0("mix_data_", model_id , ".csv")
        )

        # Load mixture data
        mix <- load_mix_data(filename=mix.filename,
                             iso_names=c("iso2H", "iso18O"),
                             factors=c("treeid", "level"),
                             fac_random=c(FALSE, FALSE),
                             fac_nested=c(FALSE, FALSE),
                             cont_effects=NULL)

        ################################################################################
        # Load source data, i.e. your:
        #    Source isotope values (trophic ecology / diet)
        #    Sediment/water source tracer values (sediment/hydrology fingerprinting)

        # 'filename': name of the CSV file with source data
        # 'source_factors': column headings of random/fixed effects you have source data by
        # 'conc_dep': TRUE or FALSE, do you have concentration dependence data in the file?
        # 'data_type': "means" or "raw", is your source data as means+SD, or do you have raw data

        # Load source data
        source <- load_source_data(filename=source.filename,
                                   source_factors=NULL,
                                   conc_dep=FALSE,
                                   data_type="raw",
                                   mix)

        ################################################################################
        # Load discrimination data, i.e. your:
        #    Trophic Enrichment Factor (TEF) / fractionation values (trophic ecology/diet)
        #    xxxxxxxx (sediment/hydrology fingerprinting)

        # 'filename' - name of the CSV file with discrimination data

        # To run on your data, replace the system.file call with the path to your file
        #discr.filename <- system.file("extdata", "wolves_discrimination.csv", package = "MixSIAR")

        discr <- load_discr_data(filename=discr.filename, mix)

        #####################################################################################
        # Make isospace plot
        # Are the data loaded correctly?
        # Is your mixture data in the source polygon?
        # Are one or more of your sources confounded/hidden?

        # 'filename' - name you'd like MixSIAR to save the isospace plot as
        #              (extension will be added automatically)
        # 'plot_save_pdf' - TRUE or FALSE, should MixSIAR save the plot as a .pdf?
        # 'plot_save_png' - TRUE or FALSE, should MixSIAR save the plot as a .png?

        # plot_data(filename=file.path(outdir, "isospace_plot"),
        #           plot_save_pdf=TRUE,
        #           plot_save_png=FALSE,
        #           mix,source,discr)

        # If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
        #   *Note 1: discrimination SD is added to the source SD (see calc_area.R for details)
        #   *Note 2: If source data are by factor (as in wolf ex), computes area for each polygon
        #             (one for each of 3 regions in wolf ex)
        if(mix$n.iso==2) calc_area(source=source,mix=mix,discr=discr) # 20769.94

        ################################################################################
        # Define your prior, and then plot using "plot_prior"
        #   RED = your prior
        #   DARK GREY = "uninformative"/generalist (alpha = 1)
        #   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

        # default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
        plot_prior(alpha.prior=1,source)

        ################################################################################
        # Write JAGS model file (define model structure)
        # Model will be saved as 'model_filename' ("MixSIAR_model.txt" is default,
        #    but may want to change if in a loop)

        # There are 3 error term options available:
        #   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
        #   2. Residual only (resid_err = TRUE, process_err = FALSE)
        #   3. Process only (resid_err = FALSE, process_err = TRUE)

        # 'model_filename': don't need to change, unless you create many different models
        # 'resid_err': include residual error in the model?
        # 'process_err': include process error in the model?

        #  *Note: If you have only 1 mix datapoint, you have no information about the
        #         mixture/consumer variability. In this case, we ues the original MixSIR
        #         error model (which does not fit a residual error term).
        #         This is the same behavior as 'siarsolo' in SIAR.

        model_filename <- file.path(outdir, paste0("MixSIAR_model", model_id, ".txt"))
        resid_err <- TRUE
        process_err <- TRUE
        write_JAGS_model(model_filename, resid_err, process_err, mix, source)

        ################################################################################
        # Run model
        # JAGS output will be saved as 'jags.1'

        # MCMC run options:
        # run <- "test"       # chainLength=1000, burn=500, thin=1, chains=3, calcDIC=TRUE
        # run <- "very short" # chainLength=10000, burn=5000, thin=5, chains=3, calcDIC=TRUE
        # run <- "short"      # chainLength=50000, burn=25000, thin=25, chains=3, calcDIC=TRUE
        # run <- "normal"     # chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE
        # run <- "long"       # chainLength=300000, burn=200000, thin=100, chains=3, calcDIC=TRUE
        # run <- "very long"  # chainLength=1000000, burn=500000, thin=500, chains=3, calcDIC=TRUE
        # run <- "extreme"    # chainLength=3000000, burn=1500000, thin=500, chains=3, calcDIC=TRUE

        # Can also set custom MCMC parameters
        # run <- list(chainLength=200000, burn=150000, thin=50, chains=3, calcDIC=TRUE)

        # Good idea to use 'test' first to check if
        #   1) the data are loaded correctly, and
        #   2) the model is specified correctly
        jags.1 <- run_model(run="test", mix, source, discr, model_filename, alpha.prior = 1)

        # After a test run works, increase the MCMC run to a value that may converge
        # jags.1 <- run_model(run="normal", mix, source, discr, model_filename, alpha.prior = 1)

        ################################################################################
        # Process JAGS output

        # Choose output options (see ?output_options for details)
        output_options <- list(summary_save = TRUE,
                               summary_name = "summary_statistics",
                               sup_post = FALSE,
                               plot_post_save_pdf = TRUE,
                               plot_post_name = "posterior_density",
                               sup_pairs = FALSE,
                               plot_pairs_save_pdf = TRUE,
                               plot_pairs_name = "pairs_plot",
                               sup_xy = TRUE,
                               plot_xy_save_pdf = FALSE,
                               plot_xy_name = "xy_plot",
                               gelman = TRUE,
                               heidel = FALSE,
                               geweke = TRUE,
                               diag_save = TRUE,
                               diag_name = "diagnostics",
                               indiv_effect = FALSE,
                               plot_post_save_png = FALSE,
                               plot_pairs_save_png = FALSE,
                               plot_xy_save_png = FALSE,
                               diag_save_ggmcmc = FALSE,
                               return_obj = TRUE)

        # Create diagnostics, summary statistics, and posterior plots
        model_dir <- file.path(outdir, model_id)
        mkdir(model_dir)
        setwd(model_dir)
        output_JAGS(jags.1, mix, source, output_options)
        setwd("/Users/johanna/Uni/ws_25_26/ecohydrology/ecohydro_project")

        get_summary(jags.1, mix, source, model_dir)
    }
}

# Colletct results
models_summary <- NULL
for (spec in unique(mix_df$species)) {
    for (day in unique(mix_df$Date)) {
        model_id <- paste(spec, as.character(as.Date(day)), sep = "_")
        curr_model_summary <- read_csv(file.path(outdir, model_id, "summary_statistics.csv"))
        curr_model_summary$Date <- as.Date(day)
        curr_model_summary$species <- spec

        if (is.null(models_summary)) {
          models_summary <- curr_model_summary
        } else {
          models_summary <- rbind(models_summary, curr_model_summary)
        }
    }
}
write_csv(models_summary, file.path(outdir, "summary_models.csv"))

# -------------------------- Plot model results -------------------------- #

models_summary <- models_summary %>%
  mutate(treeid2 = as.numeric(substr(treeid, 2, 2)))

col_precip <- "#0071A6"
col_soil <- "#B29060"
col_label <- "#BD3977"

for (day in unique(models_summary$Date)) {
  pdf(file.path(outdir, paste0(as.Date(day), "_mixsiar_contribution.pdf")), width = 10,
      height = 5)
  ggplot(
    models_summary %>% filter(Date == as.Date(day)),
    aes(x = level, y = Mean, fill = source)
  ) +
    facet_grid(species ~ treeid2, scales = "free_x", space = "free_x") +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("Precipitation" = col_precip,
                                 "Soil" = col_soil, "Label" = col_label)) +
    labs(x = "Level", y = "Contribution [%]", fill = "Source") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.spacing = unit(1, "lines")) # Adjust spacing between facets
  dev.off()
}

# -- on log scale

day <- unique(models_summary$Date)[2]
# Adjust Mean to avoid log(0)
models_summary_adjusted <- models_summary %>%
  mutate(Mean_adj = Mean + 0.001) %>%
  filter(Date == as.Date(day))

pdf(file.path(outdir, paste0(as.Date(day), "_mixsiar_contribution_log.pdf")), width = 10, height = 5)

ggplot(models_summary_adjusted, aes(x = level, y = Mean_adj, fill = source)) +
  facet_grid(species ~ treeid2, scales = "free_x", space = "free_x") +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_log10(labels = scales::percent) +
  scale_fill_manual(values = c("Precipitation" = col_precip, "Soil" = col_soil, "Label" = col_label)) +
  labs(x = "Level", y = "Contribution [Log %]", fill = "Source") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(1, "lines"))

dev.off()
