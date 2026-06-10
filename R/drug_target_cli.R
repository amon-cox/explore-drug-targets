#!/usr/bin/env Rscript

## Load required packages
library(optparse) # argument parsing for CLI
suppressMessages(library(tidyverse)) # convenient data manipulation
library(hgnc) # access HUGO Gene Nomenclature Committee symbols

## Define command-line options
options_list <- list( # specify options for an `optparse`-based command-line interface
    make_option(
        opt_str = "--target",
        type = "character",
        dest = "target",
        help = "An Ensembl gene ID (ENSG###...) or HUGO symbol of the target [required]",
        required = TRUE), # will stop if no target supplied
    make_option(
        opt_str = "--minPhase",
        type = "integer",
        dest = "minPhase",
        default = 3,
        help = "Lowest clinical trial phase to include, in addition to approved drugs [default: %default]"),
    make_option(
        opt_str = "--top",
        type = "integer",
        dest = "top",
        default = 10,
        help = "Number of top-ranked compounds to print [default: %default]"),
    make_option(
        opt_str = "--output",
        type = "character",
        dest = "output",
        default = NA,
        help = "Include a custom file name for results. Do not include the file extension. [optional]"),
    make_option(
        opt_str = "--plot",
        type = "logical",
        action = "store_true",
        dest = "plot",
        default = FALSE,
        help = "Save a scatterplot fo the ranked compounds [default: %default]"
    )
)

parser <- OptionParser( # set up the `optparse` CLI parser
    option_list = options_list,
    usage = "Rscript %prog --target <ID> [options]",
    description = "A command-line interface to rank compounds acting on a specified target by selectivity and adverse effect profile."
)

## Retrieve command-line arguments
args <- parse_args(parser)

## Retrieve HUGO gen symbols form HGNC
message("\nFetching HGNC dataset... (~10s)")
hgnc_data <- import_hgnc_dataset()
    # downloads the current HGNC dataset on each run (~10s).
    # consider caching locally and check hgnc::last_update() for frequent use

### Convert --target input
ensembl_pattern <- "^ENSG[0-9]+"

if (str_detect(args$target, ensembl_pattern)) { # check if --target is Ensembl ID
    query_ensembl <- args$target
    query_symbol <- hgnc_data |>
        filter(ensembl_gene_id == query_ensembl) |>
        pull(symbol) |> # get HUGO gene symbol
        first()
} else { # assume HUGO gene symbol
    query_symbol  <- args$target
    query_ensembl <- hgnc_data |>
        filter(symbol == str_to_upper(args$target)) |> # forces uppercase
        pull(ensembl_gene_id) |>
        first()
}

if (is.na(query_ensembl)) { # error if Ensembl ID not found or submitted
    stop(paste("Could not resolve target to an Ensembl ID:", args$target), call. = FALSE)
}

## Source the data
molecules <- read_tsv(file.path("data", "molecules.tsv.gz"), show_col_types = FALSE)
adverseEffects <- read_tsv(file.path("data", "adverseEffects.tsv.gz"), show_col_types = FALSE)

## Reformat then filter molecules by command-line arguments
molecules_filtered <- molecules |>
    mutate(
        targetsList = linkedTargets |> # parse $linkedTargets into list
            str_remove_all("\\[|\\]") |> # remove brackets
            str_split(",\\s*") |> # split on comma and spaces
            map(~ .x[.x != ""]), # drop empty strings (becomes character(0))
        nTargets = map_int(targetsList, length) # count length per element
    ) |>
    filter(
        !hasBeenWithdrawn, # exclude withdrawn drugs
        isApproved | maximumClinicalTrialPhase %in% args$minPhase:4,
            # include approved drugs and those in a specified trial phase (default 3–4)
        map_lgl(targetsList, ~ query_ensembl %in% .x) # filter down to specified target
    )

## Join with adverse effects data, then summarize and rank
molecules_ranked <- molecules_filtered |>
    mutate(chembl_id = str_remove_all(id, "CHEMBL") |> as.double()) |> # make molecules $id column compatible with adverseEffects
    left_join(adverseEffects, by = "chembl_id") |>
    group_by(id, name, nTargets) |>
    summarize(
        hasAdverseData = !all(is.na(llr)), # report if matches in adverseEffects
        max_llr = if (all(is.na(llr))) 0 else max(llr, na.rm = TRUE),
            # max LLR is most conservative safety signal
            # if no adverse effects reported, retain as a candidate
        nAdverseEvents = sum(!is.na(event)), # count reported events, not NA rows
        .groups = "drop"
    ) |>
    arrange(
        nTargets,
        max_llr,
        nAdverseEvents
    )

## Write out results and print args$top to user
if (!is.na(args$output)) {
    res_name <- args$output
} else {
    res_name <- paste0("results ", query_symbol, " ", query_ensembl, " minPhase", args$minPhase, " ", Sys.Date())
}

message("\nWriting ranked molecules to results/", res_name, ".tsv")
write_tsv(molecules_ranked, file = file.path("results", paste0(res_name, ".tsv")))

message("\nPrinting the first ", args$top ," ranked compounds for target: ", query_symbol, " (", query_ensembl, ")")
print(molecules_ranked, n = args$top)

## Plot results if requested
if (args$plot) {
    library(ggrepel) # for label repulsion

    ### ggplot of the ranked compounds
    p <- molecules_ranked |>
        ggplot(aes(x = nTargets, y = max_llr, color = nAdverseEvents, label = name)) +
            geom_point(alpha = 0.6, shape = 19) +
            scale_y_continuous(transform = "pseudo_log", breaks = c(0, 10, 100, 1000)) +
            geom_text_repel(data = head(molecules_ranked, args$top), size = 3) +
            scale_color_viridis_c(name = "Number of \nAdverse Events", direction = 1) +
            labs(
                x = "Number of Targets (fewer = more selective)",
                y = "Maximum Log-Likelihood Ratio (lower = safer)",
                title = paste("Ranked drug candidates for target:", query_symbol, " (", query_ensembl, ")")
            ) +
            theme_minimal()

    ggsave(
        filename = file.path("results", paste0(res_name, ".png")),
        plot = p,
        dpi = 300
    )
    message("\nPlot saved to results/", res_name, ".png")
}
