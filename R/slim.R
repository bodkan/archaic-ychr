# header names of the full SLiM output file
SECTION_HEADERS <- c("Populations", "Mutations", "Individuals", "Genomes")


#' Read SLiM output file generated by the sim.outputFull() method.
#'
#' @param path Path to a SLiM output file.
#'
#' @export
read_slim <- function(path) {
  readLines(path)
}


#' Parse the population section of the SLiM file into a data frame.
#'
#' @param slim_file Lines of the whole SLiM output file.
#'
#' @export
read_populations <- function(slim_file) {
  read_section_data(slim_file, "Populations") %>%
    paste0(collapse = "\n") %>%
    readr::read_delim(delim = " ", col_names = c("pop_id", "pop_size", "sex", "ratio")) %>%
    tibble::as_tibble()
}


#' Parse the mutation section of the SLiM file into a data frame.
#'
#' @param slim_file Lines of the whole SLiM output file.
#' @param mut_type,pop,time Mutation type, population of origin
#'   and generation of origin of mutations to filter from the file.
#'
#' @export
read_mutations <- function(slim_file, mut_type = "", pop = "", time = 0) {
  # read the mutation section lines and filter those mutations
  # that match the given pop. origin and mutation type
  read_section_data(slim_file, "Mutations") %>%
    stringr::str_subset(stringr::str_c(pop, " \\d+ \\d+$")) %>%
    stringr::str_subset(stringr::str_c("^\\d+ \\d+ ", mut_type)) %>%
    paste0(collapse = "\n") %>%
    readr::read_delim(
      delim = " ",
      col_names=c("mut_id", "run_id", "mut_type",
                  "pos", "s", "h", "pop_origin",
                  "gen_origin", "count"),
      progress = FALSE
    ) %>%
    dplyr::filter(gen_origin >= time)
}


#' Parse the genome section of the SLiM file into a data frame.
#'
#' @param slim_file Lines of the whole SLiM output file.
#' @param pop Population identifier.
#'
#' @export
read_genomes <- function(slim_file, pop = "") {
  # read the genomes section and filter for genomes belonging
  # to the given population
  lines <-
    read_section_data(slim_file, "Genomes") %>%
    stringr::str_subset(stringr::str_c("^", pop))
  # split each row into three parts - genome ID, auto/sex chromosome
  # and a string with mutation IDs and convert this to a 3 column tibble
  lines %>%
    stringr::str_split(" ", 3, simplify=TRUE) %>%
    tibble::as_tibble() %>%  # convert matrix into tibble
      dplyr::select(genome_id=V1, chrom_type=V2, mutations=V3) %>%
      dplyr::mutate(mutations = ifelse(mutations == "<null>", NA, mutations)) %>%
      dplyr::mutate(mutations=stringr::str_split(mutations, " ")) %>%
      tidyr::unnest(cols = c(mutations)) %>%
      dplyr::mutate(mut_id=as.integer(mutations)) %>%
      dplyr::select(-mutations)
}


#' Parse the individual section of the SLiM file into a data frame.
#'
#' @param slim_file Lines of the whole SLiM output file.
#' @param pop Population identifier.
#'
#' @export
read_individuals <- function(slim_file, pop = "") {
  read_section_data(slim_file, "Individuals") %>%
    paste0(collapse = "\n") %>%
    readr::read_delim(delim = " ",
                      col_names=c("indiv_id", "sex", "genome1_id","genome2_id"),
                      col_types = "cccc") %>%
  tibble::as_tibble() %>%
  dplyr::filter(stringr::str_detect(indiv_id, stringr::str_c("^", pop)))
}


# Read the positions of the SLiM output file section headers
# used to break up the whole file into individual chunks
# (population description, mutations, individuals and genomes).
read_section_delims <- function(slim_file) {
  # get positions of the section headers
  section_pos <- c(which(slim_file %in% stringr::str_c(SECTION_HEADERS, ":")),
                   length(slim_file) + 1)

  # put information about starts/ends of individuals section
  # into a data frame
  tibble::tibble(section=SECTION_HEADERS,
                 start=section_pos[-length(section_pos)] + 1,
                 end=section_pos[-1] - 1)
}


# Get the subset of SLiM output file contents for a given section.
read_section_data <- function(slim_file, section_name) {
  if (! section_name %in% SECTION_HEADERS)
    stop("Invalid SLiM output section specified!")

  delims <- read_section_delims(slim_file)
  pos_in_file <- dplyr::filter(delims, section == section_name)

  slim_file[pos_in_file$start:pos_in_file$end]
}


# Convert SLiM individual ID into genome IDs.
get_genome_ids <- function(indiv_id) {
  # split string such as "p3:i0" into "p3:" and "0"
  split_id <- stringr::str_split(indiv_id, "i", simplify=TRUE)

  # save both parts to separate variables
  pop_id <- split_id[, 1]
  indiv_no <- split_id[, 2] %>% as.integer

  # use those variables to construct
  stringr::str_c(pop_id, c(2 * indiv_no, 2 * indiv_no + 1))
}

#' Calculate genetic fitness for each simulated Y chromosome.
calculate_fitness <- function(f, pop = "") {
    muts <- read_mutations(f)
    genomes <- read_genomes(f, pop = pop) %>%
        dplyr::mutate(pop = str_extract(genome_id, "p[0-9]+")) %>%
        dplyr::filter(!is.na(mut_id))

    combined <- dplyr::inner_join(genomes, muts, by = "mut_id") %>%
        dplyr::select(genome_id, pop, mut_id, mut_type, s, pop_origin, gen_origin) %>%
        dplyr::mutate(s = -s)

    fitness <- dplyr::group_by(combined, genome_id, pop) %>%
        dplyr::summarise(S = sum(s)) %>%
        dplyr::mutate(fitness = exp(-S))

    fitness
}
