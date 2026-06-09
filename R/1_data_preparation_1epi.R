library(tidyverse)

# load monocle data
mon <- read.csv("raw_data/monocle-metadata.csv") %>% 
  dplyr::transmute(
    id = Lane_id,
    
    # change continent, country & region to title
    Continent = stringr::str_to_title(stringr::str_to_lower(Continent)),
    Country = stringr::str_to_title(stringr::str_to_lower(Country)),
    Region = stringr::str_to_title(stringr::str_to_lower(Region)),
    Latitude = Latitude,
    Longitude = Longitude,
    Year = Year,
    Month = Month,
    
    # epiData
    Age_years = as.numeric(Age_years),
    Less_than_5_years_old = Less_than_5_years_old,
    Clinical_manifestation = Clinical_manifestation,
    Source = Source,
    Manifestation = Manifestation,
    
    # genomic
    Pipeline_version = Pipeline_version,
    In_silico_ST = In_silico_ST,
    GPSC = GPSC,
    In_silico_serotype = In_silico_serotype,
    
    # vaccines (23F is not included in PCV-21 alt)
    Vaccine_period = Vaccine_period,
    Introduction_year = Introduction_year,
    PCV_type = PCV_type
    
  ) %>% 
  glimpse()


# load MSD data (all GPSC14), adapted from LSMS phylogenomic script
df_epi_gen_pneumo <- read.csv("raw_data/genData_pneumo_with_epiData_with_final_pneumo_decision.csv") %>% 
  dplyr::right_join(
    read.table("raw_data/qfile_filtered_19to23.txt") %>% 
      dplyr::mutate(specimen_id = V1,
                    workPoppunk_qc = "pass_qc") %>% 
      dplyr::select(specimen_id, workPoppunk_qc)
    ,
    by = "specimen_id"
    
  ) %>% 
  dplyr::filter(workPoppunk_qc == "pass_qc") %>%
  dplyr::rename(label = specimen_id) %>%
  dplyr::mutate(serotype_final_decision = factor(serotype_final_decision,
                                                 levels = c(
                                                   # VT
                                                   # "3", "6A/6B", "6A/6B/6C/6D", "serogroup 6",
                                                   # "14", "17F", "18C/18B", "19A", "19F", "23F",
                                                   "1", "3", "4", "5", "7F",
                                                   "6A", "6B", "9V", "14", "18C",
                                                   "19A", "19F", "23F",
                                                   # NVT
                                                   # "7C", "10A", "10B", "11A/11D", "13", "15A", "15B/15C",
                                                   # "16F", "19B", "20", "23A", "23B", "24F/24B", "25F/25A/38",
                                                   # "28F/28A", "31", "34", "35A/35C/42", "35B/35D", "35F",
                                                   # "37", "39", "mixed serogroups",
                                                   "serogroup 6", "6C", "7C",
                                                   "10", "10A", "10B", "11A", "13",
                                                   "15A", "15B", "15C","15B/15C", "16F",
                                                   "17F", "18A", "18B", "19B", "20", "20B",
                                                   "21", "23A", "23B", "23B1",
                                                   "24F", "24B/C/F", "24B/24C/24F", "serogroup 24",
                                                   "25B", "25F",
                                                   "28A", "31", "33B", "33G",
                                                   "34", "35A", "35B", "35C", "35F", "37",
                                                   "37F", "38", "39",
                                                   "NT")),
                # test individual serotype
                # serotype_final_decision = ifelse(serotype_final_decision == "NT", "NT", "others"),
                
                workWGS_gpsc_strain = ifelse(workWGS_gpsc_strain == "Not assigned", "not assigned",
                                             workWGS_gpsc_strain),
                vaccination_status_area = case_when(
                  area == "Lombok" |
                    area == "Sumbawa" ~ "PCV13-implemented area (Lombok & Sumbawa)",
                  TRUE ~ "Pre-implemented area (Minahasa & Sorong)"
                ),
                
                # reset logic label for AMR-MDR viz
                workWGS_AMR_logic_class_chloramphenicol = case_when(
                  workWGS_AMR_chloramphenicol != "NF" ~ stringr::str_extract(workWGS_AMR_chloramphenicol,
                                                                             "(?<=R \\().*(?=\\))"),
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_chloramphenicol = case_when(
                  workWGS_AMR_logic_class_chloramphenicol == "cat_pC194" ~ "cat (pC194)",
                  workWGS_AMR_logic_class_chloramphenicol == "cat_q" ~ "catQ",
                  TRUE ~ workWGS_AMR_logic_class_chloramphenicol
                ),
                workWGS_AMR_logic_class_chloramphenicol = factor(workWGS_AMR_logic_class_chloramphenicol),
                
                workWGS_AMR_logic_class_clindamycin = case_when(
                  workWGS_AMR_clindamycin != "NF" ~ stringr::str_extract(workWGS_AMR_clindamycin,
                                                                         "(?<=R \\().*(?=\\))"),
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_clindamycin = factor(workWGS_AMR_logic_class_clindamycin),
                
                workWGS_AMR_logic_class_erythromycin = case_when(
                  workWGS_AMR_erythromycin != "NF" ~ stringr::str_extract(workWGS_AMR_erythromycin,
                                                                          "(?<=R \\().*(?=\\))"
                  ),
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_erythromycin = case_when(
                  workWGS_AMR_logic_class_erythromycin == "mefA_10" ~ "mefA",
                  TRUE ~ workWGS_AMR_logic_class_erythromycin
                ),
                workWGS_AMR_logic_class_erythromycin = factor(workWGS_AMR_logic_class_erythromycin),
                
                workWGS_AMR_logic_class_fluoroquinolones = case_when(
                  workWGS_AMR_fluoroquinolones != "NF" ~ stringr::str_extract(workWGS_AMR_fluoroquinolones,
                                                                              "(?<=R \\().*(?=\\))"
                  ),
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_fluoroquinolones = case_when(
                  workWGS_AMR_logic_class_fluoroquinolones == "parC_S79Y" |
                    workWGS_AMR_logic_class_fluoroquinolones == "parC_D83N" ~ "parC",
                  TRUE ~ workWGS_AMR_logic_class_fluoroquinolones
                ),
                workWGS_AMR_logic_class_fluoroquinolones = factor(workWGS_AMR_logic_class_fluoroquinolones),
                
                workWGS_AMR_logic_class_tetracycline = case_when(
                  workWGS_AMR_tetracycline != "NF" ~ stringr::str_extract(workWGS_AMR_tetracycline,
                                                                          "(?<=R \\().*(?=\\))"
                  ),
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_tetracycline = case_when(
                  stringr::str_detect(workWGS_AMR_logic_class_tetracycline, "tetM") ~ "tet(M)",
                  stringr::str_detect(workWGS_AMR_logic_class_tetracycline, "tet_") ~ "tet",
                  workWGS_AMR_logic_class_tetracycline == "tetK" ~ "tet(K)",
                  TRUE ~ workWGS_AMR_logic_class_tetracycline
                ),
                workWGS_AMR_logic_class_tetracycline = factor(workWGS_AMR_logic_class_tetracycline),
                
                workWGS_AMR_logic_class_antifolates = case_when(
                  workWGS_AMR_class_antifolates != "NF" ~ stringr::str_extract(workWGS_AMR_class_antifolates,
                                                                               "(?<=R \\().*(?=\\))"
                  ),
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_antifolates = case_when(
                  workWGS_AMR_logic_class_antifolates == "folA_I100L" ~ "folA",
                  workWGS_AMR_logic_class_antifolates == "folP_57-70" ~ "folP",
                  workWGS_AMR_logic_class_antifolates == "folA_I100L & folP_57-70" ~ "folA & folP",
                  TRUE ~ workWGS_AMR_logic_class_antifolates
                ),
                workWGS_AMR_logic_class_antifolates = factor(workWGS_AMR_logic_class_antifolates),
                
                # SIR format
                workWGS_AMR_logic_class_cephalosporins = case_when(
                  workWGS_AMR_logic_class_cephalosporins ~ " Resistance",
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_cephalosporins = factor(workWGS_AMR_logic_class_cephalosporins, levels = c(" Not found", " Resistance")),
                
                workWGS_AMR_logic_class_penicillins = case_when(
                  workWGS_AMR_logic_class_penicillins ~ " Resistance",
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_penicillins = factor(workWGS_AMR_logic_class_penicillins, levels = c(" Not found", " Resistance")),
                
                workWGS_AMR_logic_class_meropenem = case_when(
                  workWGS_AMR_logic_class_meropenem ~ " Resistance",
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_meropenem = factor(workWGS_AMR_logic_class_meropenem, levels = c(" Not found", " Resistance")),
                
                workWGS_AMR_MDR_flag = ifelse(workWGS_AMR_MDR_flag == "non-MDR", " Not found", " MDR"),
                workWGS_AMR_MDR_flag = factor(workWGS_AMR_MDR_flag,
                                              levels = c(" Not found", " MDR")),
                
  ) %>% 
  glimpse()


# analyse NT from Salma (including setting up NT priority from 636 subtree) ####
nt_gps <- readxl::read_excel("raw_data/Result_NT_GPS.xlsx") %>% 
  dplyr::full_join(
    list.files("raw_data/subtree_636_nontypeables/",
               pattern = "\\.fasta$",
               full.names = FALSE
    ) |>
      tibble::tibble(file = _) %>% 
      dplyr::mutate(file = gsub(".fasta", "", file),
                    priority = "priority 17 NT")
    ,
    by = c("Sample_ID" = "file")
  ) %>% 
  dplyr::full_join(
    read.delim("raw_data/636_nontypeables_poppunk_analysedbyHarry.txt",
               header = TRUE, sep = ",") %>% 
      dplyr::mutate(analysed = "popPUNK by Harry")
    ,
    by = c("Sample_ID" = "sample")
  ) %>% 
  dplyr::mutate(
    # update analysed col
    analysed = case_when(
      !is.na(Overall_QC) ~ "gps-pipeline by Salma",
      TRUE ~ analysed
    ),
    compiled_gpsc = GPSC.x,
    compiled_gpsc = ifelse(is.na(compiled_gpsc), GPSC.y, compiled_gpsc),
    done_gpsPipeline = ifelse(is.na(Overall_QC), 0, 1),
    done_popPUNK = 1
  ) %>% 
  dplyr::filter(priority == "priority 17 NT") %>% 
  # temporary select
  dplyr::select(1:7, GPSC.x, Serotype, ST, GPSC.y, compiled_gpsc,
                done_gpsPipeline, done_popPUNK, analysed, priority) %>%
  glimpse()

files <- list.files("raw_data/subtree_636_nontypeables/",
                    pattern = "\\.fasta$",
                    full.names = FALSE
) |>
  tibble::tibble(file = _) %>% 
  dplyr::mutate(file = gsub(".fasta", "", file)) %>% 
  glimpse()
