library(seqinr)
library(babette)
library(dplyr)

EV_fasta <- read.fasta("./sequence_alignment.fasta", as.string = T, forceDNAtolower = F)

names(EV_fasta)

EV_metadata <- read.table(file = 'metadata(1).tsv', sep = '\t', header = TRUE)

EV_metadata2 <- EV_metadata %>% 
  mutate(accession_2 = paste0(accession, ".1")) %>% 
  filter(names(EV_fasta) %in% accession_2) %>% 
  select(date, country, region, accession, accession_2, strain)

for(i in 1:nrow(EV_metadata2)){
  
  date_split <- NULL
  
  date_split <- strsplit(EV_metadata2$date[i], " ")
  
  if(length(date_split[[1]]) > 1){
    
    EV_metadata2$date[i] <- date_split[[1]][1]
    
  }else{
    EV_metadata2$date[i] <- date_split[[1]]
  }
  
}

US_metadata <- EV_metadata2 %>% 
  filter(country == "USA")

US_metadata_2 <- US_metadata[-unique(c(grep("EV-D68/Homosapiens", US_metadata$strain), grep("EVD")),]

EV_output <- bbt_run(EV_fasta, 
                     site_model = create_hky_site_model(), 
                     clock_model = create_strict_clock_model(), 
                     tree_prior = create_cbs_tree_prior(), 
                     mrca_prior = NA, 
                     mcmc = create_mcmc(),
                     beast2_input_filename = create_temp_input_filename(), 
                     rng_seed = 123, 
                     beast2_output_state_filename = create_temp_state_filename(), 
                     beast2_path = "/Applications/BEAST 2.6.6/bin/beast",
                     overwrite = T, 
                     verbose = T)
