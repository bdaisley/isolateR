
assign_taxonomy <- function(folder=NULL,
                            export_csv=TRUE,
                            verbose=TRUE,
                            exclude=NULL){
  
  # function requirements------------------------------------------------------------
#checking for required packages; installing those not yet installed
if(require(dplyr)==FALSE) install.packages('dplyr')
if(require(stringr)==FALSE) install.packages('stringr')
if(require(R.utils)==FALSE) install.packages('R.utils')
if(require(rentrez)==FALSE) install.packages('rentrez')
if(require(xmlconvert)==FALSE) install.packages('xmlconvert')
if(require(reactable)==FALSE) install.packages('reactable')
if(require(reactablefmtr)==FALSE) install.packages('reactablefmtr')
if(require(pander)==FALSE) install.packages('pander')
  
#loading required packages
library(dplyr)
library(stringr)
library(R.utils)
library(rentrez)
library(xmlconvert)
library(reactable)
library(reactablefmtr)
library(pander)

#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------

#folder <- "C:/Users/Brendan Daisley/Documents/GUELPH/00-Guelph-Projects/ARTICLE---2023---Novel_BEE_microbes/sanger_sequencing_results/2023_07_06_Dylan"

#-------------------------------------------------

path <- str_replace_all(folder, '\\\\', '/')
suppressWarnings(dir.create(file.path(path, 'ncbi_database')))
files  <- dir(file.path(path, 'ncbi_database'), full.names = TRUE)

#:::::::::::::::::::::::::::::::::
#Download NCBI 16S rRNA database
#:::::::::::::::::::::::::::::::::
if(!file.exists(file.path(path, 'ncbi_database/bacteria.16SrRNA.fna.gz'))){
  download.file("https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz", 
                file.path(path, 'ncbi_database/bacteria.16SrRNA.fna.gz'), mode='wb')
}

gz_files <- str_subset(files, 'usearch')


#:::::::::::::::::::::::::::
#Download USEARCH software
#:::::::::::::::::::::::::::
usearch_files <- str_subset(files, 'usearch')

if(identical(usearch_files, character(0))){
  if (grepl('w|W', .Platform$OS.type)) {
  ## we are on Windows
  download.file("https://drive5.com/downloads/usearch11.0.667_win32.exe.gz", file.path(path, 'ncbi_database/usearch11.0.667_win32.exe.gz'), mode='wb')
} else {
  if (grepl('darwin', version$os)) {
    ## Mac
    download.file("https://drive5.com/downloads/usearch11.0.667_i86osx32.gz", file.path(path, 'ncbi_database/usearch11.0.667_i86osx32.gz'), mode='wb')
  } else {
    ## Linux
    download.file("https://drive5.com/downloads/usearch11.0.667_i86linux32.gz", file.path(path, 'ncbi_database/usearch11.0.667_i86linux32.gz'), mode='wb')
  }}}

#::::::::::::
#Unzip files
#::::::::::::
gz_files <- str_subset(files, '.gz')
if(! identical(gz_files, character(0))){
unlink("ncbi_database/*.tmp", recursive = TRUE, force = FALSE)
lapply(gz_files, R.utils::gunzip, remove = FALSE, overwrite=TRUE)
}

#::::::::::::
#Stage paths
#::::::::::::
setwd(path)
query.fasta <- "V3-V6seq.FASTA"
db.fasta <- "ncbi_database/bacteria.16SrRNA.fna"
uc.out <- "V3-V6seq.uc"
b6.out <- "V3-V6seq.b6"
#db.fasta <- file.path(path, "ncbi_database/bacteria.16SrRNA.fna")
#uc.out <- file.path(path, "V3-V6seq_results.uc")
usearch_files <- str_subset(files, 'usearch')
usearch.path <- gsub(".gz", "", usearch_files[1], fixed=TRUE)


#:::::::::::::::::::::
#Search NCBI database
#:::::::::::::::::::::
message(cat(paste0("\n", "\033[97;", 40, "m","Searching query sequences against NCBI database via '--usearch_global' function of USEARCH", "\033[0m", "\n")))

system2(usearch.path, paste(" --usearch_global ", query.fasta, " --db ", db.fasta, " --output_no_hits --blast6out ", b6.out, " --uc ", uc.out, " --id 0.7 --maxaccepts 0 --maxrejects 0 --top_hits_only --strand both --threads 1 ", sep=""), stdout="", stderr="")

#:::::::::::::::::
#Organize results
#:::::::::::::::::
message(cat(paste0("\n", "\033[97;", 40, "m","Determining closest species match", "\033[0m", "\n")))
uc.results <- read.csv(uc.out, sep="\t", header = FALSE) %>% 
  mutate(V4 = gsub("*", "51", .$V4, fixed=TRUE)) %>% mutate_at(vars(V4), as.numeric) %>%   # Fix percent ID column to make numeric
  mutate(V10 = gsub("*", "No_match", .$V10, fixed=TRUE)) %>%                               # Fix unmatched taxa column
  arrange(desc(V4)) %>% distinct(V9, .keep_all=TRUE) #%>% select(V1, V4, V9, V10) #%>% filter(V10 != "*")

#-------------------------------------------------
query.seqs <- Biostrings::readBStringSet(query.fasta)
query.seqs.df <- as.data.frame(names(query.seqs)) %>% mutate(V9 = names(query.seqs), query_seq = paste(query.seqs)) %>% mutate(length= nchar(query_seq)) %>% mutate(Ns = str_count(.$query_seq, "[Nn]")) %>% select(-1)

merged.df <- merge(uc.results, query.seqs.df, by="V9", all=TRUE) %>% arrange(V9) %>%
  mutate(NCBI_acc = str_split_fixed(.$V10, " ", 8)[,1]) %>%
  mutate(species = paste(str_split_fixed(.$V10, " ", 8)[,2], str_split_fixed(.$V10, " ", 8)[,3], sep=" ")) %>%
  dplyr::rename("ID" = "V4", "sample" = "V9", "closest_match" = "V10") %>% 
  select(sample,ID,species,length,Ns,NCBI_acc,closest_match,query_seq)

#write.csv(merged.df, file=paste("V3-V6seq_results_", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))],".csv", sep=""))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------- library(rentrez)
message(cat(paste0("\n", "\033[97;", 40, "m","Collecting higher level taxonomic rank information of closest species match", "\033[0m", "\n", "\r\n")))

#::::::::::::::::::::::::::::::::
# Step 1b) For loop
#::::::::::::::::::::::::::::::::
id.list <- 1:length(merged.df$NCBI_acc)
id.chunk <- 250
fetch.ids <- split(id.list, ceiling(seq_along(id.list)/id.chunk))

fetch.list <- list()
for (i in 1:length(fetch.ids)) {
  #setTxtProgressBar(pb.bar,i)
  #t1 <- Sys.time()
  ii <- fetch.ids[[i]]
  while(TRUE){
    r_fetch.try <- try(entrez_fetch(db="nucleotide", id=merged.df$NCBI_acc[ii], rettype="gbc")) # OR use rettype="xml" OR useing rettype="gbc", silent=TRUE)
    if(is(r_fetch.try, 'try-error'))
      cat("Failed, trying again in 10 seconds...\n")
    Sys.sleep(3)
    if(!is(r_fetch.try, 'try-error')) break
  }
  r_fetch <- r_fetch.try
  #r_fetch <- entrez_fetch(db="nucleotide", id=merged.df$NCBI_acc[ii], rettype="gbc") # OR use rettype="xml" OR useing rettype="gbc"
  r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_update-date", "INSDSeq_update_date")
  r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_create-date", "INSDSeq_create_date")
  r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_primary.accession", "INSDSeq_primary_accession")
  r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_accession.version", "INSDSeq_accession_version")
  r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_other.seqids", "INSDSeq_other_seqids")
  r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_feature-table", "INSDSeq_feature_table")
  #Convert fetched XML file to dataframe format. In this case, the 'INSDSeq' node defines individual records within the XML tree
  #xml.df <- xml_to_df(text = r_fetch, records.tags = c("//INSDSeq"), check.datatypes=FALSE, check.names=TRUE)
  fetch.list[[paste("acc", i, sep="_")]] <- lapply(i, function(x) xmlconvert::xml_to_df(text = r_fetch, records.tags = c("//INSDSeq"), check.datatypes=FALSE, check.names=TRUE))
  #fetch.list[[paste("acc_1")]] <- rbind(as.data.frame(fetch.list[[1]]), as.data.frame(fetch.list[[paste("acc", i, sep="_")]]))
  #fetch.list[[paste("acc", i, sep="_")]] <- NULL
  svMisc::progress(i, length(fetch.ids))
  #svMisc::progress(((i/length(fetch.ids))*100), 100, progress.bar=TRUE)
  #Sys.sleep(0.01)
  #t2 <- Sys.time()
  #print(difftime(t2, t1, units = "secs")[[1]])
  #if (((i/length(fetch.ids))*100) == 100) message("Done!")
}

#fetch.list.df <- bind_rows(fetch.list, .id = "column_label") #Only works if more than 2 entries in the list
fetch.list.df <- as.data.frame(fetch.list$acc_1, stringsasfactors=FALSE) #Use this if only 1 entry in the list

#::::::::::::::::::::::::::::::::
lookup.list <- setNames(fetch.list.df$INSDSeq_taxonomy, fetch.list.df$INSDSeq_accession_version)
merged.df2 <- merged.df %>% mutate(taxonomy = dplyr::recode(NCBI_acc, !!!lookup.list)) %>%
  mutate(taxonomy = gsub(" ", "", .$taxonomy)) %>%
  mutate(domain = str_split_fixed(.$taxonomy, ";", 6)[,1]) %>%
  mutate(phylum = str_split_fixed(.$taxonomy, ";", 6)[,2]) %>%
  mutate(class = str_split_fixed(.$taxonomy, ";", 6)[,3]) %>%
  mutate(order = str_split_fixed(.$taxonomy, ";", 6)[,4]) %>%
  mutate(family = str_split_fixed(.$taxonomy, ";", 6)[,5]) %>%
  mutate(genus = str_split_fixed(.$taxonomy, ";", 6)[,6])

if(export_csv==TRUE) {
  write.csv(merged.df2, file=paste("assign_taxonomy_output___", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))],".csv", sep=""))
}

message(cat(paste0("\r\n\r\n", "\033[97;", 40, "m","Done.", "\033[0m", "\n")))

#::::::::::::::::::::::::
#Make reactable output
#::::::::::::::::::::::::

merged.df3 <- merged.df2 %>% select(sample,species,ID,Ns,length,query_seq,NCBI_acc,phylum,class,order,family,genus) %>%
                             mutate(species2 = species) %>% mutate(genus= str_split_fixed(genus, ";", 2)[,1]) %>%
                             mutate(species_colors = dplyr::case_when(genus != "" ~ "green", TRUE ~ "red")) %>% 
                             mutate(species_col = ifelse(ID > 98.9, "#31a354", NA)) %>%
                             mutate(genus_col = ifelse(ID > 95.7, "#31a354", NA)) %>%
                             mutate(family_col = ifelse(ID > 89.5, "#31a354", NA)) %>%
                             mutate(class_col = ifelse(ID > 86.1, "#31a354", NA)) %>%
                             mutate(order_col = ifelse(ID > 83.4, "#31a354", NA)) %>%
                             mutate(phylum_col = ifelse(ID > 80.7, "#31a354", NA))

#grDevices::adjustcolor("#31a354", alpha.f = 0.5) #97D1A9
#grDevices::adjustcolor(merged.df3$genus_col[[index]], alpha.f = 0.5)
#which(!is.na(merged.df3$species_col))
#merged.df3$species_colors

#:::::::::::::::::::::::::::::::
#List function for table below
#:::::::::::::::::::::::::::::::
list_with_names<-function(...){dplyr::lst(...)}
uni.font <- 11

#::::::::::::::::::::::::
#Reactable
#::::::::::::::::::::::::
merged.df3_react <- reactable(merged.df3,
                            fullWidth=FALSE, searchable = TRUE, bordered = TRUE, resizable =TRUE, #width = 1200,
                            defaultPageSize=100, highlight = TRUE, showSortable = TRUE, compact=TRUE, wrap = TRUE, #highlight=TRUE,
                            #theme = reactableTheme(headerStyle = list(wrap=TRUE)),
                            #theme = fivethirtyeight(centered = TRUE, header_font_size = 11),
                            theme = reactableTheme(                  #borderColor = "#555",
                              headerStyle = list(
                                fontFamily = "sans-serif", fontWeight="bold", fontSize=uni.font,
                               "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                               "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)"))),
                            #borderColor = "#555")),
                            #defaultColDef = colDef(minWidth = 100), #
                            #defaultColDef = colDef(style = reactablefmtr::cell_style(checkseq.sub, font_size=12)),
                            style = list(minWidth = 1400),
                            #--------------------Universal settings
                            defaultColDef = colDef(vAlign = "center", style = list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                            #--------------------Column specific settings
                            columns = list(
                                species2 = colDef(name="species", vAlign = "center", headerVAlign = "bottom", minWidth=170, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                               style = function(value, index, species2){
                                                 if(is.na(merged.df3$species_col[[index]])) { 
                                                   color.x <- "black" 
                                                 } else if(!is.na(merged.df3$species_col[[index]])) {
                                                   color.x <- "white" }
                                                 list(background = merged.df3$species_col[[index]], color = color.x, 
                                                      whiteSpace = "nowrap",
                                                      fontSize=uni.font,
                                                      #fontWeight="bold",
                                                      align="center",
                                                      fontFamily = "sans-serif")}),
                                genus = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=170, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                style = function(value, index, genus){
                                                  if(is.na(merged.df3$genus_col[[index]])) { 
                                                    color.x <- "black" 
                                                  } else if(!is.na(merged.df3$genus_col[[index]])) {
                                                    color.x <- "white" }
                                                  list(background = merged.df3$genus_col[[index]], color = color.x, 
                                                       whiteSpace = "nowrap",
                                                       fontSize=uni.font,
                                                       #fontWeight="bold",
                                                       align="center",
                                                       fontFamily = "sans-serif")}),
                                family = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=170, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                               style = function(value, index, family){
                                                 if(is.na(merged.df3$family_col[[index]])) { 
                                                   color.x <- "black" 
                                                 } else if(!is.na(merged.df3$family_col[[index]])) {
                                                   color.x <- "white" }
                                                  list(background = merged.df3$family_col[[index]], color = color.x, 
                                                       whiteSpace = "nowrap",
                                                       fontSize=uni.font, 
                                                       #fontWeight="bold",
                                                       align="center",
                                                       fontFamily = "sans-serif")}),
                                order = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=170, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                               style = function(value, index, order){
                                                 if(is.na(merged.df3$order_col[[index]])) { 
                                                   color.x <- "black" 
                                                 } else if(!is.na(merged.df3$order_col[[index]])) {
                                                   color.x <- "white" }
                                                 list(background = merged.df3$order_col[[index]], color = color.x, 
                                                      whiteSpace = "nowrap",
                                                      fontSize=uni.font, 
                                                      #fontWeight="bold",
                                                      align="center",
                                                      fontFamily = "sans-serif")}),
                                class = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=170, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                style = function(value, index, class){
                                                  if(is.na(merged.df3$class_col[[index]])) { 
                                                    color.x <- "black" 
                                                  } else if(!is.na(merged.df3$class_col[[index]])) {
                                                    color.x <- "white" }
                                                  list(background = merged.df3$class_col[[index]], color = color.x, 
                                                       whiteSpace = "nowrap",
                                                       fontSize=uni.font, 
                                                       #fontWeight="bold",
                                                       align="center",
                                                       fontFamily = "sans-serif")}),
                                phylum = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=170, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                               style = function(value, index, phylum){
                                                 if(is.na(merged.df3$phylum_col[[index]])) { 
                                                   color.x <- "black" 
                                                 } else if(!is.na(merged.df3$phylum_col[[index]])) {
                                                   color.x <- "white" }
                                                 list(background = merged.df3$phylum_col[[index]], color = color.x, 
                                                      whiteSpace = "nowrap",
                                                      fontSize=uni.font, 
                                                      #fontWeight="bold",
                                                      align="center",
                                                      fontFamily = "sans-serif")}),
                                species = colDef(name="Closest match", minWidth=170),
                                species_colors = colDef(show=FALSE),
                                species_col = colDef(show=FALSE),
                                genus_col = colDef(show=FALSE),
                                family_col = colDef(show=FALSE),
                                class_col = colDef(show=FALSE),
                                order_col = colDef(show=FALSE),
                                phylum_col = colDef(show=FALSE),
                                # Gauge widgets
                                #-------------------
                                ID = colDef(name = "ID", width= 60, 
                                            cell = gauge_chart(merged.df3, #fill_color_ref = "ID_log",
                                                   #fill_color = '#1A9641',
                                                   #fill_color = c('#FDAE61', '#FFFFBF','#A6D96A','#1A9641'), opacity = 0.5,
                                                   fill_color = c('#e5f5e0', '#a1d99b', '#31a354'), opacity = 0.5,
                                                   #viridis::viridis(5)[1:4]),
                                                   #'#D7191C','#FDAE61','#FFFFBF','#A6D96A','#1A9641'),
                                                   #number_fmt = scales::number_format(scale=50), #scales::log_trans(), #scales::comma,
                                                   #bold_text = TRUE,
                                                   text_size = 10,
                                                   #background = '#555555',
                                                   min_value = 51, max_value=100,
                                                   #max_value = max(c(checkseq.sub$phred_raw,checkseq.sub$phred_trim)),
                                                   show_min_max = FALSE)),
                                Ns = colDef(name = "Ns", width= 60, 
                                            cell = gauge_chart(merged.df3, #fill_color_ref = "ID_log",
                                                               #fill_color = '#1A9641',
                                                               #fill_color = c('#FDAE61', '#FFFFBF','#A6D96A','#1A9641'), opacity = 0.5,
                                                               fill_color = c("white", "#B20A2C"), opacity = 0.5,
                                                               #viridis::viridis(5)[1:4]),
                                                               #'#D7191C','#FDAE61','#FFFFBF','#A6D96A','#1A9641'),
                                                               #number_fmt = scales::number_format(scale=50), #scales::log_trans(), #scales::comma,
                                                               #bold_text = TRUE,
                                                               text_size = 10,
                                                               #background = '#555555',
                                                               min_value = 51, max_value=100,
                                                               #max_value = max(c(checkseq.sub$phred_raw,checkseq.sub$phred_trim)),
                                                               show_min_max = FALSE)),
                                # Barplot widgets
                                #-------------------
                                length = colDef(name = "trimLength", width=150, cell = data_bars(merged.df3,
                                                                                                    text_position = "inside-base",
                                                                                                    text_size = 10,
                                                                                                    #fill_color="#6f86ab",
                                                                                                    fill_color = c('#FFF2D9','#FFE1A6','#FFCB66','#FFB627'),
                                                                                                    fill_gradient = TRUE,
                                                                                                    background = 'transparent',
                                                                                                    number_fmt = scales::comma_format(accuracy = 0.1),
                                                                                                    round_edges = FALSE, align_bars="left"))
                             )) %>% add_title(paste("|assign_taxonomy| output for project: ", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))], sep=""))
                 
merged.df3_react             
fname_html <- paste0(path, "___assign_taxonomy_output.html", sep="")
suppressMessages(reactablefmtr::save_reactable_test(merged.df3_react, fname_html))
openFileInOS(fname_html)
}


#10^1.77
#}

#list(list(color_scales(merged.df3, color_ref = 'species_col', opacity = 0.5)), "white-space: nowrap", "font-size: 16px", "font-family: Calibri")
#lapply(merged.df3, function(x){
#       cell_style(merged.df3, font_size =10,  font_color = "white", rows = which(!is.na(merged.df3$species_col)), background_color = "#1673ba")
#       cell_style(merged.df3, font_size =10,  font_color = "white", rows = which(is.na(merged.df3$species_col)), background_color = "#1673ba")
#       })


#as.list(cell_style(merged.df3, font_size =10,  font_color = "white", rows = which(is.na(merged.df3$species_col)), background_color = "#1673ba"))

#style = cell_style(merged.df3, font_size =10,  font_color = "white", rows = which(is.na(merged.df3$species_col)), background_color = "#1673ba")),
#cell_style(merged.df3, font_size =10,  font_color = "white",
#            rows = which(is.na(merged.df3$species_col)), 
#          background_color = "grey"))),
#color_scales(merged.df3, color_ref = 'species_col', opacity = 0.5)),
#cell_style(merged.df3, font_size=20))), # "white-space: nowrap", "font-size: 16px", "font-family: Calibri")),
#style = as.list(color_scales(merged.df3, color_ref = 'species_col', opacity = 0.5), "white-space: nowrap", "font-size: 16px", "font-family: Calibri")),


#style = list(whiteSpace = "nowrap", fontFamily = "sans-serif")),
#style = list_with_names(color_scales(merged.df3, color_ref = 'species_colors', opacity = 0.5),
#                       whiteSpace = "nowrap", fontSize=10, fontFamily = "sans-serif")),
#style = list_with_names(color_scales(merged.df3, color_ref = 'species_colors', opacity = 0.5))),
