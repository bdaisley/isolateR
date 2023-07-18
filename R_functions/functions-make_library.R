
make_library <- function(new_lib_csv=NULL,
                         old_lib_csv=NULL,
                         include_warnings=FALSE,
                         verbose=TRUE){
  
  # function requirements------------------------------------------------------------
  #checking for required packages; installing those not yet installed
  if(require(dplyr)==FALSE) install.packages('dplyr')
  if(require(stringr)==FALSE) install.packages('stringr')
  if(require(Biostrings)==FALSE) install.packages('Biostrings')
  if(require(htmltools)==FALSE) install.packages('htmltools')
  if(require(reactable)==FALSE) install.packages('reactable')
  if(require(reactablefmtr)==FALSE) install.packages('reactablefmtr')
  if(require(pander)==FALSE) install.packages('pander')
  if(require(crosstalk)==FALSE) install.packages('crosstalk')

  #loading required packages
  library(dplyr)
  library(stringr)
  library(Biostrings)
  library(htmltools)
  library(reactable)
  library(reactablefmtr)
  library(pander)
  library(crosstalk)

  
#------------------------------------------------- functions
    #Source "make_fasta" command to allow for CSV -> FASTA format conversion
    #make_fasta(csv_file=NULL,col_names="ID",col_seqs="Sequence")
    source("https://raw.githubusercontent.com/bdaisley/sangerseq2taxonomy/main/R_functions/make_fasta.R")
  
  #-------------------------------------------------
  #-------------------------------------------------
  #-------------------------------------------------
  #-------------------------------------------------
  #-------------------------------------------------

new_lib <- str_replace_all(new_lib_csv, '\\\\', '/')
new_lib_path <- paste(unlist(strsplit(new_lib, '/'))[1:(length(unlist(strsplit(new_lib, '/')))-1)],collapse="/")
new_lib_folder <- unlist(strsplit(new_lib_path, '/'))[length(unlist(strsplit(new_lib_path, '/')))]
new_lib_folder2 <- unlist(strsplit(new_lib_path, '/'))[(length(unlist(strsplit(new_lib_path, '/')))-1)]
new_lib_file <- read.csv(new_lib, row.names = 1)

setwd(new_lib_path)

#-------------------------------------------------
#::::::::::::
#Stage paths
#::::::::::::
message(cat(paste0("\n", "\033[0;", 32, "m","Staging files for de-replication", "\033[0m", "\n")))

if(include_warnings == TRUE){
  write.csv(new_lib_file, file=paste("make_library_output___", new_lib_folder2, ".csv", sep=""))
  new_lib_file <- read.csv(paste("make_library_output___", new_lib_folder2, ".csv", sep=""), row.names = 1)
  new_lib_file_path <- paste("make_library_output___", new_lib_folder2, ".csv", sep="")
} else {
  new_lib_file <- new_lib_file %>% filter(warning=="")
  write.csv(new_lib_file, file=paste("make_library_output___", new_lib_folder2, ".csv", sep=""))
  new_lib_file <- read.csv(paste("make_library_output___", new_lib_folder2, ".csv", sep=""), row.names = 1)
  new_lib_file_path <- paste("make_library_output___", new_lib_folder2, ".csv", sep="")
  }

make_fasta(new_lib_file_path, col_names="filename", col_seqs="query_seq") # exports as "output.fasta"
output.fasta <- "output.fasta"
out.fasta <- paste0("make_library_output___", new_lib_folder2, ".fasta", sep="")
out.drep <- paste0("make_library_output___", new_lib_folder2, ".drep", sep="")
vsearch.path <- file.path(new_lib_path,"NCBI_databases/vsearch-2.23.0.exe")

#-------------------------------------------------

#system2(vsearch.path, paste(" --derep_fulllength ", "output.fasta", " --output ", out.fasta, " --uc ", out.drep1, " --strand both --sizein --sizeout", sep=""), stdout="", stderr="")
system2(vsearch.path, paste(" --derep_prefix ", "output.fasta", " --output ", out.fasta, " --uc ", out.drep, " --sizein --sizeout", sep=""), stdout="", stderr="")
file.remove("output.fasta")
#system2(vsearch.path, paste(" --derep_fulllength ", "output.fasta", " --output ", out.fasta, " --uc ", out.uc1, sep=""), stdout="", stderr="")

#-------------------------------------------------
#:::::::::::::::::
#Organize results
#:::::::::::::::::
drep.results <- read.csv(out.drep, sep="\t", header = FALSE) %>% 
  filter(V1 != "C") %>% 
  mutate(V10 = ifelse(V10=="*", V9, V10)) %>% select(V9,V10) %>% #V10 is grouping column
  dplyr::rename("filename" = "V9", "grouping" = "V10")

#:::::::::::::::::::::::::::::::::::::::::::::
#Merge de-rep results with *new* lib file
#:::::::::::::::::::::::::::::::::::::::::::::

merged.drep <- merge(drep.results, new_lib_file, by.x="filename", by.y="filename", all=TRUE) %>%
  mutate(grouping = ifelse(is.na(grouping) == TRUE, filename, grouping))

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Pick the best representative for de-replicated sequence list
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

merged.drep1 <- merged.drep %>% select(-warning) %>% #select(-phylum,-class,-order,-family,-class,-genus, -species_col, -genus_col, -family_col) %>%
  mutate(qual_bin = ifelse(phred_trim >= 60, 6, 
                    ifelse(phred_trim >= 50 & phred_trim <60, 5,
                    ifelse(phred_trim >=40 & phred_trim <50, 4,
                    ifelse(phred_trim >=30 & phred_trim <40, 3,
                    ifelse(phred_trim >= 20 & phred_trim <30, 2, 1)))))) %>%
  mutate(override_priority = 1) %>%
  arrange(desc(qual_bin)) %>% 
  arrange(desc(length)) %>% #filter(grepl("1015|HB534", filename))
  arrange(desc(override_priority)) %>%
  group_by(grouping) %>% 
  mutate(strain_group = dplyr::first(filename)) %>%
  mutate(species = dplyr::first(species)) %>% 
  mutate(species2 = dplyr::first(species2)) %>% 
  #mutate(species2 = first(species2)) %>% 
    ungroup() %>%
  mutate(ref_strain = ifelse(grouping == filename, "yes", "no"))

#Subset only columns of interest------------------------------------------------------------------

merged.drep1.sub <- merged.drep1 %>% select(strain_group, date, filename,ID, species2, NCBI_acc, 
                                            phred_trim, Ns, length, query_seq,
                                            phylum, class, order, family, genus, species, ref_strain,
                                            phylum_col, class_col, order_col, family_col, genus_col, species_col)

#Combining old database if provided---------------------------------------------------------------


if(is.null(old_lib_csv)==FALSE){
  old_lib_file <- read.csv(old_lib_csv, row.names=1)
  message(cat(paste0("\n", "\033[0;", 32, "m","Checking old/new library CSV files have same dimensions", "\033[0m")))
  if(!ncol(old_lib_file)==ncol(merged.drep1.sub)) stop('Library file dimensions are not the same and require manual inspection.', call.=FALSE)
  message(cat(paste0("\033[0;", 32, "m","Same dimensions detected, proceeding...", "\033[0m", "\n")))
  #
  combined_lib_file <- rbind(old_lib_file, merged.drep1.sub) #old_lib_file must be at top so older seqs take priority
  #make unique strain names incase any duplicated-----------------------
  combined_lib_file$unique_filename <- paste(str_pad(1:nrow(combined_lib_file), 4, pad="0"), combined_lib_file$filename, sep="_")
  write.csv(combined_lib_file, file=paste("make_library_output_temp___", new_lib_folder2, ".csv", sep=""))
  combined_lib_file_path <- paste("make_library_output_temp___", new_lib_folder2, ".csv", sep="")

  #De-replicating-------------------------------------------------------
  make_fasta(combined_lib_file_path, col_names="unique_filename", col_seqs="query_seq") # exports as "output.fasta"
  output.fasta <- "output.fasta"
  out.fasta <- paste0("make_library_output_temp___", new_lib_folder2, ".fasta", sep="")
  out.drep <- paste0("make_library_output_temp___", new_lib_folder2, ".drep", sep="")
  vsearch.path <- file.path(new_lib_path,"NCBI_databases/vsearch-2.23.0.exe")
  
  #-------------------------------------------------
  
  #system2(vsearch.path, paste(" --derep_fulllength ", "output.fasta", " --output ", out.fasta, " --uc ", out.drep1, " --strand both --sizein --sizeout", sep=""), stdout="", stderr="")
  system2(vsearch.path, paste(" --derep_prefix ", "output.fasta", " --output ", out.fasta, " --uc ", out.drep, " --sizein --sizeout", sep=""), stdout="", stderr="")
  file.remove("output.fasta")
  #system2(vsearch.path, paste(" --derep_fulllength ", "output.fasta", " --output ", out.fasta, " --uc ", out.uc1, sep=""), stdout="", stderr="")
  #-------------------------------------------------
  #:::::::::::::::::
  #Organize results
  #:::::::::::::::::
  drep.results <- read.csv(out.drep, sep="\t", header = FALSE) %>% 
    filter(V1 != "C") %>% 
    mutate(V10 = ifelse(V10=="*", V9, V10)) %>% select(V9,V10) %>% #V10 is strain_group column
    dplyr::rename("unique_filename" = "V9", "grouping" = "V10")

  #:::::::::::::::::::::::::::::::::::::::::::::
  #Merge de-rep results with *combined* lib file
  #:::::::::::::::::::::::::::::::::::::::::::::
  regroup.list.list1 <- setNames(drep.results$grouping, drep.results$unique_filename)

  combined_lib_file_regroup <- combined_lib_file %>%
    mutate(grouping = dplyr::recode(unique_filename, !!!regroup.list.list1))
    
  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  #Pick the best representative for de-replicated sequence list
  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  merged.drep2 <- combined_lib_file_regroup %>% #%>% select(-warning) %>% #select(-phylum,-class,-order,-family,-class,-genus, -species_col, -genus_col, -family_col) %>%
    arrange(unique_filename) %>%
    group_by(grouping) %>% 
    mutate(strain_group = dplyr::first(unique_filename)) %>%
    mutate(species = dplyr::first(species)) %>% 
    mutate(species2 = dplyr::first(species2)) %>% 
    #mutate(species2 = first(species2)) %>% 
    ungroup() %>%
    mutate(ref_strain = ifelse(strain_group == unique_filename, "yes", "no")) %>%
    #Remove unique identifier characters
    mutate(strain_group = sapply(1:length(.$strain_group),function(x) str_sub(.$strain_group[x], 6, nchar(.$strain_group[x]))))


  merged.drep1.sub <- merged.drep2 %>% select(strain_group, date, filename,ID, species2, NCBI_acc, 
                                              phred_trim, Ns, length, query_seq,
                                              phylum, class, order, family, genus, species, ref_strain,
                                              phylum_col, class_col, order_col, family_col, genus_col, species_col)
  } #end of old library add-in


#merged.drep1 %>% group_by(grouping) %>% mutate(strain_group = last(filename))
#merged.drep1.sub <- merged.drep1 %>% select(grouping,strain_group) %>% mutate(same = case_when(grouping == strain_group ~ "same", TRUE ~ "other"))
#merged.drep1 <- merged.drep


###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
#:::::::::::::::::::::::::::::::
#List function for table below
#:::::::::::::::::::::::::::::::
list_with_names<-function(...){dplyr::lst(...)}
uni.font <- 8

#::::::::::::::::::::::::
#Reactable
#::::::::::::::::::::::::

library(crosstalk)

data <- SharedData$new(merged.drep1.sub)

new_library_reactable <- bscols(widths = c(1,10),
                                list(
                                  #filter_checkbox("categories", "Categories", data, ~categories, inline = TRUE),
                                  filter_slider("length", "Seq length", round=-1, ticks=FALSE, data, ~length),
                                  filter_slider("ID", "% identity", data, round=2, ticks=FALSE, ~ID),
                                  filter_checkbox("ref_strain", "Ref strain", data, ~ref_strain, inline = FALSE),
                                  filter_checkbox("date", "Date", data, ~date, inline = FALSE)
                                ),
                                htmltools::browsable(
                                  tagList(tags$button(
                                            "Expand/collapse all",
                                            onclick = "Reactable.toggleAllRowsExpanded('merged.drep1.sub')"
                                          ),
                                          tags$button("Download as CSV", onclick = "Reactable.downloadDataCSV('merged.drep1.sub')"),
                                          #new_library_reactable <- 
                                          reactable(data,
                                                    fullWidth=FALSE, searchable = TRUE, bordered = TRUE, resizable =TRUE, #width = 1200,
                                                    defaultPageSize=200, highlight = TRUE, showSortable = TRUE, compact=TRUE, #wrap = TRUE,
                                                    style = list(minWidth = 1400),
                                                    theme = reactableTheme(                  #borderColor = "#555",
                                                      headerStyle = list(
                                                        fontFamily = "sans-serif", fontWeight="bold", fontSize=uni.font,
                                                        "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                                                        "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)"))),
                                                    defaultColDef = colDef(minWidth= 75, vAlign = "center", style = list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                                                    groupBy = "strain_group",
                                                    elementId = "merged.drep1.sub",
                                              #colDef options
                                              #---------------
                                              columns = list(
                                              strain_group = colDef(name="strain_group", minWidth=130,
                                                                            grouped = JS("function(cellInfo) {
                                                                      if (cellInfo.subRows.length > 100) {
                                                                        return cellInfo.value + ' (' + cellInfo.subRows.length + ')'
                                                                      }
                                                                      return cellInfo.value
                                                                    }")),
                                                                    #V10 = colDef(style = "white-space: nowrap;", aggregate = "unique"),
                                              filename = colDef(name = "duplicate_strains", minWidth=120, aggregate = "count", style = list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                                              #species2 = colDef(minWidth=200, aggregate = "unique"),
                                              species2 = colDef(name="closest_match", minWidth=125, #aggregate = "unique"
                                                                aggregate = JS("function(values, rows){
                                                             return values[0]
                                                            }")),
                                              date = colDef(minWidth= 58, style = list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif"),
                                                            aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                              #query_seq = colDef(minWidth= 150, style = list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                                              ID = colDef(minWidth=37, aggregate = "max", format = colFormat(digits = 2)),
                                              #length = colDef(#aggregate = "mean", format = colFormat(digits = 0)),
                                              #              minWidth=60,
                                              #              aggregate = JS("function(values, rows){
                                              #                                             return values[0]
                                              #                        }"),
                                              #              cell = data_bars(merged.drep1.sub)),
                                              length = colDef(name = "trimLength", width=90,
                                                              aggregate = JS("function(values, rows){
                                                                                           return values[0]
                                                                      }"),
                                                              cell = data_bars(merged.drep1.sub,
                                                                     text_position = "inside-base",
                                                                     text_size = 10,
                                                                     fill_color = c('#FFF2D9','#FFE1A6','#FFCB66','#FFB627'),
                                                                     fill_gradient = TRUE,
                                                                     background = 'transparent',
                                                                     number_fmt = scales::comma_format(accuracy = 0.1),
                                                                     round_edges = FALSE, align_bars="left")),
                                              NCBI_acc = colDef(width=58,
                                                                aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                              phred_trim = colDef(name = "Q", width= 35,
                                                                  aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                              Ns = colDef(name = "Ns", width= 40,
                                                          aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                              query_seq = colDef(name="seq", width= 55, style = list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif"),
                                                                 aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                              phylum = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=90, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                              #cell = list(margin = list( top= 2, right= 5, bottom= 2, left= 5 )),
                                                              aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(merged.drep1.sub, color_ref="phylum_col")), #color_ref="phylum_col"
                                              class = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=90, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                              #cell = list(margin = list( top= 2, right= 5, bottom= 2, left= 5 )),
                                                              aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(merged.drep1.sub, color_ref="class_col")), 
                                              order = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=90, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                             #cell = list(margin = list( top= 2, right= 5, bottom= 2, left= 5 )),
                                                             aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                             cell = color_tiles(merged.drep1.sub, color_ref="order_col")), 
                                              family = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=100, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                             #cell = list(margin = list( top= 2, right= 5, bottom= 2, left= 5 )),
                                                             aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                             cell = color_tiles(merged.drep1.sub, color_ref="family_col")), 
                                              genus = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=100, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                             #cell = list(margin = list( top= 2, right= 5, bottom= 2, left= 5 )),
                                                             aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                             cell = color_tiles(merged.drep1.sub, color_ref="genus_col")), 
                                              species = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=135, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                             #cell = list(margin = list( top= 2, right= 5, bottom= 2, left= 5 )),
                                                             aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                             cell = color_tiles(merged.drep1.sub, color_ref="species_col")), 
                                              #species = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=170, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                              #                 style = function(value, index, species){
                                              #                   if(is.na(merged.drep1.sub$species_col[[index]])) { 
                                              #                     color.x <- "black" 
                                              #                   } else if(!is.na(merged.drep1.sub$species_col[[index]])) {
                                              #                     color.x <- "white" }
                                              #                   list(background = merged.drep1.sub$species_col[[index]], color = color.x, 
                                              #                        whiteSpace = "nowrap",
                                              #                        fontSize=uni.font, 
                                              #                        #fontWeight="bold",
                                              #                        align="center",
                                              #                        fontFamily = "sans-serif")}),
                                              #grouping = colDef(show=FALSE),
                                              ref_strain = colDef(width=57),
                                              species_col = colDef(show=FALSE),
                                              genus_col = colDef(show=FALSE),
                                              family_col = colDef(show=FALSE),
                                              class_col = colDef(show=FALSE),
                                              order_col = colDef(show=FALSE),
                                              phylum_col = colDef(show=FALSE)
                                                    )) %>% #end of reactable
                                            #add_subtitle("#", font_size = 24, font_style="normal", font_weight="normal") %>% 
                                            add_subtitle("Output table for 'make_library' command", font_size = 24, font_style="normal", font_weight="bold") %>% 
                                            add_subtitle(paste("Date generated: ", Sys.Date(), sep=""), font_size = 14, font_style="normal", font_weight="normal"), 
                                  )))

htmltools::save_html(new_library_reactable, paste("make_library_output___", new_lib_folder2, ".html", sep=""))
openFileInOS(paste("make_library_output___", new_lib_folder2, ".html", sep=""))
write.csv(merged.drep1.sub, file=paste("make_library_output___", new_lib_folder2, ".csv", sep=""))

message(cat(paste0("\033[97;", 40, "m","Export directory:", "\033[0m",
                   "\033[0;", 32, "m", " ", file.path(new_lib_path), "\033[0m","\n")))

message(cat(paste0("\033[97;", 40, "m","HTML file exported:", "\033[0m",
                   "\033[0;", 32, "m", paste(" make_library_output___", new_lib_folder2, ".html", sep=""),"\033[0m")))

if(is.null(old_lib_csv)==FALSE){
  message(cat(paste0("\033[97;", 40, "m","Combined libary exported:", "\033[0m",
                     "\033[0;", 32, "m",paste(" make_library_output___", new_lib_folder2, ".csv", sep=""),"\033[0m")))
} else {
  message(cat(paste0("\033[97;", 40, "m","New libary exported:", "\033[0m",
                     "\033[0;", 32, "m", paste(" make_library_output___", new_lib_folder2, ".csv", sep=""),"\033[0m")))
}

} #end of function


###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################





