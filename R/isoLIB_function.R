#' @title Generate new strain library or add to existing one.
#' @name isoLIB
#' @rdname isoLIB
#' @description This function creates a strain library by grouping closely related strains of interest based on sequence similarity.
#' For adding new sequences to an already-established strain library, specify the .CSV file path of the older strain library using the 'old_lib_csv" parameter.
#' @export
#' @param input Path of CSV output file from isoTAX step.
#' @param old_lib_csv Optional: Path of CSV output isoLIB file or combined isoLIB file from previous run(s)
#' @param export_html (Default=TRUE) Output the results as an HTML file
#' @param export_csv (Default=TRUE) Output the results as a CSV file.
#' @param include_warnings (Default=FALSE) Whether or not to keep sequences with poor alignment warnings from Step 2 'isoTAX' function. Set TRUE to keep warning sequences, and FALSE to remove warning sequences.
#' @param strain_group_cutoff (Default=0.995) Similarity cutoff (0-1) for delineating between strain groups. 1 = 100% identical/0.995=0.5% difference/0.95=5.0% difference/etc.
#' @param phylum_cutoff Percent cutoff for phylum rank demarcation
#' @param class_cutoff Percent cutoff for class rank demarcation
#' @param order_cutoff Percent cutoff for order rank demarcation
#' @param family_cutoff Percent cutoff for family rank demarcation
#' @param genus_cutoff Percent cutoff for genus rank demarcation
#' @param species_cutoff Percent cutoff for species rank demarcation
#' @seealso \code{\link{isoTAX}}, \code{\link{isoLIB}}
#' @return Returns an isoLIB class object. Default taxonomic cutoffs for phylum (75.0), class (78.5), order (82.0), family (86.5), genus (96.5), and species (98.7) demarcation are based on Yarza et al. 2014, Nature Reviews Microbiology (DOI:10.1038/nrmicro3330)
#' @importFrom utils download.file
#' @importFrom utils untar
#' @importFrom R.utils gunzip
#' @examples
#' #Set path to directory containing example .ab1 files
#' fpath1 <- system.file("extdata/abif_examples/rocket_salad", package = "isolateR")
#'
#' #Step 1: Run isoQC function with default settings
#' isoQC.S4 <- isoQC(input=fpath1)
#'
#' #Step 2: Run isoTAX function with default settings
#' fpath2 <- file.path(fpath1, "isolateR_output/01_isoQC_trimmed_sequences_PASS.csv")
#' isoTAX.S4 <- isoTAX(input=fpath2)
#'
#' #Step 3: Run isoLIB function with default settings
#' fpath3 <- file.path(fpath1, "isolateR_output/02_isoTAX_results.csv")
#' isoLIB.S4 <- isoLIB(input=fpath3)
#'
#' #Show summary statistics
#' isoLIB.S4

isoLIB <- function(input=NULL,
                   old_lib_csv=NULL,
                   export_html=TRUE,
                   export_csv=TRUE,
                   include_warnings=FALSE,
                   strain_group_cutoff = 0.995,
                   phylum_cutoff=75.0,
                   class_cutoff=78.5,
                   order_cutoff=82.0,
                   family_cutoff=86.5,
                   genus_cutoff=96.5,
                   species_cutoff=98.7){


  #Input checks------------------------------------------------------------

  if(!is.null(strain_group_cutoff)){
    if(!is.numeric(strain_group_cutoff)) stop("'strain_group_cutoff' is not numeric. Set between 0-1 (1=no difference i.e., identical sequences)", call.=FALSE)
    if(strain_group_cutoff <0 | strain_group_cutoff >1) stop("Wrong 'strain_group_cutoff' format. Set between 0-1 (1=no difference i.e., identical sequences)", call.=FALSE)
  }

  #------------------------------------------------- functions
  #-------------------------------------------------
  #-------------------------------------------------
  #-------------------------------------------------
  #-------------------------------------------------
  #-------------------------------------------------

  input <- str_replace_all(input, '\\\\', '/')

	new_lib <- str_replace_all(input, '\\\\', '/')
	new_lib_path <- paste(unlist(strsplit(new_lib, '/'))[1:(length(unlist(strsplit(new_lib, '/')))-1)],collapse="/")
	new_lib_file <- read.csv(new_lib) #, row.names = 1) removing row.names setting

  setwd(new_lib_path)

  #-------------------------------------------------
  #:::::::::::::::::::::::::::
  #Download VSEARCH software
  #:::::::::::::::::::::::::::

  path = new_lib_path
  vsearch.path.dl <- file.path(system.file("", package="isolateR"), "vsearch")
  suppressWarnings(dir.create(vsearch.path.dl))
  vsearch_files <- stringr::str_subset(dir(vsearch.path.dl, full.names = FALSE), 'vsearch')


  if(paste(isolateR::get_os())=="windows"){
    if(identical(vsearch_files, character(0))){
      message(cat(paste0("\n", "\033[0;", 32, "m","Operating system is ---> Windows-based <---", "\033[0m")))
      download.file("https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch-2.23.0-win-x86_64.zip", file.path(vsearch.path.dl, 'vsearch-2.23.0-win-x86_64.zip'), mode='wb')
      unzip(file.path(vsearch.path.dl,"vsearch-2.23.0-win-x86_64.zip"),  exdir=file.path(vsearch.path.dl))
      file.copy(file.path(vsearch.path.dl, "vsearch-2.23.0-win-x86_64/bin/vsearch.exe"), file.path(vsearch.path.dl, "vsearch-2.23.0.exe"), overwrite=TRUE)
      unlink(file.path(vsearch.path.dl,"vsearch-2.23.0-win-x86_64"),recursive=TRUE)
      unlink(file.path(vsearch.path.dl,"vsearch-2.23.0-win-x86_64.zip"),recursive=TRUE)
      message(cat(paste0("\n", "\033[0;", 32, "m","Download complete. VSEARCH 2.23.0 has been installed.", "\033[0m", "\n")))
      vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0.exe")
    } else {
      message(cat(paste0("\n", "\033[0;", 32, "m","Operating system is ---> Windows-based <---", "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m","VSEARCH already downloaded", "\033[0m", "\n")))
      vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0.exe")
    }
  }

  if(paste(isolateR::get_os())=="osx-mac"){
    if(identical(vsearch_files, character(0))){
      message(cat(paste0("\033[0;", 32, "m","Operating system is ---> MacOS-based <---", "\033[0m")))
      download.file("https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch-2.23.0-macos-universal.tar.gz", file.path(vsearch.path.dl, 'vsearch-2.23.0-macos-universal.tar.gz'), mode='wb')
      untar(file.path(vsearch.path.dl,"vsearch-2.23.0-macos-universal.tar.gz"),  exdir=file.path(vsearch.path.dl))
      file.copy(file.path(vsearch.path.dl, "vsearch-2.23.0-macos-universal/bin/vsearch"), file.path(vsearch.path.dl, "vsearch-2.23.0_macos"), overwrite=TRUE)
      unlink(file.path(vsearch.path.dl,"vsearch-2.23.0-macos-universal"),recursive=TRUE)
      unlink(file.path(vsearch.path.dl,"vsearch-2.23.0-macos-universal.tar.gz"),recursive=TRUE)
      message(cat(paste0("\033[0;", 32, "m","Download complete. VSEARCH 2.23.0 has been installed.", "\033[0m", "\n")))
      vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0_macos")
    } else {
      message(cat(paste0("\033[0;", 32, "m","Operating system is ---> MacOS-based <---", "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m","VSEARCH already downloaded.", "\033[0m", "\n")))
      vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0_macos")
    }
  }

  if(paste(isolateR::get_os())=="linux"){
    if(identical(vsearch_files, character(0))){
      message(cat(paste0("\033[0;", 32, "m","Operating system is ---> Linux-based <---", "\033[0m")))
      download.file("https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch-2.23.0-linux-x86_64.tar.gz", file.path(vsearch.path.dl, 'vsearch-2.23.0-linux-x86_64.tar.gz'), mode='wb')
      untar(file.path(vsearch.path.dl, "vsearch-2.23.0-linux-x86_64.tar.gz"),  exdir=file.path(vsearch.path.dl))
      file.copy(file.path(vsearch.path.dl, "vsearch-2.23.0-linux-x86_64/bin/vsearch"), file.path(vsearch.path.dl, "vsearch-2.23.0"), overwrite=TRUE)
      unlink(file.path(vsearch.path.dl,"vsearch-2.23.0-linux-x86_64"),recursive=TRUE)
      unlink(file.path(vsearch.path.dl,"vsearch-2.23.0-linux-x86_64.tar.gz"),recursive=TRUE)
      message(cat(paste0("\033[0;", 32, "m","Download complete. VSEARCH 2.23.0 has been installed.", "\033[0m", "\n")))
      vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0")
    } else {
      message(cat(paste0("\033[0;", 32, "m","Operating system is ---> Linux-based <---", "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m","VSEARCH already downloaded.", "\033[0m", "\n")))
      vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0")
    }

  }

  #::::::::::::
  #Stage paths
  #::::::::::::
  message(cat(paste0("\n", "\033[0;", 32, "m","Staging files for de-replication", "\033[0m", "\n")))

  temp_csv <- "03_isoLIB_temp.csv"
  temp_fasta <- "03_isoLIB_temp.fasta"
  temp_drep <- paste0("03_isoLIB_temp.drep", sep="")


  if(include_warnings == TRUE){
    write.csv(new_lib_file, file=paste(temp_csv, sep=""))
    new_lib_file <- read.csv(paste(temp_csv, sep=""), row.names = 1)
    new_lib_file_path <- paste(temp_csv, sep="")
  } else {
    invisible(capture.output(capture.output(new_lib_file <- tryCatch(
      {
        new_lib_file %>% filter(!grepl("Poor", warning)) %>% filter(!grepl("Fail", decision))
      },
      error = function(e) {
        new_lib_file %>% filter(!grepl("Poor", warning))
      }), type="message"), type="output"))
    if(nrow(new_lib_file)==0)stop("There are no sequences after removing warning sequences. Set 'include_warnings=TRUE' to proceed", call.=FALSE)
    write.csv(new_lib_file, file=paste(temp_csv, sep=""))
    new_lib_file <- read.csv(paste(temp_csv, sep=""), row.names = 1)
    new_lib_file_path <- paste(temp_csv, sep="")
  }

  make_fasta((file.path(getwd(), new_lib_file_path)), col_names="filename", col_seqs="seqs_trim", output=paste(temp_fasta, sep="")) # exports as "output.fasta"

  #------------------------------------------------- Trying with global searching instead...
  system2(vsearch.path, paste(" --usearch_global ", temp_fasta, " --db ", temp_fasta, " --uc_allhits --uc ", temp_drep, " --id ", strain_group_cutoff, " --query_cov 0.95 --maxaccepts 0 --maxrejects 0 --strand plus ", sep=""), stdout="", stderr="")
  #readline
  #-------------------------------------------------

  #:::::::::::::::::
  #Organize results
  #:::::::::::::::::
  drep.results <- read.csv(temp_drep, sep="\t", header = FALSE) %>%
    arrange(V2) %>%
    group_by(V9) %>%
    mutate(grouping = dplyr::first(V2)) %>%
    select(V9, grouping) %>%
    ungroup() %>%
    dplyr::rename("filename" = "V9") %>%
    distinct(filename, .keep_all=TRUE)

  unlink(file.path(temp_fasta)) #Remove temp FASTA
  unlink(file.path(temp_drep)) #Remove temp DREP
  unlink(file.path(temp_csv)) #Remove temp CSV

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
    arrange(desc(length_trim)) %>% #filter(grepl("1015|HB534", filename))
    arrange(desc(override_priority)) %>%
    group_by(grouping) %>%
    mutate(strain_group = dplyr::first(filename)) %>%
    ungroup() %>%
    mutate(ref_strain = ifelse(strain_group == filename, "yes", "no"))

  #Subset only columns of interest------------------------------------------------------------------

  merged.drep1.sub <- merged.drep1 %>% select(strain_group, date, filename, seqs_trim, phred_trim, Ns_trim, length_trim,
                                              closest_match, NCBI_acc, ID,
                                              rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, ref_strain)


  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  #Combining OLD database if provided---------------------------------------------------------
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  if(is.null(old_lib_csv)==FALSE){
    old_lib_file <- read.csv(old_lib_csv)#, row.names=1)
    message(cat(paste0("\n", "\033[0;", 32, "m","Checking old/new library CSV files have same dimensions", "\033[0m")))
    message(cat(paste0("\033[0;", 32, "m",paste("No. of columns in old library: ", ncol(old_lib_file), sep=""), "\033[0m")))
    message(cat(paste0("\033[0;", 32, "m",paste("No. of columns in old library: ", ncol(merged.drep1.sub), sep=""), "\033[0m")))
    if(!ncol(old_lib_file)==ncol(merged.drep1.sub)) stop('Library file dimensions are not the same and require manual inspection.', call.=FALSE)
    message(cat(paste0("\033[0;", 32, "m","Same dimensions detected, proceeding...", "\033[0m", "\n")))
    #

    #Check for duplicate names
    #dup.names <- merged.drep1.sub$filename[merged.drep1.sub$filename %in% old_lib_file$filename]
    #dup.names.index <- setNames(paste(str_split_fixed(dup.names, "[.]", 2)[,1], "_rep.", str_split_fixed(dup.names, "[.]", 2)[,2], sep=""), dup.names)
    #merged.drep1.sub <- merged.drep1.sub %>% dplyr::recode(filename, !!!dup.names.index) # doesn't work - can't handle periods (.)

    #combine old and new isoLIB files-----------------------
    combined_lib_file <- rbind(old_lib_file, merged.drep1.sub) #old_lib_file must be at top so older seqs take priority

    #make unique strain names incase any duplicated-----------------------
    combined_lib_file$unique_filename <- paste(str_pad(1:nrow(combined_lib_file), 4, pad="0"), combined_lib_file$filename, sep="_")

    #-------------------------------------------------------
    temp_combined_fasta <- "03_isoLIB_combined_temp.fasta"
    temp_combined_drep <- paste0("03_isoLIB_combined_temp.drep", sep="")
    temp_combined_csv <- "03_isoLIB_combined_temp.csv"
    #-------------------------------------------------------
    write.csv(combined_lib_file, file= temp_combined_csv)

    #De-replicating-------------------------------------------------------
    make_fasta((file.path(getwd(), temp_combined_csv)), col_names="unique_filename", col_seqs="seqs_trim", output=paste(temp_combined_fasta, sep="")) # exports as "output.fasta"

    #-------------------------------------------------


    #system2(vsearch.path, paste(" --derep_prefix ", "output.fasta", " --output ", out.fasta, " --uc ", temp_drep, " --sizein --sizeout", sep=""), stdout="", stderr="")
    system2(vsearch.path, paste(" --usearch_global ", temp_combined_fasta, " --db ", temp_combined_fasta, " --uc ", temp_combined_drep, " --id ", strain_group_cutoff, " --query_cov 0.95 --maxaccepts 0 --maxrejects 0 --top_hits_only --strand plus ", sep=""), stdout="", stderr="")

    #-------------------------------------------------
    #:::::::::::::::::
    #Organize results
    #:::::::::::::::::
    drep.results <- read.csv(temp_combined_drep, sep="\t", header = FALSE) %>%
      arrange(V2) %>%
      group_by(V9) %>%
      mutate(grouping = dplyr::first(V2)) %>%
      select(V9, grouping) %>%
      ungroup() %>%
      dplyr::rename("unique_filename" = "V9") %>%
      distinct(unique_filename, .keep_all=TRUE)

    unlink(file.path(temp_combined_fasta)) #Remove temp FASTA
    unlink(file.path(temp_combined_drep)) #Remove temp DREP
    unlink(file.path(temp_combined_csv)) #Remove temp CSV

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
      ungroup() %>%
      mutate(ref_strain = ifelse(strain_group == unique_filename, "yes", "no")) %>%
      #Remove unique identifier characters
      mutate(strain_group = sapply(1:length(.$strain_group),function(x) str_sub(.$strain_group[x], 6, nchar(.$strain_group[x]))))


    merged.drep1.sub <- merged.drep2 %>% select(strain_group, date, filename, seqs_trim, phred_trim, Ns_trim, length_trim,
                                                closest_match, NCBI_acc, ID,
                                                rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, ref_strain)
  } #-----------end of old library add-in


  #merged.drep1 %>% group_by(grouping) %>% mutate(strain_group = last(filename))
  #merged.drep1.sub <- merged.drep1 %>% select(grouping,strain_group) %>% mutate(same = case_when(grouping == strain_group ~ "same", TRUE ~ "other"))
  #merged.drep1 <- merged.drep


  ###########################################################################################################
  ###########################################################################################################
  ###########################################################################################################
  ###########################################################################################################
  ###########################################################################################################

  #::::::::::::::::::::::::::::::::::::::::::::::
  # Set cutoffs for taxonomic rank demarcation
  #::::::::::::::::::::::::::::::::::::::::::::::

  merged.drep1.sub <- merged.drep1.sub %>%
    mutate(phylum_cutoff= phylum_cutoff) %>%
    mutate(class_cutoff= class_cutoff) %>%
    mutate(order_cutoff= order_cutoff) %>%
    mutate(family_cutoff= family_cutoff) %>%
    mutate(genus_cutoff= genus_cutoff) %>%
    mutate(species_cutoff= species_cutoff)


  ###############
  # Output
  ###############

  #Set S4 class for outputs
  isolib.S4 <- df_to_isoLIB(merged.drep1.sub)

  #--------------------------------------------------

  #Export messages
  message(cat(paste0("\033[97;", 40, "m","'isoLIB' steps completed. Exporting files...", "\033[0m")))

  #Export HTML file----------------------------------------------------------------------
  if(export_html==TRUE){
    output_name <- "03_isoLIB_results.html"
    export_html(isolib.S4)
    #merged.df3_react <- export_html(out_df, path, output=output_name)
    message(cat(paste0("\033[97;", 40, "m","HTML results exported: ", "\033[0m",
                       "\033[0;", 32, "m", " ", file.path(path, output_name),"\033[0m", "\n")))
    } else if(export_html==FALSE){}

  #export CSV files----------------------------------------------------------------------
  if(export_csv==TRUE){
    output_name <- "03_isoLIB_results.csv"
    csv_output <- S4_to_dataframe(isolib.S4) %>% select(strain_group, date, filename,  seqs_trim, phred_trim, Ns_trim, length_trim,
                                                        closest_match, NCBI_acc, ID,
                                                        rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, ref_strain)
    write.csv(csv_output, file=paste(output_name, sep=""), row.names = FALSE)
    message(cat(paste0("\033[97;", 40, "m","CSV file exported:", "\033[0m",
                       "\033[0;", 32, "m", " ", file.path(path, output_name),"\033[0m",
                       "\033[0;", 31, "m", "  <--- Final strain library output","\033[0m")))
  } else if(export_csv==FALSE){}


  return(isolib.S4)

} #end of function

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################





