#abif_fasta2 function
#------------------------------------------------------------------------------------
abif_fasta2 <- function(folder=NULL,
                        export_html=TRUE,
                        export_csv=TRUE,
                        export_fasta=TRUE,
                        export_fasta_revcomp=FALSE,
                        verbose=TRUE,
                        exclude=NULL,
                        min_phred_score=1,
                        quality_cutoff = 25,
                        sliding_window_size = 15,
                        date=NULL,
                        files_manual=NULL){
  
  # function requirements------------------------------------------------------------
  #checking for required packages; installing those not yet installed
  suppressPackageStartupMessages({suppressWarnings({
  if(require(dplyr)==FALSE) install.packages('dplyr')
  if(require(svMisc)==FALSE) install.packages('svMisc')
  if(require(Biostrings)==FALSE){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("Biostrings", update = FALSE)
  }
  if(require(seqinr)==FALSE) install.packages('seqinr')
  if(require(stringr)==FALSE) install.packages('stringr')
  if(require(sangerseqR)==FALSE){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("sangerseqR", update = FALSE)
  }
  if(require(sangeranalyseR)==FALSE){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("sangeranalyseR", update = FALSE)
  }
  if(require(reactable)==FALSE) install.packages('reactable')
  if(require(reactablefmtr)==FALSE) install.packages('reactablefmtr')
  if(require(pander)==FALSE) install.packages('pander')
  if(require(remotes)==FALSE) install.packages('remotes')
  if(require(dataui)==FALSE){ remotes::install_github("timelyportfolio/dataui", upgrade="never")}

    
  #  if (!require("BiocManager", quietly = TRUE))
  #    install.packages("remotes", update = FALSE)
  #  remotes::install_github("timelyportfolio/dataui", upgrade="never")
  #}
  
  #loading required packages
  library(dplyr)
  library(svMisc)
  library(remotes)
  library(Biostrings)
  library(seqinr)
  library(stringr)
  library(sangerseqR)
  library(sangeranalyseR)
  library(reactable)
  library(reactablefmtr)
  library(pander)
  library(dataui)
  }) #end of suppressWarnings
    }) #end suppressPackageStartupMessages
  
  # Reading files--------------------------------------------------------------------
  # Check folder path
  if(is.null(folder)) stop('Name of folder not supplied. Please check the execution script preamble to make sure it is entered correctly', call.=FALSE)
  if(is.null(sliding_window_size)) {
    sliding_window_size = 15
    message(cat(paste0("\n", "\033[97;", 40, "m","Setting defaults for null inputs:", "\033[0m", "\n")))
    message(cat(paste0("\n", "\033[97;", 40, "m","sliding_window_size = 15", "\033[0m", "\n")))
  }
  
  #Setting folder paths for abif files----------------------------------------------------
  path <- str_replace_all(folder, '\\\\', '/')
  fname <- list.files(path)
  pattern <- 'ab1'
  abif_files <- str_subset(fname, pattern)
  if(is.null(files_manual)==FALSE){abif_files <- files_manual}

  #excluding specified files
  if(!is.null(exclude)) {
    
    # first check if files specified in exclude are in abif_files
    if(any(!exclude %in% abif_files)) {
      absent <- exclude[!exclude %in% abif_files]
      
      msg <- sprintf('%s could not be excluded because they were not found in your input folder. Only files in the input folder can be excluded.',
                     str_c(absent, collapse=', '))
      stop(msg, call.=FALSE)
    }
    else{
      abif_files <- abif_files[!abif_files %in% exclude]
      msg <- sprintf('%s have been excluded from your fasta file.',
                     str_c(exclude, collapse=', '))
      message(msg)
    }
  }
  
  #Checking overwrite inputs-------------------------------------------------
  if(!is.null(date)){
    date.formatted <- paste(str_pad(str_split_fixed(date, "_", 3)[,1], 2, pad = "0"),
                            str_pad(str_split_fixed(date, "_", 3)[,2], 2, pad = "0"),
                            str_pad(str_split_fixed(date, "_", 3)[,3], 2, pad = "0"), sep="_")
    message(cat(paste0("\n", "\033[97;", 40, "m","Setting date (YY/MM/DD) as: ",date.formatted, "\033[0m", "\n")))
  }
  
  #reading abif files----------------------------------------------------------------
  #initializing objects
  dseq <- list()
  trim.start <- c()
  trim.end <- c()
  rawsparklist <- list()
  trimsparklist <- list()
  rawlist <- list()
  trimlist <- list()
  checkseq <- c()
  seq.warnings <- c()
  date.list <- c()
  length.list <- c()
  
  
  message(cat(paste0("\n", "\033[97;", 95, "m","Trimming ", paste(length(abif_files)) ," input files", "\033[0m", "\n")))
  
  prog.bar.x <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(abif_files), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       #width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=") 
  
  #Grap max length 
  for(i in 1:length(abif_files)){
    fpath <- file.path(path, abif_files[i])
    length.list <- c(length.list, as.numeric(nchar(unlist((read.abif(fpath))@data['PBAS.2'])))+5)
  }
  length.list.max <- max(length.list)
  
  
  for(i in 1:length(abif_files)){
    
    fpath <- file.path(path, abif_files[i])
    
    # Basecalling
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    #if(basecall==FALSE ){
    #  abif <- read.abif(fpath)
    #  rawseq <- as.character(abif$Data['PBAS.1'])
    #}
    #if (basecall==TRUE){
    #abif <- makeBaseCalls(readsangerseq(fpath), ratio=0.33)
    #rawseq <- primarySeq(abif, string=TRUE)
    
    abif.pre <- read.abif(fpath)
    quality_cutoff_max <- max(unlist(abif.pre@data['PCON.1']))
    
    if(!is.null(date)){
        date.list <- c(date.list, date.formatted)
        } else {
        date.list <- c(date.list, paste(abif.pre@data['RUND.1'][[1]][[1]],
                                str_pad(abif.pre@data['RUND.1'][[1]][[2]], 2, pad = "0"),
                                str_pad(abif.pre@data['RUND.1'][[1]][[3]], 2, pad = "0"), sep="_"))
        }
    
    if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Trimming file: ", fpath, "\033[0m")))}
               
    if(verbose==FALSE){      
    invisible(capture.output(capture.output(abif1 <- tryCatch( 
      {
             SangerRead(readFeature           = "Forward Read",
                        readFileName          = fpath,
                        geneticCode           = GENETIC_CODE,
                        TrimmingMethod        = "M2",
                        M1TrimmingCutoff      = NULL,
                        M2CutoffQualityScore  = round(quality_cutoff_max*0.66),
                        M2SlidingWindowSize   = sliding_window_size,
                        baseNumPerRow         = 100,
                        heightPerRow          = 200,
                        signalRatioCutoff     = 0.33)
      },
    error = function(e) {
      message(cat(paste0("\033[97;", 41, "m","First trimming attempt failed, retrying with 'sliding window = 1'", "\033[0m")))
                  SangerRead(readFeature           = "Forward Read",
                             readFileName          = fpath,
                             geneticCode           = GENETIC_CODE,
                             TrimmingMethod        = "M2",
                             M1TrimmingCutoff      = NULL,
                             M2CutoffQualityScore  = round(quality_cutoff_max*0.66),
                             M2SlidingWindowSize   = 1,
                             baseNumPerRow         = 100,
                             heightPerRow          = 200,
                             signalRatioCutoff     = 0.33)
    }
    ), type="message"), type="output"))
    } else {abif1 <- tryCatch( 
      {
        SangerRead(readFeature           = "Forward Read",
                   readFileName          = fpath,
                   geneticCode           = GENETIC_CODE,
                   TrimmingMethod        = "M2",
                   M1TrimmingCutoff      = NULL,
                   M2CutoffQualityScore  = round(quality_cutoff_max*0.66),
                   M2SlidingWindowSize   = sliding_window_size,
                   baseNumPerRow         = 100,
                   heightPerRow          = 200,
                   signalRatioCutoff     = 0.33)
      },
      error = function(e) {
        message(cat(paste0("\033[97;", 41, "m","First trimming attempt failed, retrying with 'sliding window = 1'", "\033[0m")))
        SangerRead(readFeature           = "Forward Read",
                   readFileName          = fpath,
                   geneticCode           = GENETIC_CODE,
                   TrimmingMethod        = "M2",
                   M1TrimmingCutoff      = NULL,
                   M2CutoffQualityScore  = round(quality_cutoff_max*0.66),
                   M2SlidingWindowSize   = 1,
                   baseNumPerRow         = 100,
                   heightPerRow          = 200,
                   signalRatioCutoff     = 0.33)
      }
    )}
    
    #start=abif1@QualityReport@trimmedStartPos
    abif1.start = NULL
    #abif1.start <- abif1@QualityReport@trimmedStartPos
    abif1.start <- abif1@QualityReport@trimmedStartPos + 5 
    if(abif1.start >= abif1@QualityReport@rawSeqLength){abif1.start <- abif1@QualityReport@trimmedStartPos}
    if(abif1.start >= abif1@QualityReport@trimmedFinishPos){abif1.start <- abif1@QualityReport@trimmedStartPos}
    
    abif1.end <- nchar(paste(abif1@primarySeq))
    #if(is.null(abif1.start) == TRUE){abif1.start <- 50}
    #if(abif1.start >= abif3.end){abif1.start <- 1}
    
    abif2 <- DNAStringSet(abif1@primarySeq, start=abif1.start, end=abif1.end)
    #abif2 <- Biostrings::reverseComplement(DNAStringSet(paste(abif1@primarySeq),
    #                                                    start=abif1.start,
    #                                                    end=abif1.end))
    
    if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Searching for reverse primer endpoint: ", fpath, "\033[0m")))}
    
    if(verbose==FALSE){      
      invisible(capture.output(capture.output(abif3 <- DECIPHER::TrimDNA(DNAStringSet(paste(abif2)),
                                                       leftPatterns="",
                                                       rightPatterns="CTGCTGCCTYCCGTA",
                                                       minWidth=1,
                                                       maxDistance = 0.2,
                                                       minOverlap = 15,
                                                       allowInternal=TRUE,
                                                       threshold = 1,
                                                       maxAverageError = 1,
                                                       maxAmbiguities = 0.5,
                                                       type="both"), type="message"), type="output"))
    } else {
      abif3 <- DECIPHER::TrimDNA(DNAStringSet(paste(abif2)),
                                 leftPatterns="",
                                 rightPatterns="CTGCTGCCTYCCGTA",
                                 minWidth=1,
                                 maxDistance = 0.2,
                                 minOverlap = 15,
                                 allowInternal=TRUE,
                                 threshold = 1,
                                 maxAverageError = 1,
                                 maxAmbiguities = 0.5,
                                 type="both")
    }
    
    #This code chunk checks if trimmed sequence region looks okay before proceeding to next step
    #----------------------------------------------------------------------------------------------
    if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Removing anything beyond reverse primer: ", fpath, "\033[0m")))}
    abif3.end <- end(abif3[[1]]) + nchar("CTGCTGCCTYCCGTA")
    if(abif3.end > width(abif3[[1]])){abif3.end <- width(abif3[[1]])}
    if(width(abif3[[1]]) <= 300){abif3.end <- nchar(paste(abif2)) }
    if(abif3.end > (abif1@QualityReport@trimmedFinishPos - abif1.start)){abif3.end <- (abif1@QualityReport@trimmedFinishPos - abif1.start)}
    if(abif3.end > length.list.max){abif3.end <- length.list.max}
    if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Trimming error check: 1 (Passed)", "\033[0m")))}
    abif1.start.new <- abif1.start
    if(abif1.start <= 50 & abif3.end >= 51){abif1.start.new <- 50}
    #abif <- DNAStringSet(abif2, start=1, end = abif3.end) # equivalent to line below
    abif <- DNAStringSet(abif1@primarySeq, start=abif1.start.new, end = (abif1.start + abif3.end -1))
    if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Trimming error check: 2 (Passed)", "\033[0m")))}
    abif.scores.xx <- abif1@QualityReport@qualityPhredScores[abif1.start.new:(abif1.start + abif3.end -1)]
    if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Trimming error check: 3 (Passed)", "\033[0m")))}
    abif.scores <- c(rep(0.1, abif1.start.new), abif.scores.xx, rep(0.1, length.list.max-(((abif1.start + abif3.end) -1))))
    trim.start <- c(trim.start,abif1.start.new)
    trim.end <- c(trim.end,(abif1.start + abif3.end -1))
    rawseq <- as.character(paste(abif))
    
    if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Saving seq info: ", fpath, "\033[0m")))}
    
    if(verbose==TRUE){
      #print a message to mark progress
      msg <- sprintf('Processing "%s"', abif_files[i])
      message(msg)
    }
  

    #:::::::::::::::::::::::::::::::
    # Save seq + quality scores
    #:::::::::::::::::::::::::::::::
    #rawsparklist[[i]] <- c(abif1@QualityReport@qualityPhredScores, rep(0.1, length.list.max-length(abif1@QualityReport@qualityPhredScores)))
    rawspark.adj <- abif1@QualityReport@qualityPhredScores
    if(length(abif1@QualityReport@qualityPhredScores) >= length.list.max){rawspark.adj <- rawspark.adj[1:length.list.max]}
    rawsparklist[[i]] <- c(rawspark.adj, rep(0.1, length.list.max-length(rawspark.adj)))
    trimsparklist[[i]] <- abif.scores
    rawlist[[i]] <- paste(abif1@primarySeqRaw)
    trimlist[[i]] <- paste(abif)
    #rawlist[[i]] <- paste(Biostrings::reverseComplement(DNAStringSet(paste(abif1@primarySeq))))
    #trimlist[[i]] <- paste(Biostrings::reverseComplement(DNAStringSet(abif)))
    
    # building trim check for export-------------------------------------------------
    entry <- data.frame(date=date.list[i],
                        filename=abif_files[i],
                        #spark_data = list(sparklist[[i]]),
                        phred_raw = mean(rawsparklist[[i]]),
                        phred_trim = mean(abif.scores.xx),
                        length_raw = nchar(rawlist[[i]]),
                        length_trim = nchar(trimlist[[i]]),
                        Ns_raw = str_count(rawlist[[i]], "[Nn]"),
                        Ns_trim = str_count(trimlist[[i]], "[Nn]"),
                        raw_seq = rawlist[[i]],
                        trim_seq = trimlist[[i]])

    if(entry$length_trim < 200) {
      seq.warnings <- c(seq.warnings, paste(abif_files[i]))
      #msg <- sprintf('The file "%s" has a length of <200 bp after trimming. Manual check of sequence is recommended.', abif_files[i])
      #message(msg)
      #warning(msg, call.=FALSE)
      #
      #  seq_pass <- FALSE
      #  failreason <- paste(failreason, msg)
    }
    # building trim check for export-------------------------------------------------
    checkseq <- rbind(checkseq, entry)
    #ntseq <- trimlist[[i]]
    #svMisc::progress(((i/length(abif_files))*100),  max.value = 100, progress.bar = TRUE, init=10)
    if(verbose==FALSE){setTxtProgressBar(prog.bar.x, i)}
    Sys.sleep(0.01)
  }
  # end of file loop
  Sys.sleep(0.02)
  
  #:::::::::::::::::::::::::::::::
  #List function for table below
  #:::::::::::::::::::::::::::::::
  list_with_names<-function(...){dplyr::lst(...)}
  uni.font <- 11
  
  # Make phred scores into lists for sparkline formatting (V1=raw, V2=trim)
  sparkdf <- c()
  for(i in 1:length(rawsparklist)){
    sparkdf.x.raw <- as.data.frame(rbind(list(rawsparklist[[i]]))) %>% `colnames<-`("raw")
    sparkdf.x.trim <- as.data.frame(rbind(list(trimsparklist[[i]]))) %>% `colnames<-`("trim")
    sparkdf <- rbind(sparkdf, cbind(sparkdf.x.raw,sparkdf.x.trim))
  }
  # Add spark data to checkseq table and round long numbers to 2 decimal points
  checkseq.sub <- checkseq %>% mutate(spark_raw = sparkdf$raw) %>%
    mutate(spark_trim = sparkdf$trim) %>%
    mutate(phred_raw = round(checkseq$phred_raw, 2)) %>%
    mutate(phred_trim = round(checkseq$phred_trim, 2)) %>%
    mutate(arrcol = rep("->")) %>%
    mutate(decision = ifelse(filename %in% seq.warnings, "Fail", "Pass")) %>%
    mutate(decision = ifelse(phred_trim < min_phred_score, "Fail", decision)) %>%
    mutate(decision_col = ifelse(decision == "Fail", grDevices::adjustcolor("#B20A2C", alpha.f=0.7), grDevices::adjustcolor("#CCE6CC", alpha.f=1))) %>%
    select(date, filename, raw_seq, phred_raw, Ns_raw, length_raw, spark_raw, arrcol,
           trim_seq, phred_trim, Ns_trim, length_trim, decision, decision_col) #%>%
  #`colnames<-`(c("filename", "rawSeq", "rawQ", "rawN", "rawLength", "rawSpark", "->",
  #               "trmSeq", "trmQ", "trmN", "trmLength", "trmSpark"))
  
  checkseq.sub$spark_raw <- lapply(checkseq.sub$spark_raw, function(x) list(x))
  
  #trim.start <- c(rep(1:40))
  #trim.end <- c(rep(500:550))
  #abif1.df.spark <- sparkline(abif1@QualityReport@qualityPhredScores, width=200, fillColor=FALSE)
  

  checkseq_react <- reactable(checkseq.sub,
                              fullWidth=FALSE, searchable = TRUE, bordered = TRUE, resizable =TRUE, #width = 1200,
                              defaultPageSize=200, highlight = TRUE, showSortable = TRUE, compact=TRUE, wrap = TRUE, #highlight=TRUE,
                              #theme = reactableTheme(headerStyle = list(wrap=TRUE)),
                              theme = reactableTheme(                  #borderColor = "#555",
                                headerStyle = list(
                                  "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                                  "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)"))),
                              #borderColor = "#555")),
                              #defaultColDef = colDef(minWidth = 100), #
                              #defaultColDef = colDef(style = reactablefmtr::cell_style(checkseq.sub, font_size=12)),
                              style = list(minWidth = 1400),
                              
                              defaultColDef = colDef(style = list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif"), minWidth= 75, align="center"),
                              #defaultColDef = colDef(style = "white-space: nowrap; font-size: 16px; font-family: Courier New", minWidth= 75, align="center"),
                              #  headerStyle = list(style = "white-space: wrap;", wrap=FALSE)),
                              #--------------------
                              columns = list(
                                #--------------------
                                #arrcol = colDef(style = "white-space: nowrap; font-size: 20px; fontWeight: bold;", minWidth= 150, align="center"), #fontSize = '14px',"white-space: nowrap;"
                                # Spark lines
                                #-------------------
                                spark_raw = colDef(name = "Quality sparkline", width = 200, cell = function(value, index) {
                                  dui_sparkline(data = value[[1]], height = 25, min = min(unlist(checkseq.sub$spark_raw)), max = max(unlist(checkseq.sub$spark_raw)),
                                                margin = list( top= 2, right= 5, bottom= 2, left= 5 ),
                                                components = list(                                                    #dui_sparkpatternlines(id = "band_pattern_misc", height = 3, width = 3, stroke = "green", fillOpacity=0.9, strokeWidth = 0.7, orientation = list('diagonal')  ),
                                                  dui_sparkbandline(band = list( from = list( x = trim.start[[index]] ), to = list( x = trim.end[[index]] ) ), fill = "green", fillOpacity=0.2 ), #, fillOpacity=0.3 #use if fill not banding "url(#band_pattern_misc)"
                                                  dui_sparklineseries(strokeWidth = 0.5,stroke = "black")
                                                )
                                  )
                                }
                                ),
                                #spark_raw = colDef(width= 150, cell = function(value, index) {
                                #  sparkline(checkseq.sub$spark_raw[[index]], width=150, chartRangeMin = min(unlist(checkseq.sub$spark_raw)), chartRangeMax = max(unlist(checkseq.sub$spark_raw)), fillColor=FALSE)   #sparkline(checkseq.sub$spark_raw[[index]], width=checkseq$length_raw[[index]]/max(checkseq$length_raw,fillColor=FALSE)
                                #        }),
                                #spark_trim = colDef(width= 150, cell = function(value, index) {
                                #  sparkline(checkseq.sub$spark_trim[[index]], width=150, chartRangeMin = min(unlist(checkseq.sub$spark_trim)), chartRangeMax = max(unlist(checkseq.sub$spark_trim)), fillColor=FALSE)
                                #        }),
                                date = colDef(name = "Date", width=75, style=list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                                raw_seq = colDef(name = "rawSeq", style=list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                                trim_seq = colDef(name = "trimSeq", style=list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                                arrcol = colDef(name = ""),
                                filename = colDef(minWidth= 125, style = list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                                decision = colDef(name="Decision", vAlign = "center",  minWidth=170, # headerVAlign = "bottom",      #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                style = function(value, index, decision){
                                                  if(checkseq.sub$decision_col[[index]] == grDevices::adjustcolor("#CCE6CC", alpha.f=1)) { 
                                                    color.x <- "black" 
                                                  } else if(!checkseq.sub$decision_col[[index]] == grDevices::adjustcolor("#CCE6CC", alpha.f=1)) {
                                                    color.x <- "white" }
                                                  list(background = checkseq.sub$decision_col[[index]], color = color.x, 
                                                       whiteSpace = "nowrap",
                                                       fontSize=uni.font+6, 
                                                       #fontWeight="bold",
                                                       align="center",
                                                       fontFamily = "sans-serif")}),
                                decision_col = colDef(show=FALSE),
                                # Heatmap widgets
                                #-------------------
                                #arrcol = colDef(width=350, cell = color_scales(checkseq.sub, text_size=20)),
                                #Ns_raw = colDef(width=300, style = color_scales(checkseq.sub, colors = c("white", "#B20A2C"))),
                                #Ns_raw = colDef(width=300, style = color_scales(checkseq.sub, colors = c("white", "#B20A2C"))),
                                # Bubble widgets
                                #-------------------
                                #Ns_trim = colDef(cell = bubble_grid(checkseq.sub, colors = c("white", "#B20A2C"),
                                #                                    number_fmt = scales::number,
                                #                                    min_value=0,
                                #                                    max_value=20),
                                #                 align="center"),
                                # Gauge widgets
                                #-------------------
                                phred_raw = colDef(name = "rawQ", width= 60, cell = gauge_chart(checkseq.sub,
                                                                                                fill_color = c('#e5f5e0', '#a1d99b', '#31a354'), opacity = 0.5,
                                                                                                #viridis::viridis(5)[1:4]),
                                                                                                #'#D7191C','#FDAE61','#FFFFBF','#A6D96A','#1A9641'),
                                                                                                number_fmt = scales::comma,
                                                                                                #bold_text = TRUE,
                                                                                                text_size = 12,
                                                                                                max_value = max(c(checkseq.sub$phred_raw,checkseq.sub$phred_trim)),
                                                                                                show_min_max = FALSE)),
                                phred_trim = colDef(name = "trimQ", width= 60, cell = gauge_chart(checkseq.sub,
                                                                                                  fill_color = c('#e5f5e0', '#a1d99b', '#31a354'), opacity = 0.5,
                                                                                                  #viridis::viridis(5)[1:4]),
                                                                                                  #'#D7191C','#FDAE61','#FFFFBF','#A6D96A','#1A9641'),
                                                                                                  number_fmt = scales::comma,
                                                                                                  #bold_text = TRUE,
                                                                                                  text_size = 12,
                                                                                                  max_value = max(c(checkseq.sub$phred_raw,checkseq.sub$phred_trim)),
                                                                                                  show_min_max = FALSE)),
                                Ns_raw = colDef(name = "rawNs", width= 60, cell = gauge_chart(checkseq.sub,
                                                                                              fill_color = c("white", "#B20A2C"), opacity = 0.5,
                                                                                              #viridis::viridis(5)[1:4]),
                                                                                              #'#D7191C','#FDAE61','#FFFFBF','#A6D96A','#1A9641'),
                                                                                              number_fmt = scales::comma,
                                                                                              #bold_text = TRUE,
                                                                                              text_size = 12,
                                                                                              max_value = max(c(checkseq.sub$Ns_raw,checkseq.sub$Ns_trim)),
                                                                                              show_min_max = FALSE)),
                                Ns_trim = colDef(name = "trimNs", width= 60, cell = gauge_chart(checkseq.sub,
                                                                                                fill_color = c("white", "#B20A2C"), opacity = 0.5,
                                                                                                #viridis::viridis(5)[1:4]),
                                                                                                #'#D7191C','#FDAE61','#FFFFBF','#A6D96A','#1A9641'),
                                                                                                number_fmt = scales::comma,
                                                                                                #bold_text = TRUE,
                                                                                                text_size = 12,
                                                                                                max_value = max(c(checkseq.sub$Ns_raw,checkseq.sub$Ns_trim)),
                                                                                                show_min_max = FALSE)),
                                # Barplot widgets
                                #-------------------
                                length_raw = colDef(name = "rawLength", width=150, cell = data_bars(checkseq.sub,
                                                                                                    text_position = "inside-base",
                                                                                                    text_size = 12,
                                                                                                    #fill_color="#6f86ab",
                                                                                                    fill_color = c('#FFF2D9','#FFE1A6','#FFCB66','#FFB627'),
                                                                                                    fill_gradient = TRUE,
                                                                                                    background = 'transparent',
                                                                                                    number_fmt = scales::comma_format(accuracy = 0.1),
                                                                                                    round_edges = FALSE, align_bars="left")),
                                #checkseq.sub, text_position = "none", fill_color="red", round_edges = FALSE, align_bars="left")),
                                length_trim = colDef(name = "trimLength", width=150, cell = data_bars(checkseq.sub,
                                                                                                      text_position = "inside-base", #inside-end / outside-base / above
                                                                                                      text_size = 12,
                                                                                                      #fill_color="#6f86ab",
                                                                                                      fill_color = c('#FFF2D9','#FFE1A6','#FFCB66','#FFB627'),
                                                                                                      fill_gradient = TRUE,
                                                                                                      background = 'transparent',
                                                                                                      number_fmt = scales::comma_format(accuracy = 0.1),
                                                                                                      round_edges = FALSE, align_bars="left")) )) %>% add_title(paste("Output for project: ", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))], sep=""))  #unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))]

 
  
  #building output file--------------------------------------------------------------
  if(is.null(seq.warnings)==FALSE){
  seq.warnings.txt <- paste(seq.warnings, collapse="\r\n     ")
  message(cat(paste0("\033[97;", 40, "m","\r\nThe following sequences failed with a length of <200 bp after trimming:","\033[0m","\n     ",
                     "\033[0;", 95, "m", seq.warnings.txt,"\033[0m","\n")))
  }
  
  message(cat(paste0("\n", "\033[97;", 40, "m","Export directory:", "\033[0m",
                     "\033[0;", 32, "m", " ", file.path(path, "output"), "\033[0m","\n")))
  
  suppressWarnings(dir.create(file.path(path, "output")))
  fname <- file.path(path,"output", paste0("abif_fasta2_output_PASS+FAIL___", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))], sep=""))
  fname_pass <- file.path(path,"output", paste0("abif_fasta2_output_PASS___", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))], sep=""))
  fname_fail <- file.path(path,"output", paste0("abif_fasta2_output_FAIL___", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))], sep=""))
  checkseq.sub.pass <- checkseq.sub %>% select(-spark_raw) %>% filter(!filename %in% seq.warnings)
  checkseq.sub.fail <- checkseq.sub %>% select(-spark_raw) %>% filter(filename %in% seq.warnings)
  
  #export HTML file----------------------------------------------------------------------
  if(export_html == TRUE) {
    fname_html <- paste0(fname, ".html", sep="")
    suppressMessages(reactablefmtr::save_reactable_test(checkseq_react, fname_html))
    openFileInOS(fname_html)
    message(cat(paste0("\033[97;", 40, "m","HTML file with [PASS + FAIL] sequences exported: ", "\033[0m", 
                       "\033[0;", 32, "m", " ", file.path(path, "output", unlist(strsplit(fname_html, '/'))[length(unlist(strsplit(fname_html, '/')))]),"\033[0m", "\n")))
  }

  #export CSV files----------------------------------------------------------------------
  if(export_csv==TRUE) {
    #PASS sequences
    fname_csv_pass <- paste0(fname_pass, ".csv", sep="")
    if(isEmpty(checkseq.sub.pass[1])){
      message(cat(paste0("\033[97;", 40, "m","CSV file with [PASS] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", "No sequences passed, file not exported.","\033[0m")))
      } else {
      write.csv(file=fname_csv_pass,checkseq.sub.pass)
      message(cat(paste0("\033[97;", 40, "m","CSV file with [PASS] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, "output", unlist(strsplit(fname_csv_pass, '/'))[length(unlist(strsplit(fname_csv_pass, '/')))]),"\033[0m",
                         "\033[0;", 31, "m", "  <--- Required in Step 2: 'assign_taxonomy'","\033[0m")))
    }
    #FAIL sequences
    fname_csv_fail <- paste0(fname_fail, ".csv", sep="")
    if(isEmpty(checkseq.sub.fail[1])){
      message(cat(paste0("\033[97;", 40, "m","CSV file with [FAIL] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", "No sequences failed, file not exported.","\033[0m", "\n")))
    } else {
      write.csv(file=fname_csv_fail, checkseq.sub.fail)
      message(cat(paste0("\033[97;", 40, "m","CSV file with [FAIL] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, "output", unlist(strsplit(fname_csv_fail, '/'))[length(unlist(strsplit(fname_csv_fail, '/')))]),"\033[0m", "\n")))
    }
  }

  #export trimmed sequences in FASTA file-------------------------------------------------
  if(export_fasta == TRUE) {
    #PASS sequences
    #---------------
    fname_fasta_pass <- paste0(fname_pass, ".fasta", sep="")
    if(isEmpty(checkseq.sub.pass[1])){
      message(cat(paste0("\033[97;", 40, "m","FASTA file with [PASS] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", "No sequences failed, file not exported.","\033[0m")))
    } else {
      trimlist_pass <- as.list((unlist(trimlist))[which(!checkseq.sub$filename %in% seq.warnings)])
      abif_files_pass <- abif_files[!abif_files %in% seq.warnings]
    write.fasta(sequences=trimlist_pass, names=abif_files_pass, as.string=TRUE, nbchar = 1000,
                file.out=fname_fasta_pass)
    message(cat(paste0("\033[97;", 40, "m","FASTA file with [PASS] sequences exported: ", "\033[0m",
                       "\033[0;", 32, "m", " ", file.path(path, "output", unlist(strsplit(fname_fasta_pass, '/'))[length(unlist(strsplit(fname_fasta_pass, '/')))]),"\033[0m")))
    }
    #FAIL sequences
    #---------------
    fname_fasta_fail <- paste0(fname_fail, ".fasta", sep="")
    if(isEmpty(checkseq.sub.fail[1])){
      message(cat(paste0("\033[97;", 40, "m","FASTA file with [FAIL] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", "No sequences failed, file not exported.","\033[0m")))
    } else {
      trimlist_fail <- as.list((unlist(trimlist))[which(checkseq.sub$filename %in% seq.warnings)])
      abif_files_fail <- abif_files[abif_files %in% seq.warnings]
      write.fasta(sequences=trimlist_fail, names=abif_files_fail, as.string=TRUE, nbchar = 1000,
                  file.out=fname_fasta_fail)
      message(cat(paste0("\033[97;", 40, "m","FASTA file with [FAIL] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, "output", unlist(strsplit(fname_fasta_fail, '/'))[length(unlist(strsplit(fname_fasta_fail, '/')))]),"\033[0m", "\n")))
    }
  }
  if(export_fasta_revcomp==TRUE) {
    #PASS sequences
    #---------------
    fname_fasta_revcomp_pass <- paste0(fname_pass, "_revcomp.fasta", sep="")
    if(isEmpty(checkseq.sub.pass[1])){
      message(cat(paste0("\033[97;", 40, "m","FASTA (reverse complement) file with [PASS] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", "No sequences failed, file not exported.","\033[0m")))
    } else {
      trimlist_pass <- as.list(paste(Biostrings::reverseComplement(DNAStringSet((unlist(trimlist))[which(!checkseq.sub$filename %in% seq.warnings)]))))
      abif_files_pass <- abif_files[!abif_files %in% seq.warnings]
      suppressWarnings({
        write.fasta(sequences=trimlist_pass, names=abif_files_pass, as.string=TRUE, nbchar = 1000,
                    file.out=fname_fasta_revcomp_pass)
      })
      message(cat(paste0("\033[97;", 40, "m","FASTA (reverse complement) file with [PASS] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, "output", unlist(strsplit(fname_fasta_revcomp_pass, '/'))[length(unlist(strsplit(fname_fasta_revcomp_pass, '/')))]),"\033[0m")))
    }
    #FAIL sequences
    #---------------
      fname_fasta_revcomp_fail <- paste0(fname_fail, "_revcomp.fasta", sep="")
      if(isEmpty(checkseq.sub.pass[1])){
        message(cat(paste0("\033[97;", 40, "m","FASTA (reverse complement) file with [FAIL] sequences exported: ", "\033[0m",
                           "\033[0;", 32, "m", " ", "No sequences failed, file not exported.","\033[0m")))
      } else {
        trimlist_fail <- as.list(paste(Biostrings::reverseComplement(DNAStringSet((unlist(trimlist))[which(checkseq.sub$filename %in% seq.warnings)]))))
        abif_files_fail <- abif_files[abif_files %in% seq.warnings]
        suppressWarnings({
          write.fasta(sequences=trimlist_fail, names=abif_files_fail, as.string=TRUE, nbchar = 1000,
                      file.out=fname_fasta_revcomp_fail)
        })
        message(cat(paste0("\033[97;", 40, "m","FASTA (reverse complement) file with [FAIL] sequences exported: ", "\033[0m",
                           "\033[0;", 32, "m", " ", file.path(path, "output", unlist(strsplit(fname_fasta_revcomp_fail, '/'))[length(unlist(strsplit(fname_fasta_revcomp_fail, '/')))]),"\033[0m")))
    }
  }
}
