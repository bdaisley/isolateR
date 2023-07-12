
#abif_fasta2 function
#------------------------------------------------------------------------------------
abif_fasta2 <- function(folder=NULL,
                        export_html=TRUE,
                        export_fasta=TRUE,
                        reversecomp=TRUE,
                        verbose=TRUE,
                        exclude=NULL,
                        output='V3-V6seq.FASTA') {
  
  # function requirements------------------------------------------------------------
  #checking for required packages; installing those not yet installed
  if(require(seqinr)==FALSE) install.packages('seqinr')
  if(require(stringr)==FALSE) install.packages('stringr')
  if(require(sangerseqR)==FALSE){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("sangerseqR")
  }
  if(require(sangeranalyseR)==FALSE){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("sangeranalyseR")
  }
  if(require(reactable)==FALSE) install.packages('reactable')
  if(require(reactablefmtr)==FALSE) install.packages('reactablefmtr')
  if(require(pander)==FALSE) install.packages('pander')
  if(require(dataui)==FALSE){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("remotes")
    remotes::install_github("timelyportfolio/dataui")
  }
  
  #loading required packages
  library(seqinr)
  library(stringr)
  library(sangerseqR)
  library(sangeranalyseR)
  library(reactable)
  library(reactablefmtr)
  library(pander)
  library(dataui)
  
  # Reading files--------------------------------------------------------------------
  # Check folder path
  if(is.null(folder)) stop('Name of folder not supplied. Please check the execution script preamble to make sure it is entered correctly', call.=FALSE)
  
  # Correct any syntax issues with folder path
  path <- str_replace_all(folder, '\\\\', '/')
  
  #Load abif files
  fname <- list.files(path)
  pattern <- 'ab1'
  abif_files <- str_subset(fname, pattern)
  
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
    
    abif1 <- SangerRead(readFeature           = "Forward Read",
                        readFileName          = fpath,
                        geneticCode           = GENETIC_CODE,
                        TrimmingMethod        = "M2",
                        M1TrimmingCutoff      = NULL,
                        M2CutoffQualityScore  = 45,
                        M2SlidingWindowSize   = 5,
                        baseNumPerRow         = 100,
                        heightPerRow          = 200,
                        signalRatioCutoff     = 0.5)
    #qualityBasePlot(sangerReadF)
    
    start=abif1@QualityReport@trimmedStartPos
    abif1.start = NULL
    abif1.start <- abif1@QualityReport@trimmedStartPos
    if(is.null(abif1.start) == TRUE){abif.start <- 1}
    abif1.end <- nchar(paste(abif1@primarySeq))
    
    abif2 <- DNAStringSet(abif1@primarySeq, start=abif1.start, end=abif1.end)
    #abif2 <- Biostrings::reverseComplement(DNAStringSet(paste(abif1@primarySeq),
    #                                                    start=abif1.start,
    #                                                    end=abif1.end))
    
    abif3 <- DECIPHER::TrimDNA(DNAStringSet(paste(abif2)),
                               leftPatterns="",
                               rightPatterns="CTGCTGCCTYCCGTA",
                               minWidth=1,
                               maxDistance = 0.2,
                               minOverlap = 15,
                               allowInternal=TRUE,
                               threshold = 1,
                               maxAverageError = 1,
                               maxAmbiguities = 0.2,
                               type="both")
    
    abif3.end <- end(abif3[[1]]) + nchar("CTGCTGCCTYCCGTA")
    if(abif3.end == nchar(paste(abif2)) + nchar("CTGCTGCCTYCCGTA")){abif3.end <- abif1@QualityReport@trimmedFinishPos}
    
    
    abif <- DNAStringSet(abif2, start=1, end = abif3.end) # equivalent to line below
    #abif <- DNAStringSet(abif1@primarySeq, start=abif1.start, end = (abif1.start + abif3.end -1))
    abif.scores.xx <- abif1@QualityReport@qualityPhredScores[abif1.start:(abif1.start + abif3.end -1)]
    abif.scores <- c(rep(0, abif1.start), abif.scores.xx, rep(0.1, 810-(abif1.start + abif3.end -1)))
    
    trim.start <- c(trim.start,abif1.start)
    trim.end <- c(trim.end,(abif1.start + abif3.end -1))
    
    rawseq <- as.character(paste(abif))
    
    if(verbose==TRUE){
      #print a message to mark progress
      msg <- sprintf('Processing "%s"', abif_files[i])
      message(msg)
    }
    
    #:::::::::::::::::::::::::::::::
    # Save seq + quality scores
    #:::::::::::::::::::::::::::::::
    rawsparklist[[i]] <- c(abif1@QualityReport@qualityPhredScores, rep(0.1, 810-length(abif1@QualityReport@qualityPhredScores)))
    trimsparklist[[i]] <- abif.scores
    rawlist[[i]] <- paste(abif1@primarySeq)
    trimlist[[i]] <- paste(abif)
    #rawlist[[i]] <- paste(Biostrings::reverseComplement(DNAStringSet(paste(abif1@primarySeq))))
    #trimlist[[i]] <- paste(Biostrings::reverseComplement(DNAStringSet(abif)))
    
    # building trim check for export-------------------------------------------------
    entry <- data.frame(filename=abif_files[i],
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
      msg <- sprintf('The file "%s" has a length of <200 bp after trimming. Manual check of sequence is recommended.', abif_files[i])
      #message(msg)
      warning(msg, call.=FALSE)
      #
      #  seq_pass <- FALSE
      #  failreason <- paste(failreason, msg)
    }
    # building trim check for export-------------------------------------------------
    checkseq <- rbind(checkseq, entry)
    #ntseq <- trimlist[[i]]
  }
  # end of file loop
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
    select(filename, raw_seq, phred_raw, Ns_raw, length_raw, spark_raw, arrcol,
           trim_seq, phred_trim, Ns_trim, length_trim) #%>%
  #`colnames<-`(c("filename", "rawSeq", "rawQ", "rawN", "rawLength", "rawSpark", "->",
  #               "trmSeq", "trmQ", "trmN", "trmLength", "trmSpark"))
  
  checkseq.sub$spark_raw <- lapply(checkseq.sub$spark_raw, function(x) list(x))
  
  #trim.start <- c(rep(1:40))
  #trim.end <- c(rep(500:550))
  #abif1.df.spark <- sparkline(abif1@QualityReport@qualityPhredScores, width=200, fillColor=FALSE)
  
  checkseq_react <- reactable(checkseq.sub,
                              fullWidth=FALSE, searchable = TRUE, bordered = TRUE, resizable =TRUE, #width = 1200,
                              defaultPageSize=1000, highlight = TRUE, showSortable = TRUE, compact=TRUE, wrap = TRUE, #highlight=TRUE,
                              #theme = reactableTheme(headerStyle = list(wrap=TRUE)),
                              theme = reactableTheme(                  #borderColor = "#555",
                                headerStyle = list(
                                  "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                                  "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)"))),
                              #borderColor = "#555")),
                              #defaultColDef = colDef(minWidth = 100), #
                              #defaultColDef = colDef(style = reactablefmtr::cell_style(checkseq.sub, font_size=12)),
                              style = list(minWidth = 1400),
                              defaultColDef = colDef(style = "white-space: nowrap; font-size: 16px; font-family: Courier New", minWidth= 75, align="center"),
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
                                raw_seq = colDef(name = "rawSeq"),
                                trim_seq = colDef(name = "trimSeq"),
                                arrcol = colDef(name = ""),
                                filename = colDef(style = "white-space: nowrap;", minWidth= 150),
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
                                                                                                      round_edges = FALSE, align_bars="left")) ))
  #phred_trim = colDef(data_bars(checkseq.sub$phred_trim, text_position = "none", fill_color="#6f86ab"))
  #phred_raw = colDef(minWidth= 50, headerStyle = list(wrap=FALSE)),
  #phred_trim = colDef(minWidth= 50, headerStyle = list(wrap=FALSE)),
  #phred_raw = colDef(style = "white-space: nowrap;", minWidth= 50),
  #phred_trim = colDef(style = "white-space: nowrap;", minWidth= 50),
  #length_raw = colDef(style = "white-space: nowrap;", minWidth= 50),
  #length_trim = colDef(style = "white-space: nowrap;", minWidth= 50),
  #Ns_raw = colDef(style = "white-space: nowrap;", minWidth= 50),
  #Ns_trim = colDef(style = "white-space: nowrap;", minWidth= 50),
  #raw_seq = colDef(style = "white-space: nowrap;", minWidth= 50),
  #trim_seq = colDef(style = "white-space: nowrap;", minWidth= 50)
  #phred_raw = colDef(style = "white-space: nowrap;", maxWidth= 500)
  #phred_raw = colDef(style = "white-space: nowrap;", maxWidth= 500)
  #phred_raw = colDef(style = "white-space: nowrap;", maxWidth= 500)
  #filename = colDef(maxWidth= 500)))
  #raw_seq = colDef(style = "
  ##                  white-space: nowrap;
  #                  overflow: hidden;
  #                  text-overflow: ellipsis;
  #                  hover: visible"),
  
  
  reactablefmtr::save_reactable_test(checkseq_react, file.path(path, 'V3-V6seq_check.html'))
  
  openFileInOS(file.path(path, 'V3-V6seq_check.html'))
  
  
  
  #building output file--------------------------------------------------------------
  fname_out <- file.path(path, output)
  
  if(export_html == TRUE) {
    fexport <- str_extract(output, '.*(?=\\.)')
    fexport <- paste0(fexport, '_check.csv')
    fexport <- file.path(path, fexport)
    write.csv(file=fexport, checkseq.sub %>% select(-spark_raw))
  }
  
  if(export_html==TRUE) {
    write.csv(file=file.path(path, 'V3-V6seq_check.csv'), checkseq.sub %>% select(-spark_raw))
  }
  
  #export as one FASTA file
  if(export_fasta == TRUE) {
    write.fasta(sequences=trimlist, names=abif_files, as.string=TRUE, nbchar = 1000,
                file.out=fname_out)
  }
  
  #write.csv(file="fails.csv",failseq)
  
}
