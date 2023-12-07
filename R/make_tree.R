#' @name make_tree
#' @title Generate a phylogenetic tree from an isoLIB output file
#' @description This script will help the user make a simple phylogenetic tree from a strain library.
#' It will allow the user to colour the tree by taxonomic rank only.
#' See \link[ggtree]{ggtree} documentation for more information on customization options available.
#' @export
#' @param input Full path to isoLIB strain library output file in .CSV format.
#' @seealso \code{\link{isoLIB}}
#' @return Returns a \link[ggtree]{ggtree} class object
#' @importFrom stringr str_subset
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_squish
#' @importFrom ape as.DNAbin
#' @importFrom ape dist.dna
#' @importFrom ape njs
#' @importFrom treeio as_tibble
#' @importFrom treeio tibble
#' @importFrom full_join
#' @importFrom as.treedata
#' @importFrom scales hue_pal
#' @importFrom utils select.list
#' @importFrom ggtree ggtree
#' @importFrom ggtree geom_tippoint
#' @importFrom ggtree scale_colour_manual
#' @importFrom ggtree geom_tiplab
#' @importFrom msa msaMuscle
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

#' #Step 4: Make a tree from isoLIB output CSV file
#' fpath4 <- file.path(fpath1, "isolateR_output/03_isoLIB_results.csv")
#' make_tree(input= fpath4)


make_tree <- function(input = NULL){
  if(is.null(input)){stop("Condensed library file required.", call.=FALSE)}


  #get paths
  clean_input <- stringr::str_replace_all(input, '\\\\', '/')
  workingdir <- paste(unlist(strsplit(clean_input, '/'))[1:(length(unlist(strsplit(clean_input, '/')))-1)],collapse="/")
  proj_name <- unlist(strsplit(workingdir, '/'))[(length(unlist(strsplit(workingdir, '/')))-1)]
  setwd(workingdir)

  lib_file <- read.csv(clean_input)

  #clean library file, remove spaces from strain names
  lib_file$filename <- stringr::str_squish(lib_file$filename)
  lib_file$filename <- stringr::str_replace_all(lib_file$filename," ", "_")
  colnames(lib_file) <- stringr::str_replace_all(colnames(lib_file), "rank_", "")



  #--------------------------
  #--------------------------Do MSA on strain sequences


  seqset <- Biostrings::DNAStringSet(lib_file$seqs_trim)
  names(seqset) <- lib_file$filename

  #MSA alginment with muscle using N-J method, default params otherwise
  #params: -in noFile -cluster neighborjoining -gapOpen -400 -gapExtend -0 -maxiters 16 -seqtype dna -spn  -clwstrict
  #MUSCLE 3.8.31
  alignedSet <- msa::msaMuscle(seqset, cluster = "neighborjoining", verbose = TRUE)

  alignedBin <- ape::as.DNAbin(alignedSet)
  ddist <- ape::dist.dna(x=alignedBin)
  phy_tree <- ape::njs(ddist)

  #--------------------------
  #--------------------------Add taxonomic data to the tree

  tree.tab <- treeio::as_tibble(phy_tree)

  tax <- treeio::tibble(label = lib_file$filename,
                        P = lib_file$phylum,
                        C = lib_file$class,
                        O = lib_file$order,
                        Fa = lib_file$family,
                        G = lib_file$genus,
                        S = lib_file$species)
  tree.tab.tax <- treeio::full_join(tree.tab,tax, by = 'label')

  #convert back to tree object with taxonomy
  tree.tax <- treeio::as.treedata(tree.tab.tax)

  phyla_tips <- split(tax$label,tax$P)
  class_tips <- split(tax$label,tax$C)
  ord_tips <- split(tax$label,tax$O)
  fam_tips <- split(tax$label,tax$Fa)
  gen_tips <- split(tax$label,tax$G)
  sp_tips <- split(tax$label,tax$S)

  tree.tax.gr <- ggtree::groupOTU(tree.tax, phyla_tips, "Phylum")
  tree.tax.gr <- ggtree::groupOTU(tree.tax.gr, class_tips, "Class")
  tree.tax.gr <- ggtree::groupOTU(tree.tax.gr, ord_tips, "Order")
  tree.tax.gr <- ggtree::groupOTU(tree.tax.gr, fam_tips, "Family")
  tree.tax.gr <- ggtree::groupOTU(tree.tax.gr, gen_tips, "Genus")
  tree.tax.gr <- ggtree::groupOTU(tree.tax.gr, sp_tips, "Species")

  tax.rank <- "Order" #set default tax rank
  colours <- scales::hue_pal()(length(names(ord_tips))) #set default colour palette


  #quick function to determine new length of tax.rank if the user decides to change it
  #will allow new colour palette to be created if user changes tax rank, and as a stopper
  #if the user tries to change colours with less than the required amount of colours for the rank their tree is in
  taxlen <- function(x){
    #where x is tax.rank
    if(x == "Phylum") y <- length(names(phyla_tips))
    if(x == "Class") y <- length(names(class_tips))
    if(x == "Order") y <- length(names(ord_tips))
    if(x == "Family") y <- length(names(fam_tips))
    if(x == "Genus") y <- length(names(gen_tips))
    if(x == "Species") y <- length(names(sp_tips))
    return(y)
  }


  #--------------------------
  #--------------------------View the default tree and specify an outgroup if desired

  viewtree <- ggtree::ggtree(tree.tax.gr,
                             layout = "rectangular",
                             branch.length='none',
                             ggtree::aes(color= eval(parse(text = tax.rank)) ))+
    ggtree::geom_tippoint(size=1) +
    ggtree::geom_tiplab(size=3) +
    ggtree::scale_colour_manual(name=tax.rank,
                                values = colours,
                                na.translate=FALSE)

  #See if the user would like to set a root for their tree
  message("----Determining outgroup:")
  while(TRUE){
    message("Check the plot for a current view of the tree.You may need to expand the plot to view labels.")
    plot(viewtree)
    choice <- utils::select.list(title = "Would you like to select a different root?", choices = c("No","Yes"), graphics = FALSE, multiple =FALSE)

    if(choice == "No"){
      break
    } else {

      message("Type the strain ID as it appears in the tree that you would like to make the outgroup.")
      newroot <- readline(prompt="Enter strain ID: ")
      tree.tax <- ape::root(tree.tax, outgroup = newroot)

      #need to regroup after changing outgroup
      tree.tax.gr <- ggtree::groupOTU(tree.tax, phyla_tips, "Phylum")
      tree.tax.gr <- ggtree::groupOTU(tree.tax.gr, class_tips, "Class")
      tree.tax.gr <- ggtree::groupOTU(tree.tax.gr, ord_tips, "Order")
      tree.tax.gr <- ggtree::groupOTU(tree.tax.gr, fam_tips, "Family")
      tree.tax.gr <- ggtree::groupOTU(tree.tax.gr, gen_tips, "Genus")
      tree.tax.gr <- ggtree::groupOTU(tree.tax.gr, sp_tips, "Species")

      viewtree <- ggtree::ggtree(tree.tax.gr,
                                 layout = "rectangular",
                                 branch.length='none',
                                 ggtree::aes(color= eval(parse(text = tax.rank)) ))+
        ggtree::geom_tippoint(size=1)+
        ggtree::geom_tiplab(size=3)+
        ggtree::scale_colour_manual(name=tax.rank,
                                    values = colours,
                                    na.translate=FALSE)
    }
  }


  #--------------------------
  #--------------------------Specify taxonomic rank and colours for the tree

  #Now to make the tree customized by desired tax rank and colours
  message("----Select taxonomic rank and colours:")
  while(TRUE){
    plot(viewtree)
    choice <- select.list(title = "Would you like to make changes to the colours or taxonomic level to colour the tree by?", choices = c("Colours","Taxonomic level","No"), graphics = FALSE, multiple =FALSE)
    if(choice == "Colours"){
      message("To supply custom colours for your tree you MUST include at least as many as different taxonomic labels there are. For instance if you have 4 phyla in your tree and you would like to colour by phylum, input at least 4 colours.")
      message("Input format example: #00A9FF, #C77CFF, #FF61CC")
      colours <- readline(prompt="Enter a list of hex colours in a comma seperate list: ")
      colours <- stringr::str_replace_all(colours," ","")
      colours <- strsplit(colours,",")
      colours <- colours[[1]]

      if(length(colours) >= taxlen(tax.rank)){
        viewtree <- ggtree(tree.tax.gr, layout = "rectangular", branch.length='none', aes(color= eval(parse(text = tax.rank)) ))+ geom_tippoint(size=1)+ geom_tiplab(size=3)+scale_colour_manual(name=tax.rank, values = colours, na.translate=FALSE)
      } else {
        message(paste0("Insufficient colours supplied, user gave ",length(colours)," but ", taxlen(tax.rank)," is needed for the current taxonomic rank of the tree. Try again."))
      }


    } else if(choice == "Taxonomic level"){
      tax.choice <- select.list(title = "Which taxonomic rank should be used?", choices = c("Phylum","Class","Order","Family","Genus","Species"), graphics = FALSE, multiple =FALSE)
      tax.rank <- tax.choice

      colours <- scales::hue_pal()(taxlen(tax.rank)+1)

      viewtree <- ggtree::ggtree(tree.tax.gr,
                                 layout = "rectangular",
                                 branch.length='none',
                                 ggtree::aes(color= eval(parse(text = tax.rank)) ))+
        ggtree::geom_tippoint(size=1)+
        ggtree::geom_tiplab(size=3)+
        ggtree::scale_colour_manual(name=tax.rank,
                                    values = colours,
                                    na.translate=TRUE)

    } else {
      break
    }

  }

  #--------------------------
  #--------------------------Change xlim if required and output file to get good dimensions

  #To deal with any remaining formatting
  message("----select formatting:")
  message("If you have long tree tip labels (strain names) they may be cut abruptly in the current plot.")
  message("The following prompts will guide you in choosing plot limits and output file size to make an adequate tree. Make sure you open the output file when prompted, and select EXIT when you are happy with the results.")
  message("The output image of your tree will be in the directory where your strain library is.")

  xlim.max <- NA
  x <- 30
  y <- 30

  while(TRUE){
    #cladogram no labels or rect with labels and chose xlim if long lable names
    plot(viewtree)
    choice <- select.list(title = "Would you like to make changes to the tree plot sizing?", choices = c("Yes","No"), graphics = FALSE, multiple =FALSE)
    if(choice == "No"){
      #do nothing
    } else {
      message("You can reformat the tree multiple times, recommened values for xlim are 10-100, higher the value for many labels that are very long.")
      xlim.max <- as.numeric(readline(prompt="Input an coordinate to limit the tree plot to: "))

      viewtree <- ggtree(tree.tax.gr, layout = "rectangular", branch.length='none', aes(color= eval(parse(text = tax.rank)) ))+ geom_tippoint(size=1)+ geom_tiplab(size=3,offset = 0.05)+scale_colour_manual(name=tax.rank, values = colours, na.translate=FALSE)+xlim(0,xlim.max)

      plot(viewtree)
    }

    ggsave(paste0('library_tree____',proj_name,'.png'), units = "cm", height = y, width = x, dpi = 300)

    message(paste0("Check the output file: ",'library_tree____',proj_name,'.png'))

    while(TRUE){
      choice.out <- select.list(title = "Based on the .png would you like to change the plot limits again, change the output file dimensions, or are you happy with the current output?", choices = c("Change the plot limits","Change the output file", "Exit script"), graphics = FALSE, multiple =FALSE)
      if(choice.out == "Change the plot limits"){
        break
      } else if(choice.out == "Exit script"){
        stop("User chose to exit.", call.=FALSE)
      } else {
        message("Output dimensions are in cm, do not include units in your following inputs.")
        x <- as.numeric(readline(prompt="Enter desired width of the output image (default 30): "))
        y <- as.numeric(readline(prompt="Enter desired height of the output image (default 30): "))
        ggsave(paste0('library_tree____',proj_name,'.png'), units = "cm", height = y, width = x, dpi = 300)

        message(paste0("Check the output file: ",'library_tree____',proj_name,'.png'))
      }
    }
  }
  return(viewtree)
}
