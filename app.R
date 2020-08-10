# Code created by Matthew Elliott
# Contact: melliot1@ucsc.edu (valid until 2023)
# Phone:  231-392-1263   (just in case)


  #########################
  #########################
  #  Intro Text
  #########################
  #########################

  # I should rewrite this intro
  # This was created for me, not others to read
    
  # Informative links:
  # * https://reneshbedre.github.io/blog/expression_units.html
  # * https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
  
  # To Do:
  #     clean TSNE Plot

  
  
  #########################
  
  #########################
  #########################
  #  Head
  #########################
  #########################
  

  # Set up R Shiny
  library(shiny)
  library(shinyjs)  #library(shinyWidgets)
  library(plotly)
  library(ggthemes)
  library(ggplot2)
  library(hash)
  library(DT)
  library(RColorBrewer)
  
  # Libraries for combat seq
  library("limma")
  library("sva")
  library("edgeR")
  library("Rtsne")
  library("GISTools")
  library("plyr")
  
  # Combat Seq Packages:
  addResourcePath("www", paste(getwd() , "/www", sep="") ) # Have shiny recognize folder with files
  source("www/helper_seq.R") # Helper functions for Combat Seq
  
  
  #########################
  
  
  #########################
  #  UI
  ######################### 
  
  ui = fluidPage( style='margin-left:5px; margin-right:5px', title="COVID-19 County Tracker", div( class="container-fluid", 
                  
          ######################################################################
          #  Head HTML
          ###################################################################### 
          
          shinyjs::useShinyjs(),  # ShinyJS: https://deanattali.com/shinyjs/overview
          div( class="hidden", selectInput('page', choices = c('home','data','run','results'), selected = 'home', label=NULL  )),
          p(tags$small(" .")), # vertical white space
          h2( class="text-center", "HarmonyRNA"), 
          fluidRow(div( class="col-xs-12 col-md-offset-2 col-md-8 col-lg-offset-3 col-lg-6", HTML("<hr style='padding:0px; margin:0px'>") )),
          tags$head(tags$style(".shiny-notification  {position: fixed; top: 20% !important;left: 30% !important; width: 40% !important; font-size: 20px }")),
          
          ######################################################################
          #  Home Screen
          ######################################################################
          
          conditionalPanel( condition = "input.page == 'home'",
               # Introductory Text
               h4(class="text-center", "Harmonize RNA-Seq Count and TPM Files" ),
               p( class="text-center", style="font-size:1.1em", 
                 "RNA-Seq datasets are often combined together to run expirements across different batches or to perform metanalysis studies.
                 Correctly harmonizing RNA-Seq data requires specific statistical techniques and takes a considerable amount of work to code.
                 For this reason the", a(href="http://buttelab.ucsf.edu/", "Butte Lab"), #,inline = T, target = "_blank"),
                 "has developed", strong("HarmonyRNA,"), " an easy to use web based RNA-Seq harmonization pipeline.
                 Watch the tutorial video below. Then harmonize your own RNA-Seq data."
                 ),
               
               # Start button
               p(tags$small(" .")), # vertical white space
               h4( class="text-center", "Did you watch the tutorial?" ),
               div(class="text-center",  actionButton("startBtn", label="Start", class="btn btn-danger text-center") ), #, style="width: 100%") ),
               
               # Video Tutorial
               HTML( "<div class='row'><div class='col-xs-offset-0 col-xs-12 col-sm-offset-2 col-sm-8 col-lg-offset-3 col-lg-6'>
                          <div class='embed-responsive embed-responsive-16by9' >
                            <iframe class='embed-responsive-item'
                            src='https://www.youtube.com/embed/lm3t6yaIlV8' style='padding: 10px;'></iframe></div>
                            <hr><hr>
                          </div></div>"),
                   
               # Footer Text
               div( class="col-xs-12 col-md-7",
                 p( strong("Cite Us: "), "[Placeholder, TBD] Zalocusky KA, Kan MJ, Hu Z, Dunn P, Thomson E, Wiser J, Bhattacharya S, Butte AJ. The 10,000 Immunomes Project: Building a Resource for Human Immunology. Cell reports. 2019 Oct 9;25(2):513-22. PMID:30304689"  ),   
                 p( strong("Cite Combat-Seq: "), "ComBat-Seq: batch effect adjustment for RNA-Seq count data Yuqing Zhang, Giovanni Parmigiani, W. Evan Johnson bioRxiv 2020.01.13.904730; doi: https://doi.org/10.1101/2020.01.13.904730"  ),  
               ),
               div( class="col-xs-12 col-md-5",
                    p( strong("Data:"),  a( href="www/harmony_example.zip", "Download data"), " used in the video tutorial. Data comes from",
                       a( href="https://www.immport.org/shared/search", "ImmPort,"), "a public repository for immunology", tags$small("(Studies: SDY1412, SDY1172, & SDY1092)"),". It was originally harmonized for ",
                       a( href="https://10kimmunomes.ucsf.edu", "10k Immunomes.") ),
                    p( strong("Code:"), "Our code is open and publicly available on", a( href="https://hub.docker.com/r/pupster90/combat-seq", "Github"),
                       "and ",a(href="https://hub.docker.com/r/pupster90/combat-seq", "Dockerhub."),
                       "Both sites provide an in-depth tutorial on how to get started." ),
               ),
               #fluidRow( p("asdf asdf ") ),
               HTML("<div class='row'><div class=col-xs-12><p>.</p></div></div>"),
               div( class="col-xs-12 col-sm-offset-2 col-sm-8 col-md-offset-3 col-md-6 col-lg-offset-4 col-lg-4",
                  img(src='www/bakar_logo.png', class="img-responsive")        
               ) 
          ), # end conditional pannel
  
          ######################################################################
          #  Load Data  
          ######################################################################
          
          # Useful links for File upload:
          # https://shiny.rstudio.com/articles/upload.html
          conditionalPanel( condition = "input.page == 'data'",
                
                fluidRow( div( class="col-xs-offset-0 col-xs-12 col-sm-offset-2 col-sm-8",
                      p("."), # <- white space
                      p( style="font-size:1.1em", "Click the ", tags$b("Browse"), "button to upload a TPM or Counts file. A table of the file should then appear. Use the ", tags$b("Options"), " section to get the dataset into the correct format. When ready, click ", tags$b("Harmonize"), " to start the harmonization process.")
                    ) ), 
                
                # View uploaded datasets
                p("."), #<- whitespace
                p( strong("Datasets: "),  textOutput("og_data", inline=TRUE )  ),
                
                # Upload a dataset 
            div( id="chooseFile",    
                h3("Choose a File"),
                div( class="col-xs-10 col-sm-6 col-md-4 col-lg-3",
                  fileInput("file1", NULL, multiple = FALSE ) #, accept = c("text/csv","text/comma-separated-values,text/plain",".csv") )
                  ),
  
                # Options Button
                fluidRow( HTML('<div class="row">
                        <button type="button" class="btn btn-default btn-sm" data-toggle="collapse" data-target="#demo">Options</button>
                     </div>') ), 
                #<div class="row"><div class="col-xs-offset-2"> </div>
                
                # Warnings:
                hidden( p( id="warning_first_column", style="color:red", "The first column should be the ID's used to match datasets.")),
                hidden( p( id="warning_duplicates", style="color:red", "Some ID's appear multiple times. Duplicates will be removed." )),
                hidden( p( id="warning_numbers", style="color:red", "All columns after the Gene ID's should contain only numbers. Cannot Procceed.")),
                #hidden( p( id="warning_size", style="color:red", "Files can be at most 30MB large.")),
                
                div( class="col-xs-offset-2",
                     disabled( actionButton("dataSaveBtn", label="Next", class="btn btn-success") )
                ),
                
                # Customization Options
                fluidRow(
                div( id="demo", class="collapse",
                     fluidRow( column(12, tags$hr() )),
                     fluidRow(
                       div( class="col-sm-3 col-xs-5", radioButtons("sep", "Separator",
                                                                    choices = c(Comma = ",",
                                                                                Semicolon = ";",
                                                                                Tab = "\t"),
                                                                    selected = ",") ),
                       div( class="col-sm-3 col-xs-5", radioButtons("quote", "Quote",
                                                                    choices = c(None = "",
                                                                                "Double Quote" = '"',
                                                                                "Single Quote" = "'"),
                                                                    selected = '"') ),
                       div( class="col-sm-2 col-xs-6", checkboxInput("header", "Header", TRUE) )
                     ),
                     tags$hr()
                )
                ), # end "collapse" section 
                  
                # Output Table        #   Overflow scrolling: https://github.com/rstudio/shinydashboard/issues/40  
                fluidRow(div(  style = 'overflow-x: scroll', dataTableOutput("contents") ))
                
            ), # End of chooseFile
                
            hidden( div( id="goHarmonize",
                 fluidRow(
                   column(2, actionButton("addDataBtn", label="Add Data", class="btn btn-success", style="margin-top: 5px; margin-bottom: 10px")   ),
                   column(2, hidden( actionButton("harmonizeBtn", label="Harmonize", class="btn btn-danger", style="margin-top: 5px; margin-bottom: 10px") )   )
                 )
            ) ), # end hidden
            
            hidden( div( id="TPMorCounts",
                         fluidRow(h3("Is your data counts or TPM's?")),
                         fluidRow(
                           column(2, actionButton("countsBtn", label="Counts", class="btn btn-primary", style="margin-top: 5px; margin-bottom: 10px")   ),
                           column(2, actionButton("tpmBtn", label="TPM", class="btn btn-primary", style="margin-top: 5px; margin-bottom: 10px") )   
                         )
            ) ) # end hidden
          ),
          

          ######################################################################
          #  Run Algorithm
          ######################################################################
          
          conditionalPanel( condition = "input.page == 'run'",
                 hidden( div( id="combat_text_1",
                    fluidRow( div( class="col-xs-offset-0 col-xs-12 col-sm-offset-2 col-sm-8",
                        p("."), # <- white space
                        p( style="font-size:1.1em", "Harmony RNA-Seq finished harmonizing your datasets. The graphs below give summaries of the data before harmonization, after traditional normalization, and then the final results using Combat-Seq. Use these tools to check the algorithm's validity. Click ", tags$b("View Results"), " to visualize specific genes and to download the final dataset."  )
                        ) ), 
                    div(class="text-center", actionButton("viewResultsBtn", label="View Results", class="text-center btn btn-danger") ),
                    fluidRow( div( class='col-xs-offset-0 col-xs-12 col-sm-offset-2 col-sm-8 col-lg-offset-3 col-lg-6', tags$hr() ) ),
                    fluidRow( h3( class="text-center", "Before Harmonization") )
                  )),
                 div( class="col-xs-12 col-md-6", div( plotOutput("example_raw"), style=" max-width: 1000px; " )   ),
                 div( class="col-xs-12 col-md-6", div( plotOutput("tsne_raw"), style=" max-width: 1000px; " )   ),
                 
                 hidden( div(id="combat_text_3", 
                   fluidRow( div( class='col-xs-offset-0 col-xs-12 col-sm-offset-2 col-sm-8 col-lg-offset-3 col-lg-6', style="margin-top: -20px", tags$hr() )),
                   fluidRow( h3( class="text-center", "After Harmonization") )
                   ) ),
                   div( class="col-xs-12 col-md-6", div( plotOutput("example_harmonized"), style="max-width: 1000px; " ),  ),
                   div( class="col-xs-12 col-md-6", div( plotOutput("tsne_harmonized"), style=" max-width: 1000px; " )    )    
          ),  # end conditionalPanel
  
          ######################################################################
          #  Results
          ######################################################################
          
          conditionalPanel( condition = "input.page == 'results'", div( class="container-fluid",
              fluidRow( div( class="col-xs-offset-0 col-xs-12 col-sm-offset-2 col-sm-8",
                   p("."), # <- white space
                   p( style="font-size:1.1em",  "Below you can search for and select a gene to view results for it. Make sure to click ", tags$b("Download"), " to get the harmonized dataset on your computer. Thanks for using Harmony RNA. For questons and feedback, ", a( href="mailto:melliot1@ucsc.edu" ,"contact us.") ),  
                   p(".")
              ) ),
              fluidRow( div( plotlyOutput("main_plot"), style="border: 1px solid gray; max-width: 1000px; margin: auto; padding-top: 15px;" ) ),
              fluidRow( align="center", selectizeInput( 'select_gene', label= NULL, choices= NULL, options= list(maxOptions=15) ) ),
              div( class="text-center", style="padding-top: 10px;", downloadButton( "download_data", label="Download", class="btn btn-success" ) ),
          )),
                            
          ######################################################################
          #  Footer
          ######################################################################
          fluidRow(class="jumbotron", style="background-color: white", HTML("<h1>&nbsp</h1>") ) # vertical whitespace
  ))
  
  
  #########################
  
  #########################
  #########################
  #  SERVER
  #########################
  #########################
  
  # Define server logic required to draw a histogram
  server <- function(input, output, session) {
  
    ######################################################################
    #  Home
    ######################################################################
    
    # Values used in create plots and update page
    vals <- reactiveValues( rna_data = NULL, og_data = hash(), data_temp = NULL, # How to use hashes: https://stackoverflow.com/questions/7818970/is-there-a-dictionary-functionality-in-r/44570412
                            step = 0, batch=NULL, 
                            shared_genes=NULL, gene_names= NULL, subject_names=NULL,
                            start_raw=FALSE, start_harmonized=FALSE, start_main_plot=FALSE,
                            tsne_raw=FALSE, tsne_harmonized=FALSE, isTPM= FALSE,
                            progress=NULL, progress_step=NULL)
    
    observeEvent( input$startBtn, {
      updateSelectInput(session, "page", selected = 'data')
    })
    
    ######################################################################
    #  Combat Seq Function
    ######################################################################
    
    ComBat_seq <- function(counts, batch, group=NULL, covar_mod=NULL, full_mod=TRUE, 
                           shrink=FALSE, shrink.disp=FALSE, gene.subset.n=NULL ){  
      # Create Progress Bar
      vals$progress$inc( vals$progress_step, detail="Starting Combat-Seq" )
      
      ########  Preparation  ########  
      ## Does not support 1 sample per batch yet
      batch <- as.factor(batch)
      if(any(table(batch)<=1)){
        stop("ComBat-seq doesn't support 1 sample per batch yet")
      }
      
      ## Remove genes with only 0 counts in any batch
      keep_lst <- lapply(levels(batch), function(b){
        which(apply(counts[, batch==b], 1, function(x){!all(x==0)}))
      })
      keep <- Reduce(intersect, keep_lst)
      rm <- setdiff(1:nrow(counts), keep)
      countsOri <- counts
      counts <- counts[keep, ]
      
      # require bioconductor 3.7, edgeR 3.22.1
      dge_obj <- DGEList(counts=counts)
      
      ## Prepare characteristics on batches
      n_batch <- nlevels(batch)  # number of batches
      batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch  
      n_batches <- sapply(batches_ind, length)
      n_sample <- sum(n_batches)
      cat("Found",n_batch,'batches\n')
      
      ## Make design matrix 
      # batch
      batchmod <- model.matrix(~-1+batch)  # colnames: levels(batch)
      # covariate
      group <- as.factor(group)
      if(full_mod & nlevels(group)>1){
        cat("Using full model in ComBat-seq.\n")
        mod <- model.matrix(~group)
      }else{
        cat("Using null model in ComBat-seq.\n")
        mod <- model.matrix(~1, data=as.data.frame(t(counts)))
      }
      # drop intercept in covariate model
      if(!is.null(covar_mod)){
        if(is.data.frame(covar_mod)){
          covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), function(i){model.matrix(~covar_mod[,i])}))
        }
        covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x){all(x==1)})]
      }
      # bind with biological condition of interest
      mod <- cbind(mod, covar_mod)
      # combine
      design <- cbind(batchmod, mod)
      
      ## Check for intercept in covariates, and drop if present
      check <- apply(design, 2, function(x) all(x == 1))
      design <- as.matrix(design[,!check])
      cat("Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')
      
      ## Check if the design is confounded
      if(qr(design)$rank<ncol(design)){
        #if(ncol(design)<=(n_batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
        if(ncol(design)==(n_batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")}
        if(ncol(design)>(n_batch+1)){
          if((qr(design[,-c(1:n_batch)])$rank<ncol(design[,-c(1:n_batch)]))){stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
          }else{stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq")}}
      }
      
      ## Check for missing values in count matrix
      NAs = any(is.na(counts))
      if(NAs){cat(c('Found',sum(is.na(counts)),'Missing Data Values\n'),sep=' ')}
      
      ########  Estimate gene-wise dispersions within each batch  ########
      vals$progress$inc( vals$progress_step, detail="Estimate gene-wise dispersions for batches" )
      ## Estimate common dispersion within each batch as an initial value
      disp_common <- sapply(1:n_batch, function(i){
        if((n_batches[i] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)){ 
          # not enough residual degree of freedom
          return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=NULL, subset=nrow(counts)))
        }else{
          return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=mod[batches_ind[[i]], ], subset=nrow(counts)))
        }
      })
      
      ## Estimate gene-wise dispersion within each batch 
      genewise_disp_lst <- lapply(1:n_batch, function(j){
        if((n_batches[j] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)){
          # not enough residual degrees of freedom - use the common dispersion
          return(rep(disp_common[j], nrow(counts)))
        }else{
          return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], design=mod[batches_ind[[j]], ], 
                                        dispersion=disp_common[j], prior.df=0))
        }
      })
      names(genewise_disp_lst) <- paste0('batch', levels(batch))
      
      ## construct dispersion matrix
      phi_matrix <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
      for(k in 1:n_batch){
        phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], n_batches[k]) 
      }
      
      ########  Estimate parameters from NB GLM  ########
      cat("Fitting the GLM model\n")
      glm_f <- glmFit(dge_obj, design=design, dispersion=phi_matrix, prior.count=1e-4) #no intercept - nonEstimable; compute offset (library sizes) within function
      alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample) #compute intercept as batch-size-weighted average from batches
      new_offset <- t(vec2mat(getOffset(dge_obj), nrow(counts))) +   # original offset - sample (library) size
        vec2mat(alpha_g, ncol(counts))  # new offset - gene background expression # getOffset(dge_obj) is the same as log(dge_obj$samples$lib.size)
      glm_f2 <- glmFit.default(dge_obj$counts, design=design, dispersion=phi_matrix, offset=new_offset, prior.count=1e-4) 
      
      gamma_hat <- glm_f2$coefficients[, 1:n_batch]
      mu_hat <- glm_f2$fitted.values
      phi_hat <- do.call(cbind, genewise_disp_lst)
      
      
      ########  In each batch, compute posterior estimation through Monte-Carlo integration  ########  
      if(shrink){
        cat("Apply shrinkage - computing posterior estimates for parameters\n")
        mcint_fun <- monte_carlo_int_NB
        monte_carlo_res <- lapply(1:n_batch, function(ii){
          if(ii==1){
            mcres <- mcint_fun(dat=counts[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]], 
                               gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)
          }else{
            invisible(capture.output(mcres <- mcint_fun(dat=counts[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]], 
                                                        gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)))
          }
          return(mcres)
        })
        names(monte_carlo_res) <- paste0('batch', levels(batch))
        
        gamma_star_mat <- lapply(monte_carlo_res, function(res){res$gamma_star})
        gamma_star_mat <- do.call(cbind, gamma_star_mat)
        phi_star_mat <- lapply(monte_carlo_res, function(res){res$phi_star})
        phi_star_mat <- do.call(cbind, phi_star_mat)
        
        if(!shrink.disp){
          cat("Apply shrinkage to mean only\n")
          phi_star_mat <- phi_hat
        }
      }else{
        cat("Shrinkage off - using GLM estimates for parameters\n")
        gamma_star_mat <- gamma_hat
        phi_star_mat <- phi_hat
      }
      
      ########  Obtain adjusted batch-free distribution  ########
      mu_star <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
      for(jj in 1:n_batch){
        mu_star[, batches_ind[[jj]]] <- exp(log(mu_hat[, batches_ind[[jj]]])-vec2mat(gamma_star_mat[, jj], n_batches[jj]))
      }
      phi_star <- rowMeans(phi_star_mat)
      
      ########  Adjust the data  ########  
      #cat( paste("n_batch: ",n_batch) )
      adjust_counts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
      for(kk in 1:n_batch){
        cat(paste("Adjusting data for batch ",kk,"of",n_batch,"\n"))
        vals$progress$inc( vals$progress_step, detail=paste("Combat-Seq: Adjusting Data ",kk,"/",n_batch) )
  
        counts_sub <- counts[, batches_ind[[kk]]]
        old_mu <- mu_hat[, batches_ind[[kk]]]
        old_phi <- phi_hat[, kk]
        new_mu <- mu_star[, batches_ind[[kk]]]
        new_phi <- phi_star
        adjust_counts[, batches_ind[[kk]]] <- match_quantiles(counts_sub=counts_sub, 
                                                              old_mu=old_mu, old_phi=old_phi, 
                                                              new_mu=new_mu, new_phi=new_phi)
      }
    
      ## Add back genes with only 0 counts in any batch (so that dimensions won't change)
      adjust_counts_whole <- matrix(NA, nrow=nrow(countsOri), ncol=ncol(countsOri))
      dimnames(adjust_counts_whole) <- dimnames(countsOri)
      adjust_counts_whole[keep, ] <- adjust_counts
      adjust_counts_whole[rm, ] <- countsOri[rm, ]
      return(adjust_counts_whole)
    }
    
    ######################################################################
    #  Load Data
    ######################################################################
    
    output$og_data = renderText({ "None" })  # Initial message for current datasets
    
    # Helper: Read File
    readFile <- reactive({
      data <- read.table(input$file1$datapath, header = input$header, sep = input$sep, quote = input$quote, stringsAsFactors=FALSE)
      data
    })
    
    # Code for uploading file and formatting data
    output$contents <- renderDataTable({
      
      # input$file1 will be NULL initially. After the user selects and uploads a file, head of that data file by default, or all rows if selected, will be shown.
      req(input$file1)
      
      # when reading semicolon separated files, having a comma separator causes `read.csv` to error
      tryCatch(  { 
          df <- readFile()
          vals$data_temp = df
          
          # Warnings:
          if( !all( !duplicated( df[,1] ) ) ){ # duplicates warning
            shinyjs::show("warning_duplicates")
          }else{ shinyjs::hide("warning_duplicates") }
          
          if( !all( sapply(df[-1], is.numeric)) ){ # numbers warning
            shinyjs::show("warning_numbers")
            disable( "dataSaveBtn" ) 
            stop("Data is not Numeric")
          }else{ shinyjs::hide("warning_numbers") }
          
          if( is.numeric(df[1])  ){ # first column numbers warning
            shinyjs::show("warning_first_column")
          }else{ shinyjs::hide("warning_first_column") }
          
          # Output table, enable button
          df[2:dim(df)[2]] = round(df[2:dim(df)[2]], digits=5)  # Round digits:
          enable( "dataSaveBtn" ) # allow user to then save the data
      },
          error = function(e) { stop(safeError(e)) }  # return a safeError 
      ) # end of try/catch
      
      return( df )
    }) # end renderDataTable

    # Click Next: (Format data, add data to dictionary)
    observeEvent( input$dataSaveBtn, {
      # Save data
      vals$og_data[[input$file1$name]] = vals$data_temp
      
      # Hide warnings
      hide("warning_first_column") 
      hide("warning_duplicates")
      hide("warning_numbers")
      
      # Update interface
      output$og_data = renderText({  toString( keys(vals$og_data) )  })
      hide("chooseFile")
      if( length(keys(vals$og_data))>1 ){ shinyjs::show("harmonizeBtn") }
      shinyjs::show("goHarmonize")
      disable("dataSaveBtn")
    })  
    
    # Click Harmonize
    observeEvent( input$addDataBtn, {
      reset("file1") # reset inputFile: https://gist.github.com/bborgesr/07406b30ade8a011e59971835bf6c6f7
      hide("goHarmonize")
      shinyjs::show("chooseFile")
    })  
    
    ######################################################################
    #  Harmonize
    ######################################################################
  
    #Harmonize Button
    observeEvent( input$harmonizeBtn, {
      shinyjs::hide("goHarmonize")  
      shinyjs::show("TPMorCounts")  
    })  # ends observeEvent
    
    # Counts Button
    observeEvent( input$countsBtn, {
      updateSelectInput(session, "page", selected = 'run')
      
      vals$start_raw = TRUE
      vals$tsne_raw = TRUE
      vals$start_harmonized = TRUE
      vals$tsne_harmonized = TRUE
      vals$start_main_plot = TRUE
    })  # ends observeEvent
    
    # TPM Button
    observeEvent( input$tpmBtn, {
      updateSelectInput(session, "page", selected = 'run')
      
      vals$isTPM = TRUE
      vals$start_raw = TRUE
      vals$tsne_raw = TRUE
      vals$start_harmonized = TRUE
      vals$tsne_harmonized = TRUE
      vals$start_main_plot = TRUE
    })  # ends observeEvent
    
    # Helper function: Make Boxplot
    makeBoxplot <- function( title=FALSE ){
      small_data= isolate( as.data.frame(matrix( nrow=10000,  ncol= length(unique(vals$batch)) )) )
      names(small_data) = isolate( keys(vals$og_data) )
      for( i in isolate( unique(vals$batch) ) ){ 
        small_data[i] = isolate(   sample( vals$rna_data[, which(vals$batch==i) ] ,10000,replace=TRUE)  )
      } 
      
      # Get Color Pallete
      if( isolate(length(unique(vals$batch))) == 2 ){ # set plot colors
        color_pallete =  add.alpha( brewer.pal( n=3, name="Set1")[1:2], .5)   #c("red","blue") #c("#132B43", "#56B1F7")
      }else{
        color_pallete = add.alpha( brewer.pal(n= isolate(length(unique(vals$batch))), name = "Set1" ) ,.5)
      }
      
      if( title==TRUE){
        if( isolate(vals$isTPM) ){
          ylab = "Transcripts per Million"
        }else{
          ylab = "Counts"
        }
        return( boxplot( small_data, outline=FALSE, ylab=ylab, col=color_pallete,  main="Boxplot of Values by Dataset") )
      }else{
        return( boxplot( small_data, outline=FALSE, col=color_pallete) )
      }
    }
    
    makeTSNEPlot <- function( title=FALSE ){
      
        set.seed(200)
        tsne_results <-  Rtsne( t(isolate(vals$rna_data)), perplexity=5 )
        
        # Get Color Pallete
        if( isolate(length(unique(vals$batch))) == 2 ){ # set plot colors
          color_pallete = brewer.pal( n=3, name="Set1")[1:2]   #c("red","blue") #c("#132B43", "#56B1F7")
        }else{
          color_pallete = brewer.pal(n= isolate(length(unique(vals$batch))), name = "Set1" )
        }
        
        # Set colors for points
        my_colors= c()
        for( i in isolate(vals$batch) ){
          my_colors= c( my_colors, color_pallete[i] )
        }
        
        if( title==TRUE){
          my_title= "t-SNE Plot of Samples"
          xlab=expression('tSNE'[1])
          ylab=expression('tSNE'[2])
        }else{
          my_title= "" 
          xlab=""
          ylab=""
        }
        
        return( plot(tsne_results$Y, bg=as.factor(isolate(vals$batch)), xlab=xlab, ylab=ylab, pch = 19, 
                cex=1.5, col= my_colors, main=my_title )  )  
    }
    
    
    ######################################################################
    #  Harmonize: Raw
    ######################################################################
      
    output$example_raw <- renderPlot({
      
      if( vals$start_raw == FALSE ) return()  # only run this code after vals$start is triggered
      
      # Set up Progress Bar
      vals$progress  =  Progress$new(session)
      vals$progress_step =  1/(3 + length(keys(vals$og_data)) )
      vals$progress$inc( vals$progress_step, detail="Preparing Data" )
      
      # Backend: Create batch varaible for combat_seq
      batch=c()
      count=0
      for( key_name in keys(vals$og_data) ){  #print(length( unique_names ))
        count = count + 1
        batch = c( batch, rep( count, dim(vals$og_data[[key_name]])[2]-1  )  )
      }
      vals$batch = batch

      # Backend: Get list of shared genes # Order rows of datasets 
      shared_genes = vals$og_data[[ keys(vals$og_data)[1] ]][,1] # set to names of first dataset
      for( key_name in keys(vals$og_data) ){  #print(length( unique_names ))
        row_names = vals$og_data[[key_name]][,1] # Get row names
        shared_genes = intersect( shared_genes, row_names  ) # only keep shared names across dataset
        vals$og_data[[key_name]] = vals$og_data[[key_name]][ order(row_names) ,] # Reorder row names for next for loop
        #print( all( vals$og_data[[key_name]][ order(vals$og_data[[key_name]][,1]) ,1] ==  vals$og_data[[key_name]][  ,1] ))
      }
      #isolate( vals$shared_genes = shared_genes )
      vals$shared_genes = shared_genes
      
      # Backend: Only keep common row #remove duplicates, #combine datasets
      print("made it here 2")
      rna_data = data.frame(  gene=shared_genes, stringsAsFactors=FALSE )
      for( key_name in keys(vals$og_data) ){                  #print(length( unique_names ))    #print(key_name)
        names(vals$og_data[[key_name]])[1] = "genes_list"     #change name of genes list so we can use subset function    #print(dim(vals$og_data[[key_name]]))
        vals$og_data[[key_name]] = subset( vals$og_data[[key_name]], genes_list %in% shared_genes )  #print(dim(vals$og_data[[key_name]]))
        vals$og_data[[key_name]] = vals$og_data[[key_name]][ !duplicated(vals$og_data[[key_name]][,1]) ,] #print(dim(vals$og_data[[key_name]]))
        # Commented line below sums duplicates rather than just removing them. Probably better, but takes long to run
        #vals$og_data[[key_name]] = ddply( vals$og_data[[key_name]], "genes_list", numcolwise(sum) )  # sum duplicates
        rna_data = cbind( rna_data, vals$og_data[[key_name]][,-1] )
      }
      vals$gene_names = rna_data[,1] 
      vals$subject_names = names(rna_data)[-1] 
      updateSelectizeInput(session, 'select_gene', choices= rna_data[,1], server= TRUE) # selected= rna_data[1,1]  # update gene selection list
      
      # Front & Back: Create Plot before Formatting
      #output$test_text = renderText({ "Plot Original Data" }) 
      rna_data = as.matrix(rna_data[,-1]) #make count matrix for combat seq
      vals$rna_data = rna_data 
      makeBoxplot(title=TRUE) 
    })
  
    output$tsne_raw <- renderPlot({
      if( vals$tsne_raw == FALSE ) return()
      makeTSNEPlot( title=TRUE )
      
      # Get Color Pallete
      if( isolate(length(unique(vals$batch))) == 2 ){ # set plot colors
        color_pallete = brewer.pal( n=3, name="Set1")[1:2]   #c("red","blue") #c("#132B43", "#56B1F7")
      }else{
        color_pallete = brewer.pal(n= isolate(length(unique(vals$batch))), name = "Set1" )
      }
      legend("topright", legend=isolate( keys(vals$og_data) ), col=color_pallete, pch=19 )
    })
    
    ######################################################################
    #  Harmonize: Combat Seq
    ######################################################################
    
    # Create Box Plot
    output$example_harmonized <- renderPlot({
      req( vals$start_harmonized )
      
      #vals$progress$inc( vals$progress_step, detail="Starting Combat-Seq" )
      rna_data = isolate( round( vals$rna_data ) ) # convert data back to integers
      
      # IF TPM - turn into positive integers for combat-seq
      if( isolate(vals$isTPM) ){
        scaling_factor=colSums(rna_data)/1000000 /100 #/1000   # first make all columns add up to a milion, then multiply by large constant
        rna_data = sweep( rna_data, MARGIN = 2, STATS = scaling_factor, FUN = "/")
      }
      rna_data = round(rna_data)

      mode(rna_data) <- "integer" 
      rna_data[ rna_data<0 ] = 0 # Remove Negative Values
      rna_data <- isolate( ComBat_seq( rna_data, batch=vals$batch, group=NULL ) ) # run combat-seq
      
      # IF TPM - turn results back into tpm
      if( isolate(vals$isTPM) ){
        scaling_factor=colSums(rna_data)/1000000 
        rna_data = sweep( rna_data, MARGIN = 2, STATS = scaling_factor, FUN = "/")
      }
      vals$rna_data = rna_data

      makeBoxplot() 
    })
    
    output$tsne_harmonized <- renderPlot({
      req( vals$tsne_harmonized )

      # Show messsages for screen
      shinyjs::show("combat_text_1")
      shinyjs::show("combat_text_3")
      on.exit(vals$progress$close())
  
      makeTSNEPlot()
    })  
    
    
    
    ######################################################################
    #  View Final Results
    ######################################################################
    
    observeEvent( input$viewResultsBtn , {
      updateSelectInput(session, "page", selected = 'results')
    })
    
    # Create Main Plot
    output$main_plot <- renderPlotly({
      
      req( vals$start_main_plot )
      
      # Front & Back: Create Object to plot
      gene_num= which( vals$gene_names == input$select_gene )

      # Get batch values with names instead of numbers
      batch=c()
      for( key_name in keys( isolate(vals$og_data) ) ){  #print(length( unique_names ))
        batch = c( batch, rep( key_name, dim(vals$og_data[[key_name]])[2]-1  )  )
      }
      
      if( isolate(vals$isTPM) ){
        ylab="Transcripts per Million"
      }else{
        ylab="Counts"
      }
      
      to_graph = data.frame( Expression=vals$rna_data[gene_num,], Dataset=as.factor(batch), stringsAsFactors=FALSE )# , Subject=vals$subject_names,
      par(las=1, bty='l', lwd=2, family  = 'sans', cex.axis=1.25, cex.lab=1.25, cex.main=1.75)
      p = ggplot(data= to_graph, aes(x=Dataset, y =Expression, fill=Dataset) )  +  #, subject=Subject
        labs(y = ylab, x="" ,  title = paste( input$select_gene, "Expression" ) ) +
        geom_jitter(width = 0.15, alpha = 0.75 , stroke = .3, size=2 )  + #height = 0, opacity=.5
        geom_violin( trim= FALSE, alpha = 0.5, inherit.aes = T, show.legend=FALSE ) +
        theme_gdocs()     
      
      print("IS TPM?")
      print(isolate(vals$isTPM))
      
      ggplotly(p, tooltip = c("Expression" ))
    })  
    
    # Download Data Function
    output$download_data <- downloadHandler(
       filename = function() {
         paste("harmony_rna_data.csv")
       },
       content = function(con) {
         data = cbind( vals$gene_names, vals$rna_data )
         write.csv(data, con, row.names=FALSE )
       }
    )
    
  } # ends Server
  

  #########################
  
  #########################
  #########################
  #  RUN SHINY
  #########################
  #########################
  
  # Run the application 
  options(shiny.maxRequestSize=40*1024^2) # Allow Files up to 30mb
  options( shiny.host="0.0.0.0" )
  options( shiny.port= 8888 )
  options( shiny.launch.browser = FALSE )
  
  shinyApp(ui = ui, server = server)
  
  #########################
  
  
  
