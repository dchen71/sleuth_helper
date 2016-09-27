#Kallisto processor to help create the sleuth objects for further downstream processing or exploratory analysis 
#using the sleuth app

library(shiny)
source("scripts/directoryInput.R")
library(rhandsontable)
library(sleuth)
library(shinyjs)


##Setup biomart
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = 'ensembl.org')

#Get back gene name data for mouse
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# Define server logic required to create the sleuth object
server = (function(input, output, session) {
  # Init hiding of loading element
  hide("loading-1")
  hide("loading-2")
  hide("loading-3")
  hide("loading-4")
  hide("loading-5")
  
  #Observer for directory input
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$directory
    },
    handlerExpr = {
      if ((input$directory) > 0) {
        # condition prevents handler execution on initial app launch
        
        # launch the directory selection dialog with initial path read from the widget
        path <- choose.dir(default = readDirectoryInput(session, 'directory'))
        
        # update the widget value
        updateDirectoryInput(session, 'directory', value = path)
        
        if(is.na(path)){
          #If cancel, then don't run the rest
          return(FALSE)
        }
        
        output$inputVariables <- renderRHandsontable({
          folders = dir(readDirectoryInput(session, 'directory'))
          DF <- data.frame(sample=folders, matrix(ncol = as.numeric(input$numVar)))
          DF[sapply(DF, is.logical)] = lapply(DF[sapply(DF, is.logical)], as.character)
          
          input_data = hot_to_r(input$nameVar)
          if(length(names(DF)) == length(c("sample", input_data$'Variable Names'))){
            names(DF) = c("sample", input_data$'Variable Names')
          }
          
          if (!is.null(DF))
            rhandsontable(DF, stretchH = "all") %>%
              hot_col("sample", readOnly = TRUE) %>%
              #Validate NA
              hot_cols(validator = "
                function (value, callback) {
                setTimeout(function(){
                  callback(
                    value != 'NA'
                  );
                }, 300)
                }", allowInvalid = FALSE)
        })
        
        #Header for input df
        output$inputHeader = renderText({
          return("Input conditions")
        })
        
        #Helper for input df
        output$inputHelper = renderText({
          return("Input the conditions for each of the experimental variables. 
                 Note that there must be a minimum of 2 samples per condition for the analysis.")
        })
        
      }
    }
  )
  
  #Table to add in the variable names
  output$nameVar <- renderRHandsontable({
    DF <- data.frame(matrix(ncol=1,nrow=as.numeric(input$numVar)))
    DF[sapply(DF, is.logical),] = lapply(DF[sapply(DF, is.logical)], as.character)
    names(DF) = "Variable Names"
    rhandsontable(DF) %>%
      #Validate NA
      hot_cols(validator = "
              function (value, callback) {
              setTimeout(function(){
                callback(value != 'NA' && value != 'path' && value != 'sample');
              }, 300)
              }", allowInvalid = FALSE)
  })
  
  #Observer to begin processing kallisto objects
  observeEvent(input$startProcess, {
    show(id="loading-1")
    
    folders = dir(readDirectoryInput(session, 'directory'))
    kal_dirs <- sapply(folders, function(id) file.path(readDirectoryInput(session, 'directory'), id, "output"))
    
    s2c = hot_to_r(input$inputVariables)
    
    #Append location of files per sample
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    
    #Find the valid column names
    variable_names = hot_to_r(input$nameVar)
    variable_names = variable_names$'Variable Names'
    
    #Check to make sure there are at least 2 factors in each col
    num_factors = TRUE
    for(i in var_names){
      if(length(levels(test[,i])) < 2){
        num_factors = FALSE
        break
      }
    }
    
    #Validate level of factors to be at least 2
    if(num_factors == FALSE){
      hide(id="loading-1")
      output$createModel = renderText({return("Error: One or more variables has less than two factors")})
      return(FALSE)
    }

    if(length(grep("abundance.h5", dir(s2c$path))) == nrow(s2c)){
      #Transcript or gene level
      if(input$levelAnalysis == "trans"){
        so <- sleuth_prep(s2c, as.formula((paste("~",paste(variable_names,collapse="+")))) , target_mapping = t2g)
      } else if(input$levelAnalysis == "gene"){
        so <- sleuth_prep(s2c, as.formula((paste("~",paste(variable_names,collapse="+")))) , target_mapping = t2g,
                          aggregation_column = 'ens_gene')
      }
      
      #Likelihood test
      so <- sleuth_fit(so)
      so <- sleuth_fit(so, ~1, 'reduced')
      so <<- sleuth_lrt(so, 'reduced', 'full')
      
      output$createModel = renderText({return("Model created")})
    } else {
      output$createModel = renderText({return("Error: One or more directories does not contain Kallisto reads")})
    }
    
    hide(id="loading-1")
    
    
  })
  
  #Save event
  observeEvent(input$saveSleuth, {
    show(id="loading-2")
    save(so, file="sleuth_object.RData")
    hide(id="loading-2")
    output$completeSave = renderText({return("Object saved")})
  })
  
  #Get kallisto abundance and save as csv
  observeEvent(input$createAbun, {
    show(id="loading-3")
    write.csv(kallisto_table(so), file="kallisto_table.csv", row.names = FALSE)
    hide(id="loading-3")
    output$completeAbun = renderText({return("Abundance table created")})
  })
  
  #Extract gene table results and save
  observeEvent(input$createTable, {
    #sleuth_results
    show(id="loading-4")
    if(input$typeTest == "lrt"){
      write.csv(sleuth_gene_table(so,"reduced:full", test_type="lrt"), file="sleuth_gene_table.csv", row.names=FALSE)
    } 
    hide(id="loading-4")
    output$completeTable = renderText({return("Sleuth gene table created")})
  })
  
  #Extract wald test results and save
  observeEvent(input$createWald, {
    #sleuth_results
    show(id="loading-5")
    if(input$typeTest == "lrt"){
      write.csv(sleuth_results(so,"reduced:full", test_type="lrt"), file="sleuth_results.csv", row.names=FALSE)
    } 
    hide(id="loading-5")
    output$completeWald = renderText({return("Sleuth results created")})
  })
})

# Define UI for application
ui = (fluidPage(
  useShinyjs(),
  titlePanel("Sleuth Helper"),

  sidebarLayout(
    sidebarPanel(
      directoryInput('directory', label = 'Select directory'),
      helpText("Select the directory that contains the quantified reads from Kallisto"),
      selectInput("numVar", label = h3("Select number of variables"), 
                  choices = list("1" = 1, "2" = 2,
                                 "3" = 3, "4" = 4,
                                 "5" = 5), selected = 1),
      helpText("Select the number of condition variables to use"),
      rHandsontableOutput("nameVar"),
      helpText("Enter unique names for the condition variables. Note that the names: path and sample, are reserved names for Sleuth.")
    ),
    mainPanel(
      fluidRow(
        column(
          width = 10,
          offset = 1,
          h3(textOutput("inputHeader")),
          rHandsontableOutput("inputVariables"),
          helpText(textOutput("inputHelper")),
          #Error with not showing until actual directory has been shown
          conditionalPanel(condition = "(output.inputVariables)",
                           selectInput("levelAnalysis", label = h3("Select level of analysis"), 
                                       choices = list("Transcript" = "trans", "Gene" = "gene"), selected = "trans"),
                           selectInput("typeTest", label = h3("Select test"), 
                                       choices = list("Likelihood Ratio Test" = "lrt"), selected = "lrt"),
                           conditionalPanel(condition= "input.typeTest == 'lrt'",
                                            helpText("The likelihood ratio test is a statistical test used to compare the goodness of fit of 
                                    two models, one of which (the null model) is a special case of the other 
                                    (the alternative model). The test is based on the likelihood ratio, 
                                    which expresses how many times more likely the data are under one model 
                                    than the other. This likelihood ratio, or equivalently its logarithm, can then 
                                    be used to compute a p-value, or compared to a critical value to decide whether 
                                    to reject the null model in favour of the alternative model. When the logarithm of 
                                    the likelihood ratio is used, the statistic is known as a log-likelihood ratio 
                                    statistic, and the probability distribution of this test statistic, assuming that 
                                    the null model is true, can be approximated using Wilks’ theorem.")),
                           actionButton("startProcess", "Create Sleuth Object"),
                           br(),
                           helpText("Create model based on parameters for further examination via Sleuth"),
                           tags$img(src="spinner.gif", id="loading-1"),
                           textOutput("createModel"),
                           br()
                           ),
          conditionalPanel(condition = "(output.createModel)",
                           actionButton("saveSleuth", "Save Sleuth Object"),
                           br(),
                           helpText("Save the object for future usage in current working directory"),
                           tags$img(src="spinner.gif", id="loading-2"),
                           textOutput("completeSave"),
                           br(),
                           actionButton("createAbun", "Create Kallisto abundance table"),
                           br(),
                           helpText("Create an abundance table containing information about transcripts, length, abundance, etc and save to current directory"),
                           tags$img(src="spinner.gif", id="loading-3"),
                           textOutput("completeAbun"),
                           br(),
                           actionButton("createTable", "Create test results"),
                           br(),
                           helpText("Creates table showing which genes most significantly mapping to transcript and save in current working firectory"),
                           tags$img(src="spinner.gif", id="loading-4"),
                           textOutput("completeTable"),
                           br(),
                           actionButton("createWald", "Create test results"),
                           br(),
                           helpText("Creates table showing test results from sleuth object and save to current directory and save to current working directory"),
                           tags$img(src="spinner.gif", id="loading-5"),
                           textOutput("completeWald"),
                           br()
                           )
        )
      )
    )
  )
  

))

#Load single file shiny app
shinyApp(ui = ui, server = server)