#Kallisto processor to help create the sleuth objects for further downstream processing or exploratory analysis 
#using the sleuth app

library(shiny)
source("directoryInput.R")
library(rhandsontable)

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
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$directory
    },
    handlerExpr = {
      if ((input$directory) > 0) {
        # condition prevents handler execution on initial app launch
        
        # launch the directory selection dialog with initial path read from the widget
        path = choose.dir(default = readDirectoryInput(session, 'directory'))
        
        # update the widget value
        updateDirectoryInput(session, 'directory', value = path)
        
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
            hot_col("sample", readOnly = TRUE)
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
  
  output$nameVar <- renderRHandsontable({
    DF <- data.frame(matrix(ncol=1,nrow=as.numeric(input$numVar)))
    DF[sapply(DF, is.logical),] = lapply(DF[sapply(DF, is.logical)], as.character)
    names(DF) = "Variable Names"
    rhandsontable(DF)
  })
  
  test = function(){
    kal_dirs <- sapply(folders, function(id) file.path(readDirectoryInput(session, 'directory'), id, "output"))
    
    #Append location of files per sample
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    
    #Transcript likelihood
    so <- sleuth_prep(s2c, ~ condition , target_mapping = t2g)
    so <- sleuth_fit(so)
    so <- sleuth_fit(so, ~1, 'reduced')
    so <- sleuth_lrt(so, 'reduced', 'full')
    
    #Gene level
    soGene <- sleuth_fit(soGene)
    soGene <- sleuth_fit(soGene, ~1, 'reduced')
    soGene <- sleuth_lrt(soGene, 'reduced', 'full')
  }
  
  #setup transcript analysis
  #setup gene level analysis
  
  #choice of likelihood test or wald test
  #return back sleuth analyzed data
  
})

# Define UI for application
ui = (fluidPage(
  titlePanel("Sleuth Helper"),

  sidebarLayout(
    sidebarPanel(
      directoryInput('directory', label = 'Select directory'),
      helpText("Select the directory that contains the quantified reads from Kallisto"),
      selectInput("numVar", label = h3("Select number of variables"), 
                  choices = list("1" = 1, "2" = 2,
                                 "3" = 3, "4" = 4,
                                 "5" = 5), selected = 1),
      helpText("Select the number of condition variables"),
      rHandsontableOutput("nameVar"),
      helpText("Enter the names of the condition variables")
    ),
    mainPanel(
      fluidRow(
        column(
          width = 10,
          offset = 1,
          h3(textOutput("inputHeader")),
          rHandsontableOutput("inputVariables"),
          helpText(textOutput("inputHelper")),
          conditionalPanel(condition = "(input.directory) > 0",
                           selectInput("levelAnalysis", label = h3("Select level of analysis"), 
                                       choices = list("Transcript" = "trans", "Gene" = "gene"), selected = "trans"),
                           selectInput("typeTest", label = h3("Select test"), 
                                       choices = list("Likelihood Ratio Test" = "lrt", "Wald Test" = "wald"), selected = "lrt"),
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
                           conditionalPanel(condition= "input.typeTest == 'wald'",
                                            helpText("The Wald test is a parametric statistical test named after the Hungarian statistician 
                                    Abraham Wald. Whenever a relationship within or between data items can be expressed as 
                                    a statistical model with parameters to be estimated from a sample, the Wald test can be 
                                    used to test the true value of the parameter based on the sample estimate."))
                           ),
          actionButton("startProcess", "Process"),
          actionButton("goButton", "Go!"),
          actionButton("goButton", "Go!")
          
          #Want to show up table to create dataframe for conditions and a submit so it can save whatever it is
          #Want to choose transcript/gene level afterwards
          #Want to choose likelihood/wald or both and submit to save it
          #Probably options to get back object, or the whatever files it has
        )
      )
    )
  )
  

))

#Load single file shiny app
shinyApp(ui = ui, server = server)