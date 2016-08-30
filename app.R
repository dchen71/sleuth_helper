#Kallisto processor to help create the sleuth objects for further downstream processing or exploratory analysis 
#using the sleuth app

library(shiny)
source("directoryInput.R")
library(rhandsontable)

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
        
        output$folders = renderText({
          folders = list.dirs(readDirectoryInput(session, 'directory'), full.names = T, recursive = FALSE)
          folders
        })
        
        output$samples = renderTable({
          folders = dir(readDirectoryInput(session, 'directory'))
          return(data.frame(folders))
        })
        
        output$hot <- renderRHandsontable({
          folders = dir(readDirectoryInput(session, 'directory'))
          DF <- data.frame(GENENAME=folders, matrix(ncol = as.numeric(input$numVar)))
          if (!is.null(DF))
            rhandsontable(DF, stretchH = "all") %>%
            hot_col("GENENAME", readOnly = TRUE)
        })
        
      }
    }
  )
  
  output$nameVar <- renderRHandsontable({
    DF <- data.frame(matrix(ncol=1,nrow=as.numeric(input$numVar)))
    DF[sapply(DF, is.logical),] = lapply(DF[sapply(DF, is.logical)], as.factor)
    names(DF) = "Variable Names"
    rhandsontable(DF)
  })
  
  #use output of folders to find all the actual directories with kallisto results
  #use dir() to find the name of all the sample names
  #create conditional list(setup thing to create number of conditions to setup)
  #create table to input manually the conditions
  
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
      rHandsontableOutput("nameVar")
    ),
    mainPanel(
      fluidRow(
        column(
          width = 10,
          offset = 1,
          verbatimTextOutput("folders"),
          tableOutput("samples"),
          rHandsontableOutput("hot")
          
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