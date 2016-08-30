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
        
        output$hot <- renderRHandsontable({
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
        
      }
    }
  )
  
  output$nameVar <- renderRHandsontable({
    DF <- data.frame(matrix(ncol=1,nrow=as.numeric(input$numVar)))
    DF[sapply(DF, is.logical),] = lapply(DF[sapply(DF, is.logical)], as.character)
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
      rHandsontableOutput("nameVar"),
      helpText("Enter the names of the condition variables")
    ),
    mainPanel(
      fluidRow(
        column(
          width = 10,
          offset = 1,
          rHandsontableOutput("hot"),
          selectInput("levelAnalysis", label = h3("Select level of analysis"), 
                      choices = list("transcript" = 1, "gene" = 2), selected = 1),
          selectInput("typeTest", label = h3("Select test"), 
                      choices = list("likelihood" = 1, "wald" = 2), selected = 1),
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