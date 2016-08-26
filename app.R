#Kallisto processor to help create the sleuth objects for further downstream processing or exploratory analysis 
#using the sleuth app

library(shiny)
source("directoryInput.R")

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
      }
    }
  )
  

  
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
      directoryInput('directory', label = 'Select directory', value = '~'),
      helpText("Select the directory that contains the quantified reads from Kallisto")
    ),
    mainPanel(
      fluidRow(
        column(1),
        column(
          width = 10,
          
          # Application title
          verbatimTextOutput("folders"),
          tableOutput("samples")
          
          #Want to show up table to create dataframe for conditions and a submit so it can save whatever it is
          #Want to choose transcript/gene level afterwards
          #Want to choose likelihood/wald or both and submit to save it
          #Probably options to get back object, or the whatever files it has
        ),
        column(1)
      )
    )
  )
  

))

#Load single file shiny app
shinyApp(ui = ui, server = server)