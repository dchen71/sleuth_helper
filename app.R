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
      }
    }
  )
  

})

# Define UI for application
ui = (fluidPage(
  fluidRow(
    column(1),
    column(
      width = 10,
      
      # Application title
      titlePanel("Directory Input Demo"),
      directoryInput('directory', label = 'selected directory', value = '~')

    ),
    column(1)
  )
))

#Load single file shiny app
shinyApp(ui = ui, server = server)