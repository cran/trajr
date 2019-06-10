library(shiny)
library(trajr)

ui <- fluidPage(

  # App title ----
  titlePanel("Smooth & Shiny Trajectory"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: smoothing p parameter
      sliderInput(inputId = "smoothP",
                  label = "P:",
                  min = 3,
                  max = 11,
                  value = 3),

      # Input: smoothing n parameter
      sliderInput(inputId = "smoothN",
                  label = "N:",
                  min = 5,
                  max = 101,
                  value = 3)

    ),

    # Main panel for displaying outputs ----
    mainPanel(
      fluidRow(
        column (6,
          plotOutput(outputId = "trajPlot")
        ),
        column (6,
          plotOutput(outputId = "trajPlotSmooth")
        )
      ),
      fluidRow(
        column (6,
          plotOutput(outputId = "speedPlot")
        ),
        column (6,
          plotOutput(outputId = "speedPlotSmooth")
        )
      )
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {

  stoppedSpeed <- .01
  rawColour <- "#ffaa88"
  smoothColour <- "#aa88ff"

  trj <- trajr::TrajFromCoords(read.csv("eg.csv", stringsAsFactors = FALSE), timeCol = "Time")
  # Force n to be odd
  smoothed <- reactive({
    TrajSmoothSG(trj, input$smoothP, input$smoothN + (input$smoothN + 1) %% 2)
    })

  output$trajPlot <- renderPlot({
    plot(trj, lwd = 2, col = rawColour, main = "Raw Trajectory")
  })

  output$trajPlotSmooth <- renderPlot({
    plot(smoothed(), lwd = 2, col = smoothColour, main = sprintf("Smoothed, p = %d, n = %d", input$smoothP, input$smoothN + (input$smoothN + 1) %% 2))
  })

  output$speedPlot <- renderPlot({
    intervals <- TrajSpeedIntervals(trj, slowerThan = stoppedSpeed)
    title <- sprintf("Raw speed, < %g (%s/%s) highlighted", stoppedSpeed, TrajGetUnits(trj), TrajGetTimeUnits(trj))
    plot(intervals, col = rawColour, main = title)
  })

  output$speedPlotSmooth <- renderPlot({
    ylim <- range(TrajDerivatives(trj)$speed, na.rm = TRUE)
    intervals <- TrajSpeedIntervals(smoothed(), slowerThan = stoppedSpeed)
    title <- sprintf("Smoothed, < %g (%s/%s) highlighted",
                     stoppedSpeed, TrajGetUnits(trj), TrajGetTimeUnits(trj))
    plot(intervals, col = smoothColour, ylim = ylim, main = title)
  })
}


shinyApp(ui = ui, server = server)
