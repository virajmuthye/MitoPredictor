library(shiny)
library(readr)
library(tidyverse)
library(seqRFLP)
library(ggplot2)
library(Rcpp)
library(dplyr)
library(DT)
#library(ape)
#library(msaR)
library(shinythemes)

#read in the datasets
setwd("C:/Users/Viraj/Box Sync/MitochondrialProteome/mitopredictor")
MyTable <- read.csv("FINAL.MATRIX", header = TRUE, row.names = NULL)

#extract colum informations
col1 <- MyTable$species
col2 <- MyTable$tp
col3 <- MyTable$rc
col4 <- MyTable$mf
col5 <- MyTable$geneid
col6 <- MyTable$mito

#part 1 : the UI page
ui <- fluidPage(theme = shinytheme("cerulean"),
                
                titlePanel("Mitopredictor"),
                
                #Define the tab panels
                
                mainPanel(
                  tabsetPanel(type = "tabs",
                              tabPanel("About",
                                       br(),
                                       strong("About"),
                                       p("This is a ShinyR app for comparative analysis of mitochondrial proteomes of human, mouse, C.elegans and D.melanogaster, yeast
                                         and the query species used.")
                                       ), 
                      
                              tabPanel("Mitoproteomes",
                                       sidebarPanel(
                                         checkboxGroupInput("species",label="Select species:",choices=sort(unique(col1)), selected = sort(unique(col1))),
                                         checkboxGroupInput("mito",label="Mitochondrial / Non-mitochondrial proteins:",choices = sort(unique(col6)), selected = sort(unique(col6))),
                                         checkboxGroupInput("tp",label="TargetP prediction",choices = sort(unique(col2)), selected = sort(unique(col2))),
                                         checkboxGroupInput("rc",label="TargetP RC",choices = sort(unique(col3)), selected = sort(unique(col3))),
                                         checkboxGroupInput("mf",label="MitoFates prediction",choices = sort(unique(col4)), selected = sort(unique(col4))),
                                         br(),
                                         downloadButton(outputId = "download_data_4", label = "Download table")
                                       ),
                                       sidebarPanel(
                                         selectInput(inputId = "geneid",label="Select OG:",choices=sort(unique(col5)),multiple = FALSE),
                                         plotOutput(outputId = "plot"),
                                         width = 8
                                       ),
                                       DT::dataTableOutput(outputId = "view")
                              )
                              )
                  )
)

# part 2 Define server logic
server <- function(input, output) {

################################# presequence tabset #####################################
  output$view <- DT::renderDataTable({
    
    req(input$species)
    
    table_select <- MyTable %>% filter(species %in% input$species)  %>% filter(geneid %in% input$geneid) %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf)
    
    DT::datatable(data = table_select,
                  options = list(pageLength =500 ), 
                  rownames = FALSE)
    
  })
  
  df_subset_4 <- reactive({
    MyTable[col1 %in% input$species,] %>% filter(geneid %in% input$geneid) %>% filter(species %in% input$species) %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(mf %in% input$mf) 
      })
  
  df_subset_5 <- reactive({
    MyTable %>% filter(species %in% input$species) %>% filter(geneid %in% input$geneid) %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf)
  })
  df_subset_6 <- reactive({
    MyTable %>% filter(species %in% input$species) %>% filter(geneid %in% input$geneid)  %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf)
  })
  
  df_subset_7 <- reactive({
    MyTable %>% filter(species %in% input$species) %>% filter(geneid %in% input$geneid)  %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf)
  })
  
  output$plot = renderPlot({
    p <- ggplot(df_subset_4(),aes(factor(species)))
    p + geom_bar(fill = "#e6add8") + stat_count(aes(label=..count..), vjust=5, geom="text", position="identity", size=10) + theme_minimal()+theme(axis.text.x=element_text(size=rel(2)),axis.text.y=element_text(size=rel(2)), axis.title.x=element_text(size=rel(1.5)),axis.title.y=element_text(size=rel(1))) +
      xlab("Species") 
  })
  
  output$download_data_4 <- downloadHandler(
    filename = function() {
      paste("download.csv")
    },
    content = function(file) {
      write.csv(df_subset_7(), file, quote = FALSE)
    }
  ) 
}

# part 3 Bind ui and server together
shinyApp(ui=ui, server=server)