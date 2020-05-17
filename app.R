library(shiny)
library(tidyverse)
library(e1071)
library(survival)
library(ggfortify)
library(survminer)
library(shinythemes)
# Helper functions

load("eplet_risk_model.RData")

model_names = c("GenderCode",
                "agetxn",
                "donorage",
                "X276L",
                "X66IT",
                "X76ED",
                "X76L",
                "X62GK2",
                "X267PE",
                "X185I",
                "X142M3",
                "X82LR",
                "X164FR",
                "X70RT",
                "X185T",
                "X116I",
                "X17RS",
                "X157I",
                "X173K",
                "X130A",
                "X65GK")


# Define UI for application that draws a histogram
ui <- navbarPage(theme = shinytheme("yeti"),
                 
                 # Application title
                 title = "Kidney Transplant Risk Calculator",
                 
                 # I have no idea what this code means but it removes a ghost tab
                 # https://github.com/rstudio/shiny/issues/827
                 
                 id = "page-nav",
                 tags$head(tags$style(HTML(
                   "#page-nav > li:first-child { display: none; }"
                 ))),
                 
                 # Sidebar with a slider input for number of bins 
                 sidebarLayout(
                   sidebarPanel(
                     textInput("patient_name", 
                               "Patient name:"),
                     
                     numericInput("patient_age",
                                  "Age of patient:",
                                  value = 1,
                                  min = 0,
                                  max = 120),
                     
                     radioButtons("patient_gender",
                                  "Gender of patient",
                                  c("Female" = 0, 
                                    "Male" = 1),
                                  inline = TRUE),
                     
                     numericInput("donor_age",
                                  "Age of donor:",
                                  value = 1,
                                  min = 0,
                                  max = 120),
                     
                     checkboxGroupInput("eplets",
                                        "Mismatched eplets:",
                                        c("X276L" = "X276L",
                                          "X66IT" = "X66IT",
                                          "X76ED" = "X76ED",
                                          "X76L" = "X76L",
                                          "X62GK2" = "X62GK2",
                                          "X267PE" = "X267PE",
                                          "X185I" = "X185I",
                                          "X142M3" = "X142M3",
                                          "X82LR" = "X82LR",
                                          "X164FR" = "X164FR",
                                          "X70RT" = "X70RT",
                                          "X185T" = "X185T",
                                          "X116I" = "X116I",
                                          "X17RS" = "X17RS",
                                          "X157I" = "X157I",
                                          "X173K" = "X173K",
                                          "X130A" = "X130A",
                                          "X65GK" = "X65GK")),
                     
                     actionButton("submit_button",
                                  "Submit",
                                  class = "btn-primary")
                     ),
                   
                   # Show a plot of the generated distribution
                   mainPanel(
                     tabsetPanel(type = "tabs",
                                 tabPanel("Risk Calculator",
                                          plotOutput("pred_quantile_plot")),
                                 
                                 tabPanel("Survival Analysis",
                                          plotOutput("survPlot"))),
                     
                     downloadButton("report", 
                                    "Generate report",
                                    class = "btn-primary")
                     )
                   )
                 )

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  patient_data = eventReactive(input$submit_button, {
    
    eplet_binary = ifelse(model_names[4:length(model_names)] %in% input$eplets, 1, 0) %>% 
      data.frame() %>% 
      t() 
      
      patient_dat = data.frame(x1 = as.numeric(input$patient_gender),
                               x2 = as.numeric(input$patient_age),
                               x3 = as.numeric(input$donor_age))
      
      patient_dat = cbind(patient_dat, eplet_binary)
      
      colnames(patient_dat) = model_names
      
      patient_dat
    })
    
    output$pred_quantile_plot = renderPlot({
      
      prediction  = predict(svm_mod, patient_data(), probability=TRUE)
      svm_probs = attr(prediction, "probabilities")
      
      rej_prob = svm_probs[1,2]
      
      ind = min(which(qq_plot_df$rej_prob > rej_prob))
      
      quantile_patient = qq_plot_df$quantile[ind]
      
      qq_plot_df$outcome = ifelse(qq_plot_df$outcome == 0, "Stable", "Reject")
      
      qq_plot_df %>% ggplot(aes(x = rej_prob, y = quantile, colour = outcome)) +
        geom_point() +
        scale_colour_manual(values = c("red", "green")) +
        geom_point(data = data.frame(rej_prob,
                                     quantile_patient), 
                   mapping = aes(x = rej_prob, y = quantile_patient),
                   colour = "black",
                   size = 8,
                   shape = 4,
                   stroke = 2) + 
        labs(x = "Rejection Probability", 
             y = "Rejection Quantile", 
             colour = "Rejection Outcome") +
        ylim(0, 1) + 
        xlim(0, 1) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0)
    })
    
    output$survPlot = renderPlot({
      survfit(km_fit, patient_data(), conf.int=.9) %>% 
        autoplot() + 
        labs(x = "Time (years)",
             y = "Survival Probability") +
        ylim(0, 1) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0)
    })
    
    output$report <- downloadHandler(
      filename = "report.pdf",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        
        tempReport <- file.path("report.Rmd")
        file.copy("report.pdf", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(name = input$patient_name, 
                       age = input$patient_age, 
                       gender = input$patient_gender, 
                       donor_age = input$donor_age, 
                       eplets = input$eplets)

        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, 
                          output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
