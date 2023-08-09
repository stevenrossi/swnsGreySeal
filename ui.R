pageWithSidebar(

  headerPanel('Bayesian Population Model for Southwest Nova Scotia Grey Seal'),
  sidebarPanel(
    
    br(),"Fixed parameters",br(),
    numericInput('I0', 'I[0]', 10, min = 1, max = 10000),

    hr(),"Priors",br(),
    numericInput('rMu', 'rMu', 0.2, min= 0, max = 10),
    numericInput('rSD', 'rSD', 0.5, min = 0, max = 1e6),
    numericInput('KMu', 'kMu', 5, min = 1, max = 9),
    numericInput('KSD', 'KSD', 10, min = 1, max = 9),

    hr(),"MCMC controls",br(),
    numericInput('nIter', '# of iterations', 5e3, min = 100, max = 1e20),
    numericInput('nChain', '# of chains', 3, min = 3, max = 100),
    numericInput('adapt_delta', 'adapt_delta', 0.95, min = 0.5, max = 0.99999999),
    numericInput('max_treedepth', 'max_treedepth', 10, min = 5, max = 15),

    hr(),"Projection controls",br(),
    numericInput('finalYear', 'Final year', 2023, min = 2021, max = 3000),
    numericInput('pupAdultRatio', 'Pup:adult ratio', 0.26, min = 0, max = 1),
    width=2,  

  ),
  mainPanel(

      tabsetPanel(type = "tabs",
                  selected="Time-series",
                  tabPanel("Model Info", includeMarkdown("modelInfo.md")),
                  tabPanel("Time-series", plotOutput('timeSeries',width="100%")),
                  tabPanel("Prior/Posterior", plotOutput('priorPost',width="100%")),
                  tabPanel("Diagnostics", tableOutput('diagnos')),
      )

  )
)


