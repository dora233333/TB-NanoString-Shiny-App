library(shiny)
library(ggplot2)
library(uwot)
library(NMF)
library(mclust)

# Load the TB NanoString data
TBnanostring <- readRDS(
  "TBnanostring.rds"
)

# Separate TB status from the gene-expression variables
tb_status <- TBnanostring$TB_Status

gene_matrix <- as.matrix(
  TBnanostring[, -1]
)

# ------------------------------------------------------------
# Compute dimension reductions once when the app starts
# ------------------------------------------------------------

# PCA
pca_model <- prcomp(
  gene_matrix,
  center = TRUE,
  scale. = FALSE
)

pca_coordinates <- data.frame(
  Dim1 = pca_model$x[, 1],
  Dim2 = pca_model$x[, 2]
)

# UMAP
set.seed(2026)

umap_result <- umap(
  gene_matrix,
  n_neighbors = 15,
  n_components = 2,
  metric = "euclidean"
)

umap_coordinates <- data.frame(
  Dim1 = umap_result[, 1],
  Dim2 = umap_result[, 2]
)

# Shift each gene so that NMF receives non-negative data
gene_minimums <- apply(
  gene_matrix,
  2,
  min
)

nmf_matrix <- sweep(
  gene_matrix,
  MARGIN = 2,
  STATS = gene_minimums,
  FUN = "-"
)

# NMF
set.seed(2026)

nmf_model <- nmf(
  nmf_matrix,
  rank = 2,
  nrun = 20,
  seed = 2026
)

nmf_basis <- basis(
  nmf_model
)

nmf_coordinates <- data.frame(
  Dim1 = nmf_basis[, 1],
  Dim2 = nmf_basis[, 2]
)

# Store the three dimension-reduction results
dimension_reductions <- list(
  PCA = pca_coordinates,
  UMAP = umap_coordinates,
  NMF = nmf_coordinates
)

# Define colors for as many as nine clusters
cluster_colors <- c(
  "#E41A1C",
  "#377EB8",
  "#4DAF4A",
  "#984EA3",
  "#FF7F00",
  "#A65628",
  "#F781BF",
  "#999999",
  "#66C2A5"
)

# ------------------------------------------------------------
# User interface
# ------------------------------------------------------------

ui <- fluidPage(
  titlePanel(
    "TB NanoString Dimension Reduction and K-means Clustering"
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      selectInput(
        inputId = "method",
        label = "Dimension reduction method",
        choices = c(
          "PCA",
          "UMAP",
          "NMF"
        ),
        selected = "PCA"
      ),
      
      numericInput(
        inputId = "clusters",
        label = "Number of k-means clusters",
        value = 2,
        min = 1,
        max = 9,
        step = 1
      ),
      
      helpText(
        "Points are colored by k-means cluster and shaped by TB status. Black crosses show cluster centers."
      )
    ),
    
    mainPanel(
      width = 9,
      
      fluidRow(
        column(
          width = 8,
          
          plotOutput(
            outputId = "cluster_plot",
            height = "550px"
          )
        ),
        
        column(
          width = 4,
          
          h3(
            "Cluster results"
          ),
          
          h4(
            "Clusters versus TB status"
          ),
          
          tableOutput(
            outputId = "cluster_table"
          ),
          
          h4(
            "Adjusted Rand Index"
          ),
          
          textOutput(
            outputId = "ari_text"
          )
        )
      )
    )
  )
)

# ------------------------------------------------------------
# Server
# ------------------------------------------------------------

server <- function(input, output, session) {
  
  # Select the requested two-dimensional representation
  selected_data <- reactive({
    dimension_reductions[[input$method]]
  })
  
  # Run k-means on the selected two dimensions
  clusters <- reactive({
    set.seed(2026)
    
    kmeans(
      selected_data(),
      centers = input$clusters,
      nstart = 25
    )
  })
  
  # Combine coordinates, cluster assignments, and TB status
  plot_data <- reactive({
    data.frame(
      selected_data(),
      Cluster = factor(
        clusters()$cluster
      ),
      TB_Status = tb_status
    )
  })
  
  # Create the dimension-reduction plot
  output$cluster_plot <- renderPlot({
    centers_data <- as.data.frame(
      clusters()$centers
    )
    
    colnames(centers_data) <- c(
      "Dim1",
      "Dim2"
    )
    
    ggplot(
      plot_data(),
      aes(
        x = Dim1,
        y = Dim2,
        color = Cluster,
        shape = TB_Status
      )
    ) +
      geom_point(
        size = 3,
        alpha = 0.8
      ) +
      geom_point(
        data = centers_data,
        aes(
          x = Dim1,
          y = Dim2
        ),
        inherit.aes = FALSE,
        color = "black",
        shape = 4,
        size = 7,
        stroke = 2
      ) +
      scale_color_manual(
        values = cluster_colors
      ) +
      scale_shape_manual(
        values = c(
          "TB" = 16,
          "LTBI" = 17
        )
      ) +
      labs(
        title = paste(
          input$method,
          "with K-means Clustering"
        ),
        subtitle = paste(
          "Number of clusters:",
          input$clusters
        ),
        x = paste0(
          input$method,
          " dimension 1"
        ),
        y = paste0(
          input$method,
          " dimension 2"
        ),
        color = "K-means cluster",
        shape = "TB status"
      ) +
      theme_minimal(
        base_size = 14
      ) +
      theme(
        legend.position = "bottom"
      )
  })
  
  # Cross-tabulate cluster assignments and TB status
  output$cluster_table <- renderTable(
    {
      cluster_status_table <- table(
        Cluster = clusters()$cluster,
        TB_Status = tb_status
      )
      
      as.data.frame.matrix(
        cluster_status_table
      )
    },
    rownames = TRUE
  )
  
  # Calculate the adjusted Rand index
  output$ari_text <- renderText({
    current_ari <- adjustedRandIndex(
      clusters()$cluster,
      tb_status
    )
    
    round(
      current_ari,
      3
    )
  })
}

# Launch the application
shinyApp(
  ui = ui,
  server = server
)