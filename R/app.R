simple_theme_grid <- ggplot2::theme_bw() +
  ggplot2::theme(
    panel.border = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "grey90"),
    panel.grid.minor = ggplot2::element_line(colour = "grey95"),
    axis.line = ggplot2::element_line(colour = "black")
  )
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

max_points <- 10
pdat <- readRDS(file="../data/prodat.med.rds")
DEdat <- readRDS(file="../data/res.rds")
DEdat$"-log10(P.Value)" <- -log10(DEdat$P.Value)

plotVolcano <- function(res, bins=80, xmax=NULL, ymax=NULL, marginal.histograms=FALSE, text.size=12, show.legend=TRUE) {
  
  tr <- attr(res, "transform.fun")
  conds <- attr(res, "conditions")
  xlab <- ifelse(is.null(tr), "FC", paste(tr, "FC"))
  tit <- paste(conds, collapse=":")
  id <- names(res)[1]
  
  g <- ggplot(res, aes(logFC, -log10(P.Value)))
  g <- g + geom_point()
  g <- g + simple_theme_grid
  g <- g + geom_vline(colour='red', xintercept=0) +
    theme(text = element_text(size=text.size)) +
    labs(x=xlab, y="-log10 P", title=tit)
  
  
  if(!is.null(xmax)) g <- g + scale_x_continuous(limits = c(-xmax, xmax), expand = c(0, 0))
  if(!is.null(ymax) ) g <- g + scale_y_continuous(limits = c(0, ymax), expand = c(0, 0))
  
  if(marginal.histograms) g <- ggExtra::ggMarginal(g, size=10, type = "histogram", xparams=list(bins=100), yparams=list(bins=50))
  return(g)
}

manage.pkg<-function(package_name){
  if(!package_name%in%installed.packages()){
    install.packages(package_name)
  }
  library(package_name,character.only = TRUE)
}

# Load required packages
lapply(c("shiny","ggplot2","dplyr","DT","gplots"), function(x) manage.pkg(x))

#######################################################################

ui <- shinyUI(fluidPage(
  
  # Application title
  titlePanel("PlotVolcano Shiny app"),
  
  fluidRow(
    column(5, plotOutput("plotVolcano", height = "700px", width = "100%", brush = "plot_brush",hover="plot_hover")),
    column(7,
           fluidRow(htmlOutput("proteinInfo")),
           fluidRow(
             column(4,
                    radioButtons("intensityScale","Intesity Scale:",choices = c("Linear scale" = "","Log scale"="Log"),inline = TRUE)
             )
           ),
           fluidRow(
             column(4,
                    fluidRow(plotOutput("jitterPlot", height = "400px",width = "100%")),
                    fluidRow(htmlOutput("gap")),
                    fluidRow(tableOutput("significanceTable"))
             ),
             column(4,
                    fluidRow(tableOutput("replicateTable"))
             ),
             column(4,
                    fluidRow(plotOutput("heatMap",height = "700px" ,width = "100%"))
             )
           )
    ),
    fluidRow(
      column(6, htmlOutput("proteinTable"))
    )
  ),
  
  # Show main protein table
  fluidRow(
    column(width = 12,
           DT::dataTableOutput("allProteinTable"))
  )
)
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #function to fetch selected proteins from Volcano plot or table
  selectProtein <- function(data,max_hover=1){
    # print('selectProtein method')
    sel = -1
    tab_idx <- as.numeric(input$allProteinTable_rows_selected)
    if(!is.null(input$plot_brush)){
      brushed <- na.omit(brushedPoints(DEdat, input$plot_brush))
      sel <-as.numeric(rownames(brushed))
    }else if(!is.null(input$plot_hover)){
      near <- nearPoints(DEdat,input$plot_hover,threshold = 20,maxpoints = max_hover)
      sel <- as.numeric(rownames(near))
    }else if(length(tab_idx)>0){
      sel <- tab_idx
      return(sel)
    }
  }
  
  #ProteinInfo
  output$proteinInfo <- renderUI({
    # print('proteinInfo method')
    sel <- selectProtein(pdat$tab)
    n <- length(sel)
    if (n == 1 && sel > 0){
      name <- paste0('<H3>', as.character(rownames(pdat$tab)[sel]),'</H3>')
      descr<-""
      HTML(paste0(name,descr,'<hr/>'))
    }else if (n > 1 && n <= max_points){
      HTML(paste0('<H3>','selection of ', n, ' proteins', '</H3><hr/>'))
    }else if (n > max_points){
      HTML(paste0('<H3>','only ',max_points,' points can be selected', '</H3><hr/>'))
    }
  })
  
  output$gap <- renderUI({HTML('<br/>')})
  
  # replicateTable
  output$replicateTable <- renderTable({
    # print('replicateTable method')
    sel <- selectProtein(pdat$tab)
    i <- pdat$tab[sel,]
    if (input$intensityScale == 'Log'){
      i <- log10(i)
    }
    if(length(sel) > 1 && sel > 0 && length(sel) <= max_points){
      data.frame(Sample=colnames(pdat$tab),Intensity=colMeans(i,na.rm = TRUE))
    }
    if (length(sel) == 1 && sel > 0){
      data.frame(Sample=colnames(pdat$tab),Intensity=i)
    }
  },digits = 1, width = "80px"
  )
  
  #significanceTable
  output$significanceTable <- renderTable({
    sel <- selectProtein(pdat$tab)
    if(length(sel) == 1 && sel > 0){
      data.frame(Contrast="1112-BMO", pvalue=as.numeric(DEdat$P.Value[sel]))
    }
  },digits = 4,width = "50px"
  )
  
  #heatMap
  output$heatMap <- renderPlot({
    # print('heatMap method')
    sel<-selectProtein(pdat$tab)
    if(length(sel) > 1 && sel > 0 && length(sel) <= max_points){
      d <- as.matrix(pdat$tab[sel,])
      mean <- rowMeans(d,na.rm = TRUE)
      mean[mean == 0] <- 1
      d <- d/mean
      d[is.nan(d)] <- NA
      row.names(d) <- rownames(pdat$tab)[sel]
      heatmap.2(d, na.rm=TRUE, dendrogram = "row",key=FALSE,keysize = 1,lhei = c(1,100),Colv = FALSE,srtRow = -35,cexRow = 1.0,na.color = "black")
    }
  })
  
  #jitterPlot
  output$jitterPlot <- renderPlot({
    # print('jitterPlot method')
    sel <- selectProtein(pdat$tab)
    if(length(sel)>0 && sel > 0 && length(sel) <= max_points){
      dataIntensity<-pdat$tab[sel,]
      if (input$intensityScale == 'Log'){
        dataIntensity<-log10(dataIntensity)
      }
      if (length(sel) == 1){
        i <- dataIntensity
      }else{
        i <- colMeans(dataIntensity,na.rm = TRUE)
      }
      s <- sapply(as.data.frame(dataIntensity),function(x) sd(x,na.rm = TRUE)/sqrt(length(x)))
      n <- length(sel)
      p <- data.frame(
        intensity= i,
        lo=i-s,
        up=i+s,
        condition=factor(pdat$metadata$condition,levels=unique(pdat$condition)),
        replicates=factor(pdat$metadata$replicate)
      )
      #shape
      p$shape <- rep(21,length(p$intensity))
      p$shape[which(p$intensity==0)] <- 24
      pd <- position_dodge(width = 0.4)
      #colorblind friendly definition
      cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
      ggplot(p, aes(x=condition, y=intensity, ymin=lo, ymax=up,colour=replicates,  shape= shape, fill=replicates)) +
        theme(text = element_text(size=20), legend.position = "bottom",legend.direction = "horizontal") +
        {if (input$intensityScale == '') ylim(0, NA)} +
        geom_point(position=pd, size=4) +
        {if(n > 1) geom_errorbar(position=pd, width = 0.1)} +
        scale_shape_identity() +  # necessary for shape mapping
        scale_fill_manual(values=cbPalette) +
        {if (input$intensityScale == 'Log') labs(x = 'Condition', y = 'Log Intensity') else labs(x = 'Condition', y = 'Intensity')}
    }
  })
  #Volcano plot
  output$plotVolcano <- renderPlot({
    # print('plotVolcano method')
    tab_idx <- as.numeric(input$allProteinTable_rows_selected)
    pVol <- plotVolcano(DEdat)
    if (length(tab_idx) > 0){
      pVol <- pVol + geom_point(data=DEdat[tab_idx,],size=3,color='red')
    }
    pVol
  })
  
  #AllProteinTable
  output$allProteinTable <-DT::renderDataTable({
    # print('allProteinTable method')
    d <- data.frame(ProteinId=DEdat$protein,mean_1112=formatC(DEdat$mean_1112),mean_BMO=formatC(DEdat$mean_BMO))
    datatable(
      d, class = 'cell-border strip hover'
    ) %>% formatStyle(0, cursor = 'pointer')
  })
}

# Run the application
shinyApp(ui = ui, server = server)

