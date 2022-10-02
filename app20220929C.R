#--------------------------------------------------------------------#
#
## Dose response curves via drc ####
# R code, in a shiny wrapper
#
#
# Francis Crick Institute, London, UK
#
# Edward Carr
## 29 September 2022
#
#--------------------------------------------------------------------#
library(openxlsx)
library(tidyverse)
library(drc)
library(broom)
library(khroma)
library(svglite)
#--------------------------------------------------------------------#
server <- function(input, output) {
    
  datamAb <- reactive({
    
    # input$file.mAb will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file.mAb)
    
    # when reading semicolon separated files,
    # having a comma separator causes `read.csv` to error
    tryCatch(
      {
        df <- read.csv(input$file.mAb$datapath,
                       header = input$header,
                       sep = input$sep,
                       quote = input$quote)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
  })
  
  
    output$contentsmAb <- renderTable({
      
      if(input$disp == "head") {
        return(head(datamAb()))
      }
      else {
        return(datamAb())
      }
      
    })
    
    ####
    
    datavoc <- reactive({
      
      # input$file.voc will be NULL initially. After the user selects
      # and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.
      
      req(input$file.voc)
      
      # when reading semicolon separated files,
      # having a comma separator causes `read.csv` to error
      tryCatch(
        {
          df <- read.csv(input$file.voc$datapath,
                         header = input$header,
                         sep = input$sep,
                         quote = input$quote)
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )
      
    })
    
    
    
    output$contentsvoc <- renderTable({
      
        if(input$disp == "head") {
          return(head(datavoc()))
        }
        else {
          return(datavoc())
        }
        
      })
    #####
    
    dataSPC <- reactive({
      
      # input$file.mAb will be NULL initially. After the user selects
      # and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.
      
      req(input$file.spc)
      
      # when reading semicolon separated files,
      # having a comma separator causes `read.csv` to error
      tryCatch(
        {
          df.spc <- read.csv(input$file.spc$datapath,
                             header = input$header,
                             sep = input$sep,
                             quote = input$quote)
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )
      
      
    })
    
    
    output$contentsspc <- renderTable({
      
      
      if(input$disp == "head") {
        return(head(dataSPC()))
      }
      else {
        return(dataSPC())
      }
      
    })

    ######
    

    dataXLSX <- reactive({
      
      
      # input$file.xlsx will be NULL initially. After the user selects
      # and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.
      
      req(input$file.xlsxs)
      
      xl <- input$file.xlsxs$datapath
      # read each xlsx
      # and combine
      dat <- lapply(xl, function(x) {
        xx <- read.xlsx(xlsxFile = x)
        xx <- as_tibble(xx)
        if(!c("well.fail") %in% colnames(xx)) { xx <- xx %>% mutate(well.fail = " ")}
        return(xx)
      })
      
      dat <- do.call(rbind, dat)
      
    })
    output$contentsxlsx <- renderTable({
      
      dat <- dataXLSX()
        return(dat)
      
    })
    
  output$distPlot <- renderPlot({
    cowplot::plot_grid(
    ggplot(data=datamAb(), aes(x=mAb_plotting_order, y=mAb_name_in_raw_data))+geom_point(),
    ggplot(data=datavoc(), aes(x=VOC_plotting_order, y=VOC_name_in_raw_data))+geom_point(),
    ncol=2)
  })
  
  
  dataReady <- reactive({
    
    req(input$file.xlsxs)
    
    ## Access raw data ####
    dat <- dataXLSX()
    
    ## Access mAb dictionary ####
    mAb.dict <- datamAb()
    mAb.dict <- mAb.dict %>% arrange(mAb_plotting_order)
    
    ## Access VOC dictionary ####
    voc.dict <- datavoc()
    voc.dict <- voc.dict %>% arrange(VOC_plotting_order)
    
    ## Access SPC ####
    spc.data <- dataSPC()
    spc.data <- spc.data %>% filter(!is.na(average))
    ## ensure numeric
    spc.data <- spc.data %>% mutate(average = as.numeric(average)) %>%
      mutate(sd = as.numeric(sd))
    ## make sure spc.data$mAb is a factor and ordered as desired
    ## [calling this via geom_line + geom_ribbon -> fortified to a factor...
    ## and therefore alphabetical]
    spc.data <- spc.data %>%
      filter(mAb %in% mAb.dict$mAb_name_to_use)
    spc.data <- spc.data %>%
      mutate(mAb=factor(mAb, levels = (mAb.dict$mAb_name_to_use)))
    
    ## label monoclonals ####
    dat$mAb <- dat$Well.Content
    
    ## a for loop is acceptably quick to add labels to mAbs ##
    for(i in 1:nrow(mAb.dict)){
      dat[dat$mAb==mAb.dict$mAb_name_in_raw_data[i],]$mAb <-mAb.dict$mAb_name_to_use[i]
    }
    
    ## remove monoclonals not listed ####
    dat <- dat %>%
      filter(mAb %in% mAb.dict$mAb_name_to_use)
    
    ## fix factor levels
    mAb.dict <- mAb.dict %>% arrange(mAb_plotting_order)
    dat <- dat %>%
      mutate(mAb = factor(mAb, levels = mAb.dict$mAb_name_to_use))
    
    #--------------------------------------------------------------------#
    ## filter wells that fail QC ####
    dat <- dat %>% filter(well.fail != 'fail')
    #--------------------------------------------------------------------#
    ## rename variants and set factor order ####
    
    ## a for loop is acceptably quick ##
    for(i in 1:nrow(voc.dict)){
      dat[dat$variant==voc.dict$VOC_name_in_raw_data[i],]$variant <-
        voc.dict$VOC_name_to_use[i]
    }
    
    ## remove variants not listed ####
    dat <- dat %>%
      filter(variant %in% voc.dict$VOC_name_to_use)
    
    ## fix factor levels
    dat <- dat %>%
      mutate(variant = factor(variant, levels = voc.dict$VOC_name_to_use))
    
    readydat <- list()
    readydat$mAb.dict <- mAb.dict
    readydat$spc.data <- spc.data
    readydat$voc.dict <- voc.dict
    readydat$dat <-dat
    return(readydat)
  })
  
  output$quickPlot <- renderPlot({
    
    readydat <- dataReady()
    dat <- readydat$dat
    #--------------------------------------------------------------------#
    ## quick plot ####
    dat %>% #filter(str_detect(mAb, "NbD2")) %>%
      mutate(well.content.concentration = 
               case_when(well.content.concentration==0 ~ 10^-2,
                         well.content.concentration>0 ~ well.content.concentration)) %>%
      ggplot( aes(x=log10(well.content.concentration), 
                  y = Percentage_infected, 
                  col=variant)) + 
      geom_point(alpha=0.3, shape=20, size=2) + 
      ylim(-40, 200) + 
      geom_smooth() +
      khroma::scale_color_muted() +
      facet_grid(variant~mAb) +
      labs(title ="Quick plot, geom_smooth curves") +
      theme_bw()
    
    
  }, width = 1200, height= 1000)
  
  
  
  dataDR <- reactive({
    
    readydat <- dataReady()
    dat <- readydat$dat
    
    #--------------------------------------------------------------------#
    ## 4 parameter curve fits ####
    #--------------------------------------------------------------------#
    #--------------------------------------------------------------------#
    ## Set up data to use for prediction ####
    ## note this is mg/mL
    ## (and gets * 1000 later to go to ng/mL, which then logs to mostly +ve values)
    newdata <- expand.grid(well.content.concentration=10^seq(-5,5, length=100)) 
    
    ## Set up a prediction function
    pred <- function(x,  ...){
      z <- predict(x, newdata,  interval="confidence")
      z <- as_tibble(z)
      z <- z %>% mutate(well.content.concentration=newdata$well.content.concentration)
      return(z)
    }
    
    ## Set up a ECx function (where ECx = EC10, EC50 etc)
    edX <- function(x, n, level, ...){
      z <-ED(x, c(n), interval="tfls", display = F, level=level,
             lref=0.0001,uref=100,...)
      z <- as_tibble(z)
      return(z)
    }
    ## DRC fits and predictions #####
    ## Nest by monoclonal, so can purrr::map the DRC fits
    nested <- dat %>%
      nest(data = -c(variant,mAb))
    ## Map DRC fit + map predictions
    set.seed(42)
    results <- nested %>%
      mutate(drc = map(data, ~ drm(Percentage_infected ~ well.content.concentration,
                                   data = .x,
                                   fct = LL.4(), lowerl=c(
                                     # set slope to >0 (so flat lines are not fitted)
                                     0.1,
                                     
                                     0,0,0),
                                   upperl=c(
                                     # set slope < 1.5
                                     1.5,
                                     # set upper bound for bottom of curve = 50
                                     # (this stops a curve being fitted where all dilutions are 100%)
                                     # reduced to 20, to stop bottom of curves being fitted @ 50
                                     #
                                     20,
                                     # upper bound for response = 110
                                     110,NA)
      ))) %>% 
      mutate(ed50 = map(drc, ~ edX(.x,n=50,level=0.95))) %>%
      mutate(ed10 = map(drc, ~ edX(.x,n=10,level=0.95))) %>%
      mutate(ed20 = map(drc, ~ edX(.x,n=20,level=0.95))) %>%
      mutate(ed80 = map(drc, ~ edX(.x,n=80,level=0.95))) %>%
      mutate(ed90 = map(drc, ~ edX(.x,n=90,level=0.95))) %>%
      mutate(p = map(drc, ~ pred(.x, newdata=newdata))) %>%
      mutate(tidied = map(drc, tidy))
    #--------------------------------------------------------------------#
    
    #--------------------------------------------------------------------#
    ## Move negatives off the y axis by adding tiny amount to 0 ###
    dat <- dat %>% mutate(well.content.concentration = 
                            case_when(well.content.concentration==0 ~ 10^-3.4,
                                      well.content.concentration>0 ~ well.content.concentration))
    
    out <- list()
    out$dat <- dat
    out$results <- results
    return(out)
  })
  
  DRPlotter <- reactive({
    
    out <- dataDR()
    results <- out$results
    dat <- out$dat
    #--------------------------------------------------------------------#
    readydat <- dataReady()
    spc.data <- readydat$spc.data
    
    #--------------------------------------------------------------------#
    ## PLOTS ####
    DR <- results %>% 
      dplyr::select(variant,mAb,p) %>%
      unnest(p) %>%
      ## well.content.concentration is mcg/mL ##
      ## easier reading axis in ng/mL ##
      ## -> all x axis calls are * 1000.
      ggplot(aes(x=well.content.concentration*1000, y=Prediction, col = variant)) +
      geom_vline(data = results %>% 
                   dplyr::select(variant,mAb,ed50) %>%
                   unnest(ed50) %>%
                   ## do not plot extreme EC50s, with low conf ##
                   filter(log10(Estimate*1000) < 4.5),
                 aes(xintercept=log10(Estimate*1000), col=variant),
                 show.legend = F) +
      geom_pointrange(data = spc.data %>%
                        filter(timepoint=="cmax") %>%
                        dplyr::select(mAb, average,sd),
                      aes(xmin = log10((average-sd-sd)*1000),
                          xmax = log10((average+sd+sd)*1000),
                          x=log10(average*1000),
                          y= -10), col = "darkgrey",
                      fatten = 1,
                      inherit.aes = F,
                      show.legend = F) +
      # geom_rect(data = spc.data %>%
      #             filter(timepoint=="cmax") %>%
      #             dplyr::select(mAb, average,sd),
      #           aes(xmin = log10((average-sd-sd)*1000),
      #               xmax = log10((average+sd+sd)*1000),
      #               ymin=-Inf, ymax=Inf),fill="grey90",alpha=0.5, inherit.aes = F) +
      # geom_vline(data = spc.data %>%
      #              filter(timepoint=="c28") %>%
      #              dplyr::select(mAb, average),
      #            aes(xintercept = log10(average*1000)), linetype=3) +
      geom_pointrange(data = spc.data %>%
                   filter(timepoint=="c28") %>%
                   dplyr::select(mAb, average,sd),
                   aes(xmin = log10((average-sd-sd)*1000),
                       xmax = log10((average+sd+sd)*1000),
                       x=log10(average*1000),
                       y= -20),
                   fatten = 1.5,
                   inherit.aes = F) +
      # geom_vline(data = spc.data %>%
      #              filter(timepoint=="cmax") %>%
      #              dplyr::select(mAb, average),
      #            aes(xintercept = log10(average*1000)), linetype=2)+
      # geom_ribbon(aes(x=log10(well.content.concentration*1000),
      #                 y=Prediction, ymin=Lower, ymax=Upper, 
      #                 col=NULL,group=mAb), alpha=0.2, show.legend = F) + 
      geom_line(aes(x=log10(well.content.concentration*1000), y=Prediction)) +
      ## plot the raw data too:
      geom_point(data = dat,
                 aes(x=log10(well.content.concentration*1000), 
                     y=Percentage_infected, 
                     col = variant),
                 alpha=0.3, shape=20, size=1) +
      scale_color_muted(drop=F) + 
      scale_y_continuous(name= "% of maximal infection",
                         limits= c(-25, 150),
                         breaks = c(0,25,50,75,100),
                         labels = c(0,25,50,75,100)) + 
      scale_x_continuous(name= bquote(~log[10]~'[mAb]  ng/mL'),
                         limits = c(-1, 6),
                         breaks = c(-1,0,1,2,3,4,5),
                         labels = c(-1,0,1,2,3,4,5)) + 
      facet_grid(variant~mAb) +
      guides(colour = guide_legend(nrow = 1)) +
      theme_bw(base_family = "Helvetica Neue Thin") +
      theme(panel.grid.minor = element_blank(),
            axis.title.x = element_text(size=12),
            axis.text.x = element_text(size=8),
            axis.title.y = element_text(size=12),             strip.text = element_text(size=10),
            axis.text.y = element_text(size=8),
            strip.background = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal")
    
    return(DR)
    
  })
  
  output$DRPlot <- renderPlot({
    print(DRPlotter())
  }, width = 1200, height = 1000)
  
  output$downloadDRPlot <- downloadHandler(
    filename = function() { paste(input$dataset, '.svg', sep='') },
    content = function(file) {
      ggsave(file,DRPlotter(), height = 9, width = 8)
    }
  )

  
  InSilicoC28sera <- reactive({
    
    readydat <- dataReady()
    spc.data <- readydat$spc.data
    dat <- readydat$dat
    #--------------------------------------------------------------------#
    ## Set up data to use for prediction ####
    newdata <- spc.data %>%
      filter(timepoint == "c28") %>%
      mutate(dil1_well.content.concentration = average/40) %>%
      mutate(dil2_well.content.concentration = average/160) %>%
      mutate(dil3_well.content.concentration = average/640) %>%
      mutate(dil4_well.content.concentration = average/2560) %>%
      dplyr::select(mAb, dplyr::contains("dil")) %>%
      pivot_longer(cols = dplyr::contains("dil"), names_to = "dil", values_to = "well.content.concentration") %>%
      dplyr::select(well.content.concentration) %>% 
      as.data.frame(.)
    
    
    ## Set up a prediction function
    pred <- function(x,  ...){
      z <- predict(x, newdata,  interval="confidence")
      z <- as_tibble(z)
      z <- z %>% mutate(well.content.concentration=newdata$well.content.concentration)
      return(z)
    }
    ## DRC fits and predictions #####
    ## Nest by monoclonal, so can purrr::map the DRC fits
    nested <- dat %>%
      nest(data = -c(variant,mAb))
    ## Map DRC fit + map predictions
    set.seed(42)
    results <- nested %>%
      mutate(drc = map(data, ~ drm(Percentage_infected ~ well.content.concentration,
                                   data = .x,
                                   fct = LL.4(), lowerl=c(
                                     # set slope to >0 (so flat lines are not fitted)
                                     0.1,
                                     
                                     0,0,0),
                                   upperl=c(
                                     # set slope < 1.5
                                     1.5,
                                     # set upper bound for bottom of curve = 50
                                     # (this stops a curve being fitted where all dilutions are 100%)
                                     # reduced to 20, to stop bottom of curves being fitted @ 50
                                     #
                                     20,
                                     # upper bound for response = 110
                                     110,NA)
      ))) %>% 
      mutate(p = map(drc, ~ pred(.x, newdata=newdata))) %>%
      mutate(tidied = map(drc, tidy))
    
    outInSilico <- list()
    outInSilico$results <- results
    outInSilico$spc.data <- spc.data
    return(outInSilico)
    
  })
  
  
  EC50Plotter <- reactive({
    
    out <- dataDR()
    results <- out$results
    dat <- out$dat
    #--------------------------------------------------------------------#
    readydat <- dataReady()
    spc.data <- readydat$spc.data
    
    #--------------------------------------------------------------------#
    ## PLOTS ####
    ec50 <- results %>%
      dplyr::select(variant,mAb,ed50) %>%
      unnest(ed50) %>%
      mutate(variantForX = factor(variant, levels = rev(levels(variant)))) %>%
      mutate(Lower = ifelse(Lower<0, NA, Lower)) %>%
      filter(log10(Estimate*1000)<4) %>%
      ggplot(aes(x=variantForX, y=log10(Estimate*1000),
                 ymin=log10(Lower*1000),
                 ymax=log10(Upper*1000), col=variant)) +
      geom_errorbar() +
      geom_point(shape=20, alpha=1, size=2) +
      scale_color_muted(drop=F) + 
      scale_y_continuous(name= bquote(~log[10]~'[EC50]  ng/mL'),
                         limits = c(-1, 6),
                         breaks = c(-1,0,1,2,3,4,5),
                         labels = c(-1,0,1,2,3,4,5)) +
      facet_wrap(mAb~.,nrow = 1) +
      # labs(title ="EC50 with 95% CI") +
      scale_color_muted(drop=F) + 
      geom_pointrange(data = spc.data %>%
                        filter(timepoint=="cmax") %>%
                        dplyr::select(mAb, average,sd),
                      aes(ymin = log10((average-sd-sd)*1000),
                          ymax = log10((average+sd+sd)*1000),
                          y=log10(average*1000),
                          x= 1.2), col = "darkgrey",
                      fatten = 1.5,
                      inherit.aes = F,
                      show.legend = F) +
      geom_pointrange(data = spc.data %>%
                        filter(timepoint=="c28") %>%
                        dplyr::select(mAb, average,sd),
                      aes(ymin = log10((average-sd-sd)*1000),
                          ymax = log10((average+sd+sd)*1000),
                          y=log10(average*1000),
                          x= 0.8),
                      fatten = 1.5,
                      inherit.aes = F) +
      # geom_rect(data = spc.data %>%
      #             filter(timepoint=="cmax") %>%
      #             dplyr::select(mAb, average,sd),
      #           aes(ymin = log10((average-sd-sd)*1000),
      #               ymax = log10((average+sd+sd)*1000),
      #               xmin=-Inf, xmax=Inf),fill="grey90",alpha=0.5, inherit.aes = F) +
      # geom_hline(data = spc.data %>%
      #              filter(timepoint=="c28") %>%
      #              dplyr::select(mAb, average),
      #            aes(yintercept = log10(average*1000)), linetype=3) +
      # geom_hline(data = spc.data %>%
      #              filter(timepoint=="cmax") %>%
      #              dplyr::select(mAb, average),
      #            aes(yintercept = log10(average*1000)), linetype=2)+
      theme_bw(base_family = "Helvetica Neue Thin") +
      labs(x=NULL) + 
      theme(legend.position = "none",
            strip.background = element_blank(),
            panel.grid.minor = element_blank(),axis.title.x = element_text(size=12),
            axis.text.x = element_text(size=8),
            axis.title.y = element_text(size=12),
            strip.text = element_text(size=10),
            axis.text.y = element_text(size=8), strip.text.y = element_text(angle=0)) +
      coord_flip()# +
      # scale_x_discrete(breaks =
      #                    rev(ec50vocs),
      #                  limits = rev(ec50vocs), name = "")
    ec50
    
  })
  
  output$EC50Plot <- renderPlot({
    print(EC50Plotter())
  }, width = 1200, height = 300)
  
  output$downloadEC50Plot <- downloadHandler(
    filename = function() { paste(input$dataset, '.svg', sep='') },
    content = function(file) {
      ggsave(file,EC50Plotter(),  height = 4.5, width = 8)
    }
  )
  
  
  EC20Plotter <- reactive({
    
    out <- dataDR()
    results <- out$results
    dat <- out$dat
    #--------------------------------------------------------------------#
    readydat <- dataReady()
    spc.data <- readydat$spc.data
    
    #--------------------------------------------------------------------#
    ## PLOTS ####
    ec20 <- results %>%
      dplyr::select(variant,mAb,ed20) %>%
      unnest(ed20) %>%
      mutate(variantForX = factor(variant, levels = rev(levels(variant)))) %>%
      mutate(Lower = ifelse(Lower<0, NA, Lower)) %>%
      filter(log10(Estimate*1000)<4) %>%
      ggplot(aes(x=variantForX, y=log10(Estimate*1000),
                 ymin=log10(Lower*1000),
                 ymax=log10(Upper*1000), col=variant)) +
      geom_errorbar() +
      geom_point(shape=20, alpha=1, size=2) +
      scale_color_muted(drop=F) + 
      scale_y_continuous(name= bquote(~log[10]~'[EC20]  ng/mL'),
                         limits = c(-1, 6),
                         breaks = c(-1,0,1,2,3,4,5),
                         labels = c(-1,0,1,2,3,4,5)) +
      facet_wrap(mAb~.,nrow = 1) +
      # labs(title ="EC50 with 95% CI") +
      scale_color_muted(drop=F) + 
      geom_pointrange(data = spc.data %>%
                        filter(timepoint=="cmax") %>%
                        dplyr::select(mAb, average,sd),
                      aes(ymin = log10((average-sd-sd)*1000),
                          ymax = log10((average+sd+sd)*1000),
                          y=log10(average*1000),
                          x= 1.2), col = "darkgrey",
                      fatten = 1.5,
                      inherit.aes = F,
                      show.legend = F) +
      geom_pointrange(data = spc.data %>%
                        filter(timepoint=="c28") %>%
                        dplyr::select(mAb, average,sd),
                      aes(ymin = log10((average-sd-sd)*1000),
                          ymax = log10((average+sd+sd)*1000),
                          y=log10(average*1000),
                          x= 0.8),
                      fatten = 1.5,
                      inherit.aes = F) +
      # geom_rect(data = spc.data %>%
      #             filter(timepoint=="cmax") %>%
      #             dplyr::select(mAb, average,sd),
      #           aes(ymin = log10((average-sd-sd)*1000),
      #               ymax = log10((average+sd+sd)*1000),
      #               xmin=-Inf, xmax=Inf),fill="grey90",alpha=0.5, inherit.aes = F) +
      # geom_hline(data = spc.data %>%
      #              filter(timepoint=="c28") %>%
      #              dplyr::select(mAb, average),
      #            aes(yintercept = log10(average*1000)), linetype=3) +
      # geom_hline(data = spc.data %>%
      #              filter(timepoint=="cmax") %>%
      #              dplyr::select(mAb, average),
      #            aes(yintercept = log10(average*1000)), linetype=2)+
      theme_bw(base_family = "Helvetica Neue Thin") +
      labs(x=NULL) + 
      theme(legend.position = "none",
            strip.background = element_blank(),
            panel.grid.minor = element_blank(),axis.title.x = element_text(size=12),
            axis.text.x = element_text(size=8),
            axis.title.y = element_text(size=12),             strip.text = element_text(size=10),
            axis.text.y = element_text(size=8), strip.text.y = element_text(angle=0)) +
      coord_flip()# +
    # scale_x_discrete(breaks =
    #                    rev(ec50vocs),
    #                  limits = rev(ec50vocs), name = "")
    ec20
    
  })
  
  output$EC20Plot <- renderPlot({
    print(EC20Plotter())
  }, width = 1200, height = 300)
  
  output$downloadEC20Plot <- downloadHandler(
    filename = function() { paste(input$dataset, '.svg', sep='') },
    content = function(file) {
      ggsave(file,EC20Plotter(), height = 4.5, width = 8)
    }
  )
  
  EC80Plotter <- reactive({
    
    out <- dataDR()
    results <- out$results
    dat <- out$dat
    #--------------------------------------------------------------------#
    readydat <- dataReady()
    spc.data <- readydat$spc.data
    
    #--------------------------------------------------------------------#
    ## PLOTS ####
    ec80 <- results %>%
      dplyr::select(variant,mAb,ed80) %>%
      unnest(ed80) %>%
      mutate(variantForX = factor(variant, levels = rev(levels(variant)))) %>%
      mutate(Lower = ifelse(Lower<0, NA, Lower)) %>%
      filter(log10(Estimate*1000)<4) %>%
      ggplot(aes(x=variantForX, y=log10(Estimate*1000),
                 ymin=log10(Lower*1000),
                 ymax=log10(Upper*1000), col=variant)) +
      geom_errorbar() +
      geom_point(shape=20, alpha=1, size=2) +
      scale_color_muted(drop=F) + 
      scale_y_continuous(name= bquote(~log[10]~'[EC80]  ng/mL'),
                         limits = c(-1, 6),
                         breaks = c(-1,0,1,2,3,4,5),
                         labels = c(-1,0,1,2,3,4,5)) +
      facet_wrap(mAb~.,nrow = 1) +
      # labs(title ="EC50 with 95% CI") +
      scale_color_muted(drop=F) + 
      geom_pointrange(data = spc.data %>%
                        filter(timepoint=="cmax") %>%
                        dplyr::select(mAb, average,sd),
                      aes(ymin = log10((average-sd-sd)*1000),
                          ymax = log10((average+sd+sd)*1000),
                          y=log10(average*1000),
                          x= 1.2), col = "darkgrey",
                      fatten = 1.5,
                      inherit.aes = F,
                      show.legend = F) +
      geom_pointrange(data = spc.data %>%
                        filter(timepoint=="c28") %>%
                        dplyr::select(mAb, average,sd),
                      aes(ymin = log10((average-sd-sd)*1000),
                          ymax = log10((average+sd+sd)*1000),
                          y=log10(average*1000),
                          x= 0.8),
                      fatten = 1.5,
                      inherit.aes = F) +
      # geom_rect(data = spc.data %>%
      #             filter(timepoint=="cmax") %>%
      #             dplyr::select(mAb, average,sd),
      #           aes(ymin = log10((average-sd-sd)*1000),
      #               ymax = log10((average+sd+sd)*1000),
      #               xmin=-Inf, xmax=Inf),fill="grey90",alpha=0.5, inherit.aes = F) +
      # geom_hline(data = spc.data %>%
      #              filter(timepoint=="c28") %>%
      #              dplyr::select(mAb, average),
      #            aes(yintercept = log10(average*1000)), linetype=3) +
      # geom_hline(data = spc.data %>%
      #              filter(timepoint=="cmax") %>%
      #              dplyr::select(mAb, average),
      #            aes(yintercept = log10(average*1000)), linetype=2)+
      theme_bw(base_family = "Helvetica Neue Thin") +
      labs(x=NULL) + 
      theme(legend.position = "none",
            strip.background = element_blank(),
            panel.grid.minor = element_blank(),axis.title.x = element_text(size=12),
            axis.text.x = element_text(size=8),
            axis.title.y = element_text(size=12),             strip.text = element_text(size=10),
            axis.text.y = element_text(size=8), strip.text.y = element_text(angle=0)) +
      coord_flip()# +
    # scale_x_discrete(breaks =
    #                    rev(ec50vocs),
    #                  limits = rev(ec50vocs), name = "")
    ec80
    
  })
  
  
  EC90Plotter <- reactive({
    
    out <- dataDR()
    results <- out$results
    dat <- out$dat
    #--------------------------------------------------------------------#
    readydat <- dataReady()
    spc.data <- readydat$spc.data
    
    #--------------------------------------------------------------------#
    ## PLOTS ####
    ec90 <- results %>%
      dplyr::select(variant,mAb,ed90) %>%
      unnest(ed90) %>%
      mutate(variantForX = factor(variant, levels = rev(levels(variant)))) %>%
      mutate(Lower = ifelse(Lower<0, NA, Lower)) %>%
      filter(log10(Estimate*1000)<4) %>%
      ggplot(aes(x=variantForX, y=log10(Estimate*1000),
                 ymin=log10(Lower*1000),
                 ymax=log10(Upper*1000), col=variant)) +
      geom_errorbar() +
      geom_point(shape=20, alpha=1, size=2) +
      scale_color_muted(drop=F) + 
      scale_y_continuous(name= bquote(~log[10]~'[EC80]  ng/mL'),
                         limits = c(-1, 6),
                         breaks = c(-1,0,1,2,3,4,5),
                         labels = c(-1,0,1,2,3,4,5)) +
      facet_wrap(mAb~.,nrow = 1) +
      # labs(title ="EC50 with 95% CI") +
      scale_color_muted(drop=F) + 
      geom_pointrange(data = spc.data %>%
                        filter(timepoint=="cmax") %>%
                        dplyr::select(mAb, average,sd),
                      aes(ymin = log10((average-sd-sd)*1000),
                          ymax = log10((average+sd+sd)*1000),
                          y=log10(average*1000),
                          x= 1.2), col = "darkgrey",
                      fatten = 1.5,
                      inherit.aes = F,
                      show.legend = F) +
      geom_pointrange(data = spc.data %>%
                        filter(timepoint=="c28") %>%
                        dplyr::select(mAb, average,sd),
                      aes(ymin = log10((average-sd-sd)*1000),
                          ymax = log10((average+sd+sd)*1000),
                          y=log10(average*1000),
                          x= 0.8),
                      fatten = 1.5,
                      inherit.aes = F) +
      # geom_rect(data = spc.data %>%
      #             filter(timepoint=="cmax") %>%
      #             dplyr::select(mAb, average,sd),
      #           aes(ymin = log10((average-sd-sd)*1000),
      #               ymax = log10((average+sd+sd)*1000),
      #               xmin=-Inf, xmax=Inf),fill="grey90",alpha=0.5, inherit.aes = F) +
      # geom_hline(data = spc.data %>%
      #              filter(timepoint=="c28") %>%
      #              dplyr::select(mAb, average),
      #            aes(yintercept = log10(average*1000)), linetype=3) +
      # geom_hline(data = spc.data %>%
      #              filter(timepoint=="cmax") %>%
      #              dplyr::select(mAb, average),
      #            aes(yintercept = log10(average*1000)), linetype=2)+
      theme_bw(base_family = "Helvetica Neue Thin") +
      labs(x=NULL) + 
      theme(legend.position = "none",
            strip.background = element_blank(),
            panel.grid.minor = element_blank(),axis.title.x = element_text(size=12),
            axis.text.x = element_text(size=8),
            axis.title.y = element_text(size=12),             strip.text = element_text(size=10),
            axis.text.y = element_text(size=8), strip.text.y = element_text(angle=0)) +
      coord_flip()# +
    # scale_x_discrete(breaks =
    #                    rev(ec50vocs),
    #                  limits = rev(ec50vocs), name = "")
    ec90
    
  })
  
  output$EC80Plot <- renderPlot({
    print(EC80Plotter())
  }, width = 1200, height = 300)
  
  output$downloadEC80Plot <- downloadHandler(
    filename = function() { paste(input$dataset, '.svg', sep='') },
    content = function(file) {
      ggsave(file,EC80Plotter(), height = 4.5, width = 8)
    }
  )
  
  
  output$EC90Plot <- renderPlot({
    print(EC90Plotter())
  }, width = 1200, height = 300)
  
  output$downloadEC90Plot <- downloadHandler(
    filename = function() { paste(input$dataset, '.svg', sep='') },
    content = function(file) {
      ggsave(file,EC90Plotter(), height = 4.5, width = 8)
    }
  )
  
  
  
  
  output$downloadData <- downloadHandler(
    filename = function(){
      paste0("EC50_",lubridate::now(),".csv")
    },
    
    content = function(file) {
      out <- dataDR()
      results <- out$results
      
      results %>%
        dplyr::select(variant,mAb,ed50) %>%
        unnest(ed50) %>%
        mutate(across(!c(variant, mAb), ~ .x * 1000)) %>%
        mutate(across(!c(variant, mAb), ~ format(.x, digits=3,scientific=F,nsmall=1))) %>%
      write.csv(., file, row.names = FALSE)
    })

  
  output$downloadData80 <- downloadHandler(
    filename = function(){
      paste0("EC80_",lubridate::now(),".csv")
    },
    
    content = function(file) {
      
      out <- dataDR()
      results <- out$results
      
      results %>%
        dplyr::select(variant,mAb,ed80) %>%
        unnest(ed80) %>%
        mutate(across(!c(variant, mAb), ~ .x * 1000)) %>%
        mutate(across(!c(variant, mAb), ~ format(.x, digits=3,scientific=F,nsmall=1))) %>%
        write.csv(., file, row.names = FALSE)
    })



output$downloadData20 <- downloadHandler(
  filename = function(){
    paste0("EC20_",lubridate::now(),".csv")
  },
  
  content = function(file) {
    
    out <- dataDR()
    results <- out$results
    
    results %>%
      dplyr::select(variant,mAb,ed20) %>%
      unnest(ed20) %>%
      mutate(across(!c(variant, mAb), ~ .x * 1000)) %>%
      mutate(across(!c(variant, mAb), ~ format(.x, digits=3,scientific=F,nsmall=1))) %>%
      write.csv(., file, row.names = FALSE)
  })



output$downloadData90 <- downloadHandler(
  filename = function(){
    paste0("EC90_",lubridate::now(),".csv")
  },
  
  content = function(file) {
    
    out <- dataDR()
    results <- out$results
    
    results %>%
      dplyr::select(variant,mAb,ed90) %>%
      unnest(ed90) %>%
      mutate(across(!c(variant, mAb), ~ .x * 1000)) %>%
      mutate(across(!c(variant, mAb), ~ format(.x, digits=3,scientific=F,nsmall=1))) %>%
      write.csv(., file, row.names = FALSE)
  })



output$downloadDataInSilico <- downloadHandler(
  filename = function(){
    paste0("InSilicoC28_",lubridate::now(),".csv")
  },
  
  content = function(file) {
    
    outInSilico <- InSilicoC28sera()
    results <- outInSilico$results

    spc.data <- outInSilico$spc.data

    spc.data <- spc.data %>%
      dplyr::filter(timepoint == "c28") %>%
      dplyr::mutate(dil1_well.content.concentration = average/40) %>%
      dplyr::mutate(dil2_well.content.concentration = average/160) %>%
      dplyr::mutate(dil3_well.content.concentration = average/640) %>%
      dplyr::mutate(dil4_well.content.concentration = average/2560) %>%
      dplyr::select(mAb, dplyr::contains("dil")) %>%
      tidyr::pivot_longer(cols = dplyr::contains("dil"), names_to = "dil", values_to = "well.content.concentration") %>%
      dplyr::mutate(to_keep = "y") %>%
      dplyr::select(mAb, well.content.concentration, to_keep)

    results %>%
      dplyr::select(variant, mAb, p) %>% unnest(p) %>% 
      dplyr::left_join(spc.data) %>%
      dplyr::filter(to_keep=="y") %>%
      write.csv(., file, row.names = FALSE)
  })

}


#--------------------------------------------------------------------#

ui <- fluidPage(
  navbarPage("Dose:response curves",
             
             #--------
             tabPanel("mAb dictionary import",
                      # Sidebar panel for inputs ----
                      sidebarPanel(
                        fileInput("file.mAb", h3("mAb CSV dictionary"),
                                  multiple = FALSE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                # Input: Select quotes ----
                                radioButtons("quote", "Quote",
                                             choices = c(None = "",
                                                         "Double Quote" = '"',
                                                         "Single Quote" = "'"),
                                             selected = '"'),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Select number of rows to display ----
                                radioButtons("disp", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head")
                      ),
                      # Main panel for displaying outputs ----
                      mainPanel(
                        
                        # Output: Data file ----
                        tableOutput("contentsmAb")
                        
                      )
                      
             ),
             #--------
             
             #--------
             tabPanel("VOC dictionary import",
                      # Sidebar panel for inputs ----
                      sidebarPanel(
                        fileInput("file.voc", h3("VOC CSV dictionary"),
                                  multiple = FALSE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        # Horizontal line ----
                        tags$hr(),
                        
                        # Input: Checkbox if file has header ----
                        checkboxInput("header", "Header", TRUE),
                        
                        # Input: Select separator ----
                        radioButtons("sep", "Separator",
                                     choices = c(Comma = ",",
                                                 Semicolon = ";",
                                                 Tab = "\t"),
                                     selected = ","),
                        
                        # Input: Select quotes ----
                        radioButtons("quote", "Quote",
                                     choices = c(None = "",
                                                 "Double Quote" = '"',
                                                 "Single Quote" = "'"),
                                     selected = '"'),
                        
                        # Horizontal line ----
                        tags$hr(),
                        
                        # Input: Select number of rows to display ----
                        radioButtons("disp", "Display",
                                     choices = c(Head = "head",
                                                 All = "all"),
                                     selected = "head")
                      ),
                      # Main panel for displaying outputs ----
                      mainPanel(
                        
                        # Output: Data file ----
                        tableOutput("contentsvoc")
                        
                      )
                      
             ),
             #--------
             #--------
             tabPanel("SPC concentrations import",
                      # Sidebar panel for inputs ----
                      sidebarPanel(
                        fileInput("file.spc", h3("SPC concentrations CSV"),
                                  multiple = FALSE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        # Horizontal line ----
                        tags$hr(),
                        
                        # Input: Checkbox if file has header ----
                        checkboxInput("header", "Header", TRUE),
                        
                        # Input: Select separator ----
                        radioButtons("sep", "Separator",
                                     choices = c(Comma = ",",
                                                 Semicolon = ";",
                                                 Tab = "\t"),
                                     selected = ","),
                        
                        # Input: Select quotes ----
                        radioButtons("quote", "Quote",
                                     choices = c(None = "",
                                                 "Double Quote" = '"',
                                                 "Single Quote" = "'"),
                                     selected = '"'),
                        
                        # Horizontal line ----
                        tags$hr(),
                        
                        # Input: Select number of rows to display ----
                        radioButtons("disp", "Display",
                                     choices = c(Head = "head",
                                                 All = "all"),
                                     selected = "head")
                      ),
                      # Main panel for displaying outputs ----
                      mainPanel(
                        
                        # Output: Data file ----
                        tableOutput("contentsspc")
                        
                      )
                      
             ),
             #--------
             #--------
             tabPanel("Neutralisation data import",
                      # Sidebar panel for inputs ----
                      sidebarPanel(
                        fileInput("file.xlsxs", h3("Upload ALL xlsx files"),
                                  multiple = TRUE,
                                  accept = c(".xlsx")),
                        # Horizontal line ----
                        tags$hr()#,
                        
                        # # Input: Select number of rows to display ----
                        # radioButtons("disp", "Display",
                        #              choices = c(Head = "head",
                        #                          All = "all"),
                        #              selected = "head")
                      ),
                      # Main panel for displaying outputs ----
                      mainPanel(
                        
                        # Output: Data file ----
                        tableOutput("contentsxlsx")
                        
                      )
                      
             ),
             #--------
             tabPanel("QC plot",
                      # Main panel for displaying outputs ----
                      mainPanel(
                        
                        # Output: Data file ----
                        plotOutput("distPlot")
                        
                      )
                      
             ),
             #--------
             #--------
             tabPanel("Quick plot",
                      # Main panel for displaying outputs ----
                      mainPanel(
                        
                        # Output: Data file ----
                        plotOutput("quickPlot")
                        
                      )
                      
             ),
             #--------
             tabPanel("Dose-response plot",
                      sidebarPanel(downloadButton('downloadDRPlot', 'Download Plot')),
                      # Main panel for displaying outputs ----
                      mainPanel(
                        
                        # Output: Data file ----
                        plotOutput("DRPlot")
                        
                        
                      )
                      
             ),
             #--------
             #--------
             tabPanel("EC50 plot",
                      sidebarPanel(downloadButton('downloadEC50Plot', 'Download Plot')),
                      # Main panel for displaying outputs ----
                      mainPanel(
                        
                        # Output: Data file ----
                        plotOutput("EC50Plot")
                        
                      )
                      
             ),
             #--------
             tabPanel("EC80 plot",
                      sidebarPanel(downloadButton('downloadEC80Plot', 'Download Plot')),
                      # Main panel for displaying outputs ----
                      mainPanel(
                        
                        # Output: Data file ----
                        plotOutput("EC80Plot")
                        
                      )
                      
             ),
             #--------
             tabPanel("EC20 plot",
                      sidebarPanel(downloadButton('downloadEC20Plot', 'Download Plot')),
                      # Main panel for displaying outputs ----
                      mainPanel(
                        
                        # Output: Data file ----
                        plotOutput("EC20Plot")
                        
                      )
                      
             ),
             #--------
             tabPanel("EC90 plot",
                      sidebarPanel(downloadButton('downloadEC90Plot', 'Download Plot')),
                      # Main panel for displaying outputs ----
                      mainPanel(
                        
                        # Output: Data file ----
                        plotOutput("EC90Plot")
                        
                      )
                      
             ),
             #--------
             tabPanel("Write out ECx CSV",
                      downloadButton('downloadData', 'Download EC50'),
                      downloadButton('downloadData80', 'Download EC80'),
                      downloadButton('downloadData20', 'Download EC20'),
                      downloadButton('downloadData90', 'Download EC90'),
                      downloadButton('downloadDataInSilico', "Download InSilico EC50s for SPC-C28")
                      )
                      
             )
             #--------

                                
                      # )#,
                      # fileInput("file.VOC", h3("VOC dictionary"),
                      #           multiple = FALSE,
                      #           accept = c("text/csv",
                      #                      "text/comma-separated-values,text/plain",
                      #                      ".csv")),
                      # fileInput("file.spc", h3("SPC concentrations"),
                      #           multiple = FALSE,
                      #           accept = c("text/csv",
                      #                      "text/comma-separated-values,text/plain",
                      #                      ".csv")),
                      # ),
             # tabPanel("XLSX import"#,
                      # sidebarLayout(
                      #   sidebarPanel(
                      #     sliderInput("obs", "Number of observations:", min = 10, max = 500, value = 100)
                      #   ),
                      #   mainPanel(plotOutput("distPlot"))
                      # )
             #          ),
             # tabPanel("Quick plot",
             #          mainPanel(plotOutput("distPlot"))),
             # tabPanel("DR plot"),
             # tabPanel("EC50 table")
  # )
)

shinyApp(ui = ui, server = server)
