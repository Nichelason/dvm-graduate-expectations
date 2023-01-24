###################################################################################################
#####  Code written by Scott Hetzel MS, biostatistician at University of Wisconsin - Madison  #####
###################################################################################################
#####  We acknowledge the intellectual and technical contributions of Scott Hetzel of the
#####  Biostatistics and Epidemiology Research Design Core to the development of this code
#####  for manuscript publication. This work was funded by Institutional Clinical and
#####  Translational Science Award UL1 TR002373.
###################################################################################################
library(nlme)

# LOAD the Analysis Data Set - you will need to update the file path
dat <- read.csv("AnalysisDataSet.csv",header=TRUE,
                stringsAsFactors=TRUE)

# Add columns for grouping schemes
dat$Companion <- factor(ifelse(dat$Emphasis=="Companion Animal","Companion Animal",
                               ifelse(dat$Emphasis=="Non-Technical",NA,"Non-Companion Animal")))

# SumFunc code
SumFunc <- function(data,
                    completed="True",
                    task=NULL,
                    group_by=NULL,
                    comparison=c("left","right","bottom","bottom-left","bottom-right"),
                    time=c("Pre","Post"),
                    numDec=2,
                    pvalue=FALSE,
                    pTest=c("t-test","wilcoxon","lme"),
                    pAdjust=c("none","holm","bonferroni"),
                    printTable=TRUE,
                    printData=FALSE,
                    csvFileName=NULL)
{
  ###########  FORMATTING CONFIRMATIONS  #################################
  # make sure completed is either "True" or "False"
  if(!(completed %in% c("True","False")))
    stop("'completed' must be either 'True' or 'False'")
  # make sure task is one of 'Technical' or 'Non-Technical'
  if(!is.null(task))
  {
    if(!(task[1] %in% c("Technical","Non-Technical")))
      stop("'task' must be one of 'Technical' or 'Non-Technical'")
    task <- task[1]
  }
  # make sure group_by is a column in data
  if(!is.null(group_by))
  {
    if(!(group_by[1] %in% colnames(data)))
      stop("'group_by' must be a column in data")
    group_by <- group_by[1]
  }
  # make sure comparison is one of the possible inputs
  if(length(comparison)>1)
  {
    print("First entry of 'comparison' was used")
    comparison <- comparison[1]
  }
  if(!(comparison %in% c("left","right","bottom","bottom-left","bottom-right")))
    stop("'comparison' must be one of 'left', 'right', 'bottom', 'bottom-left', or 'bottom-right'")
  # check if 'time' is needed and whether the value is valid
  if(comparison %in% c("bottom-left","bottom-right"))
    time <- NULL
  else
  {
    if(length(time) > 1)
    {
      print("First entry of 'time' was used")
      time <- time[1]
    }
    if(!(time %in% c("Pre","Post")))
      stop("'time' must be one of 'Pre' or 'Post'")
  }
  # check if pvalue is logical
  if(!is.logical(pvalue))
    stop("'pvalue' must be TRUE or FALSE")
  # check if pTest is needed and reduce to first letter for determining Test method
  if(pvalue)
  {
    if(is.null(pTest))
      pTest <- "t"
    else
      pTest <- substr(pTest[1],1,1)
    if(!(pTest %in% c("t","w","l")))
      stop("'pTest' must length of 1 and first entry must start with 't', 'w', or 'l'")
  }
  ###############################################################################
  ############  DATA PREP  ######################################################
  # subset data based on 'completed'
  if(completed=="True")
    data <- data[data$Finished=="True",]
  # check task and subset if needed
  if(!is.null(task))
    data <- data[data$Task %in% task,]
  # subset based on comparison
  if(comparison == "left")  # 'left' means WVMA expectations vs Student Pre or Post expectations
    data <- data[data$Type=="Expectation" & data$Group %in% c("WVMA",paste("Student",time)),]
  else if(comparison == "right") # right means WVMA expectations vs Student Pre or Post abilities
    data <- data[data$Group=="WVMA" | (data$Group==paste("Student",time) & data$Type=="Ability"),]
  else if(comparison == "bottom") # bottom means Student Expectations Pre/Post vs Ability Pre/Post
    data <- data[data$Group==paste("Student",time),]
  else if(comparison == "bottom-left") # bottom-left means Student Expectations Pre vs Post
    data <- data[data$Type=="Expectation" & data$Group!="WVMA",]
  else if(comparison == "bottom-right") # bottom-right means Student Ability Pre vs Post
    data <- data[data$Type=="Ability" & data$Group!="WVMA",]
  # remove any data rows that are NA in group_by
  data <- data[!is.na(data[[group_by]]),]
  # remove any data rows that are NA in Response
  data <- data[!is.na(data$Response),]
  # re-factor pertinent variables to remove any blank levels caused by the subsetting
  c.names <- c("RecipientEmail","VarName","Group","Task","Emphasis","Category","Type","Question")
  for(i in 1:length(c.names))
    data[[c.names[i]]] <- factor(data[[c.names[i]]])
  
  ###############################################################################
  ###########  TABLE CREATION  ##################################################
  mat <- matrix("",nrow=(nlevels(data[[group_by]])+2),ncol=8)
  levs <- factor(data$Group:data$Type)
  data$combn <- data[[group_by]]:levs
  colnames(mat) <- c("","","",levels(levs)[1],"","",levels(levs)[2],"")
  mat[1,] <- c("Group","Questions",rep(c("N","Mean (SD)","Median (Range)"),2))
  mat[2:(nrow(mat)-1),1] <- levels(data[[group_by]])
  mat[2:(nrow(mat)-1),2] <- tapply(data$VarName,data[[group_by]],function(x){length(unique(x))})
  # if WVMA is involved in the subset of data (i.e. comparison = left or right) the ResponseId must be used
  # for the response grouping to calculate N.  If WVMA is not involved then RecipientEmail must be used
  if(!(comparison %in% c("left","right")))
    mat[2:(nrow(mat)-1),c(3,6)] <- matrix(tapply(data$RecipientEmail,data$combn,function(x){length(unique(x))}),
                                          ncol=2,byrow=TRUE)
  else
    mat[2:(nrow(mat)-1),c(3,6)] <- matrix(tapply(data$ResponseId,data$combn,function(x){length(unique(x))}),
                                          ncol=2,byrow=TRUE)
  
  meanVal <- list()
  medianVal <- list()
  means <- c()
  sds <- c()
  medians <- c()
  mins <- c()
  maxes <- c()
  for(j in 1:nlevels(data$combn))
  { # if WVMA is involved in the subset of data (i.e. comparison = left or right) the ResponseId must be used
    # for the response grouping to calculate mean and median.  If WVMA is not involved then RecipientEmail must be used
    if(!(comparison %in% c("left","right")))
    {
      meanVal[[j]] <- tapply(data$Response[data$combn==levels(data$combn)[j]],
                             factor(data$RecipientEmail[data$combn==levels(data$combn)[j]]),mean)
      medianVal[[j]] <- tapply(data$Response[data$combn==levels(data$combn)[j]],
                               factor(data$RecipientEmail[data$combn==levels(data$combn)[j]]),mean)
    }
    else
    {
      meanVal[[j]] <- tapply(data$Response[data$combn==levels(data$combn)[j]],
                             factor(data$ResponseId[data$combn==levels(data$combn)[j]]),mean)
      medianVal[[j]] <- tapply(data$Response[data$combn==levels(data$combn)[j]],
                               factor(data$ResponseId[data$combn==levels(data$combn)[j]]),mean)
    }
    names(meanVal)[j] <- levels(data$combn)[j]
    names(medianVal)[j] <- levels(data$combn)[j]
    means[j] <- round(mean(meanVal[[j]]),numDec)
    sds[j] <- round(sd(meanVal[[j]]),numDec)
    medians[j] <- round(median(medianVal[[j]]),numDec)
    mins[j] <- round(min(medianVal[[j]]),numDec)
    maxes[j] <- round(max(medianVal[[j]]),1)
  }
  
  mat[2:(nrow(mat)-1),c(4,7)] <- matrix(paste0(means," (",sds,")"),ncol=2,byrow=TRUE)
  mat[2:(nrow(mat)-1),c(5,8)] <- matrix(paste0(medians," (",mins," - ",maxes,")"),ncol=2,byrow=TRUE)
  mat[nrow(mat),1] <- paste0("completed=",completed,", task=",task,", group_by=",group_by)
  if(pvalue)
  {
    pv <- c()
    if(pTest=="t")
      if(comparison == "bottom") # do a paired analysis
        for(k in 1:nlevels(data[[group_by]])) # this is based on knowing that meanVal list is stacked by group_by level
          pv[k] <- t.test(meanVal[[(k*2-1)]],meanVal[[(k*2)]],paired=TRUE)$p.value
    else
      for(k in 1:nlevels(data[[group_by]])) # this is based on knowing that meanVal list is stacked by group_by level
        pv[k] <- t.test(meanVal[[(k*2-1)]],meanVal[[(k*2)]])$p.value
      else if(pTest=="w")
        if(comparison == "bottom") # do a paired analysis
          for(k in 1:nlevels(data[[group_by]])) # this is based on knowing that meanVal list is stacked by group_by level
            pv[k] <- wilcox.test(medianVal[[(k*2-1)]],medianVal[[(k*2)]],paired=TRUE)$p.value
      else
        for(k in 1:nlevels(data[[group_by]]))
          pv[k] <- wilcox.test(medianVal[[(k*2-1)]],medianVal[[(k*2)]])$p.value
        else if(pTest=="l")
          for(k in 1:nlevels(data[[group_by]]))
          {
            resp <- as.numeric(c(medianVal[[(k*2-1)]],medianVal[[(k*2)]]))
            fact <- factor(rep(c(names(medianVal)[(k*2-1)],names(medianVal)[(k*2)]),
                               c(length(medianVal[[(k*2-1)]]),length(medianVal[[(k*2)]]))))
            ids <- as.character(c(names(medianVal[[(k*2-1)]]),names(medianVal[[(k*2)]])))
            anv <- anova(lme(resp ~ fact, random=~1|ids))
            pv[k] <- anv[2,4]
          }
        if(pAdjust %in% c("holm","bonferroni"))
          pv <- p.adjust(pv,method=pAdjust)
        pv <- ifelse(pv<0.0005,"< 0.001",round(pv,3))
        mat <- cbind(mat,c(paste(ifelse(pAdjust %in% c("holm","bonferroni"),"Adj",""),
                                 ifelse(pTest=="t","T-test p-value",
                                        ifelse(pTest=="w","Wilcoxon p-value","LME p-value"))),pv,""))
  }
  if(!is.null(csvFileName))
    write.csv(mat,file=csvFileName,row.names=FALSE)
  if(printTable)
    return(mat)
  if(printData)
    return(data)
}



##########  There are 2 examples of SumFunc runs below  #################        

#dat$Companion <- factor(ifelse(dat$Emphasis=="Companion Animal","Companion Animal",
#                               ifelse(dat$Emphasis=="Non-Technical",NA,"Non-Companion Animal")))

#mat1 <- SumFunc(data=dat,
#                completed="False",
#                task=NULL,
#                group_by="Companion",
#                comparison="bottom-right",
#                time=NULL,
#                numDec=1,
#                pvalue=TRUE,
#                pTest="l",
#                pAdjust="bonferroni",
#                printTable=TRUE,
#                printData=FALSE)

#tab2 <- SumFunc(data=dat,
#        completed="False",
#        task=NULL,
#        group_by="Task",
#        comparison="bottom-right",
#        time=NULL,
#        numDec=2,
#        pvalue=TRUE,
#        pTest="w",
#        pAdjust="bonferroni",
#        printTable=FALSE,
#        printData=TRUE,
#        csvFileName=NULL)

##########  Shiny App Setup  ################# 
library(shiny)
library(shinyjs)
library(DT)
library(htmltools)

# Get column names for groupings
grouping_cols <- colnames(dat)
grouping_cols <- grouping_cols[c(7:9, 11:ncol(dat))]

# Make into a more descriptive list of options
groupings_explained <- grouping_cols
groupings_explained[1] <- "Technical or Non-Technical"
groupings_explained[2] <- "Emphasis Area"
groupings_explained[5] <- "Companion Animal or Non-Companion Animal"

ui <- fluidPage(
  shinyjs::useShinyjs(),
  titlePanel("Graduate Expectations Explorer"),
  
  sidebarLayout(
    
    sidebarPanel(
      # completed
      checkboxInput("completed", "Include only respondents who completed a survey", FALSE),
      
      # task
      radioButtons("task", "Select item type", choices=c("Both", "Technical", "Non-Technical"), selected = "Both"),
      
      # group_by
      radioButtons("group_by", "Group items by", selected="Emphasis", choiceNames=groupings_explained, choiceValues=grouping_cols),
      
      # comparison
      radioButtons("comparison", "Compare", choiceNames=c("Student Expectations vs. Practitioner Expectations","Student Confidence vs. Practitioner Expectations","Student Expectations vs. Student Confidence","Student Expectations (Pre vs. Post)","Student Confidence (Pre vs. Post)"), choiceValues=c("left","right","bottom","bottom-left","bottom-right"), selected="left"),
      
      # time
      radioButtons("time", "Use pre or post? (For non-'Pre vs. Post' comparisons only.)", choices=c("Pre","Post"), selected="Post"),
      
      # numDec
      sliderInput("numDec", "Select number of decimals in output", min=1, max=5, value=2, step=1, round=TRUE, ticks=FALSE),
      
      # pvalue
      # select by adding 'none' as an option to pTest below
      
      # pTest
      radioButtons("pTest", "Select a statistical test", choices=c("none", "t-test","wilcoxon", "lme"), selected="none"),
      
      # pAdjust
      radioButtons("pAdjust", "Select a method of p-value correction", choices=c("none","holm","bonferroni"), selected="none"),
      
      # printTable
      # always print
      
      # Show result
      actionButton("doBtn", "Show result"),
      
      # Download button - only show if results are shown
      downloadButton("downloadData", label = "Download"),
    ),
    
    mainPanel(
      DTOutput("result"),
    )
  )
)

server <- function(input, output) {
  # Disable download button
  shinyjs::disable("downloadData")
  
  # Press button
  observeEvent(input$doBtn, {
    # Transform form inputs as needed for function inputs
    if(input$completed == TRUE) completed = "True" else completed = "False"
    if(input$task == "Both") task = NULL else task = input$task
    if(input$pTest == "none") pvalue = FALSE else pvalue = TRUE
    
    # Execute function
    mat <- SumFunc(dat,
                   completed=completed,
                   task=task,
                   group_by=input$group_by,
                   comparison=input$comparison,
                   time=input$time,
                   numDec=input$numDec,
                   pvalue=pvalue,
                   pTest=input$pTest,
                   pAdjust=input$pAdjust,
                   printTable=TRUE,
                   csvFileName="last_result.csv")
    
    # Create custom formatting for rendering table output
    sketch = htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(colspan = 2, ''),
          th(colspan = 3, colnames(mat)[4]),
          th(colspan = 3, colnames(mat)[7]),
          th(colspan = 2, '')
        ),
        tr(
          lapply(mat[1,], th)
        )
      )
    ))
    
    # Conditionally add last column if present
    if(ncol(mat) > 8){
      sketch = htmltools::withTags(table(
        class = 'display',
        thead(
          tr(
            th(colspan = 2, ''),
            th(colspan = 3, colnames(mat)[4]),
            th(colspan = 3, colnames(mat)[7]),
            th(colspan = 2, '') # include extra column for p value
          ),
          tr(
            lapply(mat[1,], th)
          )
        )
      ))
    } else {
      sketch = htmltools::withTags(table(
        class = 'display',
        thead(
          tr(
            th(colspan = 2, ''),
            th(colspan = 3, colnames(mat)[4]),
            th(colspan = 3, colnames(mat)[7]),
          ),
          tr(
            lapply(mat[1,], th)
          )
        )
      ))
    }
    
    
    # Populate table
    output$result <- DT::renderDataTable(mat[2:(nrow(mat)-1),], container = sketch, server=FALSE)
    show("result")
    
    # Attach file to download button
    output$downloadData <- downloadHandler(
      filename <- "last_result.csv",
      
      content <- function(file) {
        file.copy("last_result.csv", file)
      },
      contentType = "text/csv"
    )
    
    # Enable download button
    shinyjs::enable("downloadData")
  })
  
  # Hide table if input parameters are changed (anything but results button)
  observeEvent(ignoreInit = TRUE, c(
    input$completed,
    input$task,
    input$group_by,
    input$comparison,
    input$time,
    input$numDec,
    input$pTest,
    input$pAdjust
  ),
  {
    hide("result")
    shinyjs::disable("downloadData")
  })
  
}


shinyApp(ui = ui, server = server)
