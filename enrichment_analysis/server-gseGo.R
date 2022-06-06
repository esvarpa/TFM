
myValues = reactiveValues()

observe({
  gseGoReactive()
})

gseGoReactive <- eventReactive(input$initGo,{
  
  withProgress(message = "Processing , please wait",{
    
    isolate({
      
      # eliminar notificaciones si existen
      removeNotification("errorNotify")
      removeNotification("errorNotify1")
      removeNotification("errorNotify2")
      removeNotification("warnNotify")
      removeNotification("warnNotify2")
      
      shiny::validate(need(tryCatch({
        df <- inputDataReactive()$data
        # log2 fold change 
        original_gene_list <- df[[input$log2fcColumn]]
        
        # nombre del vector
        names(original_gene_list) <- df[[input$geneColumn]]
        
        # omitir valores Na 
        gene_list<-na.omit(original_gene_list)
        
        # ordenar la lista en orden decreciente (requerido para clusterProfiler)
        gene_list = sort(gene_list, decreasing = TRUE)
        
        myValues$gene_list = gene_list
        
        
        # 
        
        setProgress(value = 0.3, detail = "Performing GSE analysis, please wait ...")
        
        orgDb.obj = eval( parse(text = input$organismDb, keep.source=FALSE))
        
        go_gse <- gseGO(geneList=gene_list, 
                     ont = input$ontology, 
                     keyType = input$keytype, 
                     minGSSize = input$minGSSize, 
                     maxGSSize = input$maxGSSize, 
                     pvalueCutoff = 0.05, 
                     verbose = T, 
                     OrgDb = orgDb.obj, 
                     pAdjustMethod = input$pAdjustMethod
                     )
        
        if(nrow(go_gse) < 1)
        {
          showNotification(id="warnNotify", "No gene can be mapped ...", type = "warning", duration = NULL)
          showNotification(id="warnNotify2", "Tune the parameters and try again.", type = "warning", duration = NULL)
          return(NULL)
        }
        
        updateNumericInput(session, "showCategory_bar", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_dot", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_enrichmap", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_goplot", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_cnet", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
        
        updateNumericInput(session, "showCategory_bar_kegg", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_dot_kegg", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_enrichmap_kegg", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
        #updateNumericInput(session, "showCategory_goplot", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_cnet_kegg", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
        
        
        updateSelectizeInput(session,'pubmedTerms', choices=go_gse@result$Description)
        
        ## KEGG gse
        
        # Convierta IDs de genes para la funci?n gseKEGG
        # Se pierden alguna genes porque no todas las identificaciones se convierten
        myValues$convWarningMessage = capture.output(ids<-bitr(names(original_gene_list), fromType = input$keytype, toType = "ENTREZID", OrgDb=input$organismDb), type = "message")
        
        # elimine IDS duplicados (aqu? uso "ENSEMBL", pero deber?a ser lo que se seleccion? como keyType)
        dedup_ids = ids[!duplicated(ids[c(input$keytype)]),]
        
        # Cree un nuevo marco de datos df2 que tenga solo los genes que se mapearon con ?xito usando la funci?n bitr anterior
        #df2 = df[df$X %in% dedup_ids$ENSEMBL,]
        #df2 = df[df$X %in% dedup_ids[,1],]
        df2 = df[df[[input$geneColumn]] %in% dedup_ids[,1],]
        
        # Crear una nueva columna en df2 con los ID de ENTREZ correspondientes
        df2$Y = dedup_ids$ENTREZID
        
        # crear un vector del universo de genes
        kegg_gene_list <- df2[[input$log2fcColumn]]
        
        # nombre del vector con ENTREZ ids
        names(kegg_gene_list) <- df2$Y
        
        # omitir valores NA 
        kegg_gene_list<-na.omit(kegg_gene_list)
        
        # ordenar la lista en orden decreciente (requerido para clusterProfiler)
        kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
        
        
        setProgress(value = 0.6, detail = "Performing KEGG enrichment analysis, please wait ...")
        
        organismsDbKegg = c("org.Hs.eg.db"="hsa","org.Mm.eg.db"="mmu","org.Rn.eg.db"="rno",
                            "org.Sc.sgd.db"="sce","org.Dm.eg.db"="dme","org.At.tair.db"="ath",
                            "org.Dr.eg.db"="dre","org.Bt.eg.db"="bta","org.Ce.eg.db"="cel",
                            "org.Gg.eg.db"="gga","org.Cf.eg.db"="cfa","org.Ss.eg.db"="ssc",
                            "org.Mmu.eg.db"="mcc","org.EcK12.eg.db"="eck","org.Xl.eg.db"="xla",
                            "org.Pt.eg.db"="ptr","org.Ag.eg.db"="aga","org.Pf.plasmo.db"="pfa",
                            "org.EcSakai.eg.db"="ecs")
        
        kegg_gse <- gseKEGG(geneList=kegg_gene_list, 
                            organism=organismsDbKegg[input$organismDb],
                            minGSSize = input$minGSSize, 
                            maxGSSize = input$maxGSSize, 
                            pvalueCutoff = input$pvalCuttoff,
                            pAdjustMethod = input$pAdjustMethod,
                            keyType = "ncbi-geneid")
        
        myValues$organismKegg = organismsDbKegg[input$organismDb]
        
        updateSelectInput(session, "geneid_type", choices = gene.idtype.list, selected = input$keytype)
        updateSelectizeInput(session,'pathwayIds', choices=kegg_gse@result$ID)
        
      }, error = function(e) {
        myValues$status = paste("Error: ",e$message)
        
        showNotification(id="errorNotify", myValues$status, type = "error", duration = NULL)
        showNotification(id="errorNotify1", "Make sure the right organism was selected", type = "error", duration = NULL)
        showNotification(id="errorNotify2", "Make sure the corresponding required columns are correctly selected", type = "error", duration = NULL)
        return(NULL)
      }
      
      ), 
      "Error merging files. Check!"))
      
      
    })
    
    shinyjs::show(selector = "a[data-value=\"pubmedTab\"]")
    shinyjs::show(selector = "a[data-value=\"wordcloudTab\"]")
    shinyjs::show(selector = "a[data-value=\"pathviewTab\"]")
    shinyjs::show(selector = "a[data-value=\"keggPlotsTab\"]")
    shinyjs::show(selector = "a[data-value=\"goplotsTab\"]")
    shinyjs::show(selector = "a[data-value=\"gseKeggTab\"]")
    shinyjs::show(selector = "a[data-value=\"gseGoTab\"]")
    
    return(list('go_gse'=go_gse, 'kegg_gse' = kegg_gse))
    
  })
})



output$gseGoTable <- renderDataTable({
  gseGo <- gseGoReactive()
  
  if(!is.null(gseGo)){
    resultDF = gseGo$go_gse@result
    
    DT::datatable(resultDF, options = list(scrollX = TRUE, columnDefs = list(list(visible=input$showAllColumns, targets= 10:12 )) ))
  }
  
})

output$downloadgseGoCSV <- downloadHandler(
  filename = function()  {paste0("gsego",".csv")},
  content = function(file) {
    write.csv(gseGoReactive()$go_gse@result, file, row.names=TRUE)}
)

output$gseGoAvailable <-
  reactive({
    return(!is.null(gseGoReactive()$go_gse))
  })
outputOptions(output, 'gseGoAvailable', suspendWhenHidden=FALSE)


output$gseKEGGTable <- renderDataTable({
  gseKEGG <- gseGoReactive()
  
  if(!is.null(gseKEGG)){
    resultDF = gseKEGG$kegg_gse@result
    
    DT::datatable(resultDF, options = list(scrollX = TRUE, columnDefs = list(list(visible=input$showAllColumns_kegg, targets= 10:11 ))))
  }
  
},
options = list(scrollX = TRUE))

output$downloadgseKEGGCSV <- downloadHandler(
  filename = function()  {paste0("gseKEGG",".csv")},
  content = function(file) {
    write.csv(gseGoReactive()$kegg_gse@result, file, row.names=TRUE)}
)

output$gseKEGGAvailable <-
  reactive({
    return(!is.null(gseGoReactive()$kegg_gse))
  })
outputOptions(output, 'gseKEGGAvailable', suspendWhenHidden=FALSE)


output$warningText <- renderText({
  
  outputText = myValues$convWarningMessage
  if(length(outputText) == 3)
    outputText[3] = paste0('<strong>',outputText[3],'</strong>')
  
  paste("<p>",outputText,"</p>")
})

observeEvent(input$gotoGoPlots, {
  GotoTab('goplotsTab')
})

observeEvent(input$gotoKeggPlots, {
  GotoTab('keggPlotsTab')
})

observeEvent(input$gotoPathview, {
  GotoTab('pathviewTab')
})

observeEvent(input$gotoPubmed, {
  GotoTab('pubmedTab')
})