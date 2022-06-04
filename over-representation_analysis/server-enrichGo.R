
myValues = reactiveValues()

observe({
  enrichGoReactive()
})

enrichGoReactive <- eventReactive(input$initGo,{
  
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
        
        # omitir valores NA 
        gene_list<-na.omit(original_gene_list)
        
        # ordenar la lista en orden decreciente (requerido para clusterProfiler
        gene_list = sort(gene_list, decreasing = TRUE)
        
        myValues$gene_list = gene_list
        
        
        sig_genes_df = df[df[,input$padjColumn] < input$padjCutoff,]
        sig_genes_df = na.omit(sig_genes_df)
        
        # A partir de resultados significativos, filtrar el cambio log2fold
        genes <- sig_genes_df[[input$log2fcColumn]]
        
        # nombre del vector
        names(genes) <- sig_genes_df[[input$geneColumn]]
        
        # omitir valores NA
        genes <- na.omit(genes)
        
        # filtrar en cambio mínimo log2fold (PARÁMETRO)
        genes <- names(genes)[abs(genes) > input$logfcCuttoff]
        
        
        setProgress(value = 0.3, detail = "Performing Go enrichment analysis, please wait ...")
        orgDb.obj = eval( parse(text = input$organismDb, keep.source=FALSE))
        
        go_enrich <- enrichGO(gene = genes,
                              universe = names(gene_list),
                              OrgDb = orgDb.obj, 
                              keyType = "ENSEMBL",
                              minGSSize = input$minGSSize, 
                              maxGSSize = input$maxGSSize,
                              readable = T,
                              ont = input$ontology,
                              pvalueCutoff = input$pvalCuttoff, 
                              qvalueCutoff = input$qvalCuttoff,
                              pAdjustMethod = input$pAdjustMethod)
        
        
        if(nrow(go_enrich) < 1)
        {
          showNotification(id="warnNotify", "No gene can be mapped ...", type = "warning", duration = NULL)
          showNotification(id="warnNotify2", "Tune the parameters and try again.", type = "warning", duration = NULL)
          return(NULL)
        }
        
        updateNumericInput(session, "showCategory_bar", max = nrow(go_enrich@result) , min = 0, value = ifelse(nrow(go_enrich@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_dot", max = nrow(go_enrich@result) , min = 0, value = ifelse(nrow(go_enrich@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_enrichmap", max = nrow(go_enrich@result) , min = 0, value = ifelse(nrow(go_enrich@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_goplot", max = nrow(go_enrich@result) , min = 0, value = ifelse(nrow(go_enrich@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_cnet", max = nrow(go_enrich@result) , min = 0, value = ifelse(nrow(go_enrich@result) > 0, 5,0))
        
        updateNumericInput(session, "showCategory_bar_kegg", max = nrow(go_enrich@result) , min = 0, value = ifelse(nrow(go_enrich@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_dot_kegg", max = nrow(go_enrich@result) , min = 0, value = ifelse(nrow(go_enrich@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_enrichmap_kegg", max = nrow(go_enrich@result) , min = 0, value = ifelse(nrow(go_enrich@result) > 0, 5,0))
        #updateNumericInput(session, "showCategory_goplot", max = nrow(go_enrich@result) , min = 0, value = ifelse(nrow(go_enrich@result) > 0, 5,0))
        updateNumericInput(session, "showCategory_cnet_kegg", max = nrow(go_enrich@result) , min = 0, value = ifelse(nrow(go_enrich@result) > 0, 5,0))
        
        
        ## KEGG enrich
        
        # Conviertir ID de genes para la función enrichKEGG
        # Se pierden algunos genes aquí porque no todas las identificaciones se convertirán
        myValues$convWarningMessage = capture.output(ids<-bitr(names(original_gene_list), fromType = input$keytype, toType = "ENTREZID", OrgDb=input$organismDb), type = "message")
        
        # elimine IDS duplicados (aquí uso "ENSEMBL", pero debería ser lo que se seleccionó como keyType)
        dedup_ids = ids[!duplicated(ids[c(input$keytype)]),]
        
        # Cree un nuevo marco de datos df2 que tenga solo los genes que se mapearon con éxito usando la función bitr anterior
        df2 = df[df[[input$geneColumn]] %in% dedup_ids[,input$keytype],]
        
        # Crear una nueva columna en df2 con los ID de ENTREZ correspondientes
        df2$Y = dedup_ids$ENTREZID
        
        # Crear un vector del universo de genes.
        kegg_gene_list <- df2[[input$log2fcColumn]]
        
        # Vector de nombre con ID de ENTREZ
        names(kegg_gene_list) <- df2$Y
        
        # omitir valores NA 
        kegg_gene_list<-na.omit(kegg_gene_list)
        
        # ordenar la lista en orden decreciente (requerido para clusterProfiler)
        kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
        
        myValues$kegg_gene_list = kegg_gene_list
        
        
        kegg_sig_genes_df = df2[df2[,input$padjColumn] < input$padjCutoff,]
        kegg_sig_genes_df = na.omit(kegg_sig_genes_df)
        
        # A partir de resultados significativos, filtrar el cambio log2fold
        kegg_genes <- kegg_sig_genes_df[[input$log2fcColumn]]
        
        # ¡Nombre el vector con la ID CONVERTIDA!
        names(kegg_genes) <- kegg_sig_genes_df$Y
        
        # omitir valores NA
        kegg_genes <- na.omit(kegg_genes)
        
        # filtro en cambio log2fold (PARÁMETRO)
        kegg_genes <- names(kegg_genes)[abs(kegg_genes) > input$logfcCuttoff]
        
        
        setProgress(value = 0.6, detail = "Performing KEGG enrichment analysis, please wait ...")
        
        organismsDbKegg = c("org.Hs.eg.db"="hsa","org.Mm.eg.db"="mmu","org.Rn.eg.db"="rno",
                            "org.Sc.sgd.db"="sce","org.Dm.eg.db"="dme","org.At.tair.db"="ath",
                            "org.Dr.eg.db"="dre","org.Bt.eg.db"="bta","org.Ce.eg.db"="cel",
                            "org.Gg.eg.db"="gga","org.Cf.eg.db"="cfa","org.Ss.eg.db"="ssc",
                            "org.Mmu.eg.db"="mcc","org.EcK12.eg.db"="eck","org.Xl.eg.db"="xla",
                            "org.Pt.eg.db"="ptr","org.Ag.eg.db"="aga","org.Pf.plasmo.db"="pfa",
                            "org.EcSakai.eg.db"="ecs")
        
        kegg_enrich <- enrichKEGG(gene=kegg_genes, 
                                  universe=names(kegg_gene_list),
                                  organism=organismsDbKegg[input$organismDb], 
                                  pvalueCutoff =input$pvalCuttoff,
                                  qvalueCutoff = input$qvalCuttoff,
                                  keyType = "ncbi-geneid",
                                  minGSSize = input$minGSSize, 
                                  maxGSSize = input$maxGSSize)
        
        myValues$organismKegg = organismsDbKegg[input$organismDb]
        
        updateSelectInput(session, "geneid_type", choices = gene.idtype.list, selected = input$keytype)
        updateSelectizeInput(session,'pathwayIds', choices=kegg_enrich@result$ID)
        
      }, error = function(e) {
        
        myValues$status = paste("Error: ",e$message)
        
        showNotification(id="errorNotify", myValues$status, type = "error", duration = NULL)
        showNotification(id="errorNotify1", "Make sure the right organism was selected", type = "error", duration = NULL)
        showNotification(id="errorNotify2", "Make sure the corresponding required columns are correctly selected", type = "error", duration = NULL)
        return(NULL)
      }
      
      ), 
      "Error. Check!"))
      
      
    })
    
    
    
    shinyjs::show(selector = "a[data-value=\"wordcloudTab\"]")
    shinyjs::show(selector = "a[data-value=\"pathviewTab\"]")
    shinyjs::show(selector = "a[data-value=\"keggPlotsTab\"]")
    shinyjs::show(selector = "a[data-value=\"goplotsTab\"]")
    shinyjs::show(selector = "a[data-value=\"enrichKeggTab\"]")
    shinyjs::show(selector = "a[data-value=\"enrichGoTab\"]")
    
    return(list('go_enrich'=go_enrich, 'kegg_enrich' = kegg_enrich))
    
  })
})

output$enrichGoTable <- renderDataTable({
  enrichGo <- enrichGoReactive()
  
  if(!is.null(enrichGo)){
    resultDF = enrichGo$go_enrich@result
    if(isFALSE(input$showGeneidGo))
      resultDF = resultDF[,-which(names(resultDF) == "geneID")]
    
    DT::datatable(resultDF, options = list(scrollX = TRUE))
  }
  
},
options = list(scrollX = TRUE))

output$downloadEnrichGoCSV <- downloadHandler(
  filename = function()  {paste0("enrichgo",".csv")},
  content = function(file) {
    write.csv(enrichGoReactive()$go_enrich@result, file, row.names=TRUE)}
)

output$enrichGoAvailable <-
  reactive({
    return(!is.null(enrichGoReactive()$go_enrich))
  })
outputOptions(output, 'enrichGoAvailable', suspendWhenHidden=FALSE)


output$enrichKEGGTable <- renderDataTable({
  enrichKEGG <- enrichGoReactive()
  
  if(!is.null(enrichKEGG)){
    resultDF = enrichKEGG$kegg_enrich@result
    if(isFALSE(input$showGeneidKegg))
      resultDF = resultDF[,-which(names(resultDF) == "geneID")]
    
    DT::datatable(resultDF, options = list(scrollX = TRUE))
  }
  
},
options = list(scrollX = TRUE))

output$downloadEnrichKEGGCSV <- downloadHandler(
  filename = function()  {paste0("enrichKEGG",".csv")},
  content = function(file) {
    write.csv(enrichGoReactive()$kegg_enrich@result, file, row.names=TRUE)}
)

output$enrichKEGGAvailable <-
  reactive({
    return(!is.null(enrichGoReactive()$kegg_enrich))
  })
outputOptions(output, 'enrichKEGGAvailable', suspendWhenHidden=FALSE)


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

observeEvent(input$gotoWordcloud, {
  GotoTab('wordcloudTab')
})

observeEvent(input$gotoWordcloud1, {
  GotoTab('wordcloudTab')
})