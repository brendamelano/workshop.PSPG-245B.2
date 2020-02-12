
k2 <- function (df, tl=""){
    kable( df , format = "html", booktabs = T, caption = tl) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"))
}

samples = c("TCGA-C8-A3M7-01","TCGA-GM-A2D9-01","TCGA-BH-A1FL-01", "TCGA-E2-A14Z-01","TCGA-BH-A1FC-01","TCGA-AC-A2QJ-01", "TCGA-EW-A1P8-01","TCGA-E2-A1LE-01","TCGA-E2-A1LK-01","TCGA-AC-A2FE-01")

url = 'https://search.vumc.org/?query='
raw = 0 
if ( raw == 1){
    cosmic.70 = "./external/hg19_cosmic70.txt"
    
    cosmic.70 = read.table(cosmic.70, header=F,sep="\t",stringsAsFactors = FALSE,na.strings=".",  quote = "")
    
    colnames ( cosmic.70) = c("chr","start","end","ref","alt","description" )
    
    
    cosmic.70$igv = paste ( cosmic.70$chr, cosmic.70$start, cosmic.70$end, cosmic.70$ref, cosmic.70$alt)
    
    
    cancer.list = read.csv("https://www.dropbox.com/s/naheek0wicegf77/cancer.list.csv?dl=1")
    
    
    drug  = read.csv( "./external/civic.txt" )
    
    
    library("openxlsx")
    
    drug <- read.xlsx("./external/civic.xlsx" , sheet="civic", colNames = TRUE)
    drug$match = paste ( drug$gene, drug$variant)
    
    saveRDS( list (cosmic.70=cosmic.70, cancer.list = cancer.list, drug = drug, 
                   studies=studies, mycaselist=mycaselist , mygeneticprofile=mygeneticprofile, mutations = mutations
                   
                   ), "workshop.rds")
    
}else{

temp = readRDS("workshop.rds")
cosmic.70 = temp$cosmic.70
cancer.list = temp$cancer.list
drug = temp$drug
cancer.gene = as.character ( cancer.list$gene )


}




venn.this <- function (data1, cp = c("#a6cee3","#fdbf6f","#b2df8a"), type= 3, dgg=0, title1 = "", title2="", cexL = 3, cexC=3, titlesize=25, titlecol= "steelblue" ) {
    
    if ( type==3){
        cat = names (data1 )
        a1 = data.table(unlist ( data1[1] )  )
        a2 = data.table(unlist ( data1[2] ) )
        a3 = data.table(unlist ( data1[3] ) )
        
        n12 = fintersect(a1,a2)
        n23 = fintersect(a2,a3)
        n13 = fintersect(a1,a3)
        n123 = fintersect(n12,a3)
        grid.newpage()
        v3= draw.triple.venn( 
            length ( a1$V1), 
            length ( a2$V1), 
            length ( a3$V1), 
            length (n12$V1), 
            length ( n23$V1), 
            length ( n13$V1), 
            length ( n123$V1), 
            category = c( cat[1], cat[2], cat[3] ), 
            fill = cp, 
            rotation.degree = dgg, 
            euler           = F,
            scaled          = FALSE
            , cex = cexL
            , cat.cex = cexC
            #,cat.pos         = c(0, 0, 0)
        )
        
        #v3 = grid.arrange(gTree(children=v3), top=textGrob(title1, gp=gpar(fontsize=titlesize, col=titlecol) )
        #                         , bottom=textGrob(title2, gp=gpar(fontsize=titlesize, col=titlecol) ))
        
        # get all intersection and union 
        
        n.all = n123$V1
        # unique to a1 
        u.a1 = setdiff(a1$V1, unique ( c(a2$V1, a3$V1) ) )
        u.a2 = setdiff(a2$V1, unique ( c(a1$V1, a3$V1) ) )
        u.a3 = setdiff(a3$V1, unique ( c(a1$V1, a2$V1) ) )
        
        v3 =  as.ggplot( grobTree(v3) ) + ggtitle ( title1 )
        
        main.cmp =  qpcR:::cbind.na ( sort ( u.a1, decreasing=T) 
                                      , sort ( u.a2, decreasing=T) 
                                      , sort ( u.a3, decreasing=T) 
                                      , sort ( n.all, decreasing=T) 
                                      
        )
        main.cmp = data.frame ( main.cmp , stringsAsFactors = F)
        main.cmp[is.na(main.cmp)] <-  '.' 
        
        
        colnames ( main.cmp) = c( cat[1], cat[2], cat[3], "all.three" )
        
        
        return ( list ( venn = v3, main.cmp=main.cmp ) )
        
    } else if ( type == 2){
        
        grid.newpage()
        names.cat = names ( data1)
        ## data1 is one row with n1, n2 and int as colnames
        v3 = draw.pairwise.venn(             area1           = data1[,1],
                                             area2           = data1[,2],
                                             cross.area      = data1[,3],
                                             category        = c(names.cat[1], names.cat[2]),
                                             fill            = cp[1:2],
                                             lty             = "blank",
                                             cex             = 2,
                                             cat.cex         = 2,
                                             cat.pos         = c(180, 160),
                                             cat.dist        = - .03,
                                             #cat.just        = list(c(-1, -1), c(1, 1)),
                                             ext.pos         = 30,
                                             ext.dist        = -0.05,
                                             ext.length      = 0.85,
                                             ext.line.lwd    = 2,
                                             ext.line.lty    = "dashed", 
                                             rotation.degree = dgg # this is to flip left to right cat
        )
        
     
        
        
        
    }
    
    
}



frozen = readRDS("freeze.2020")

#
