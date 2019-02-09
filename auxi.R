
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




#
