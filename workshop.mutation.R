library ( cgdsr)
library ( ggplot2)
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")


# Get list of cancer studies at server
studies = getCancerStudies(mycgds)
# lets take a look at it. 
View ( head ( studies, 10) )
# what columns are in here, this gives you an idea what the dataframe consist of.  
names ( studies )
# lets find how many unique sets are in here 
dim ( studies )
# so it looks like 240 studies 

# Find the dataset you want. 
# According to the syllabus we are looking for TCGA study on breast cancer from Cell, 2015 
# so lets search for it using the grepl function 
breast = studies[ grepl("Breast", studies$name, ignore.case = T ), ]
View ( breast )

### BONUS EXERCISE: Level 1:  TRY TO FIND SOME OTHER CANCER GROUP, eg Lung Cancer

# from there we find what we are looking for and its called, brca_tcga_pub2015

brca.study = "brca_tcga_pub2015"

# lets see if you sample is in here. 
mycaselist = getCaseLists(mycgds,brca.study)
# first lets study what is in the case list.  
View ( head ( mycaselist, 10) ) 
# here we see that brca_tcga_pub2015_3way_complete is probably the best to use because it includes
# only All Complete Tumors vs something like All tumor samples with methylation data
case.list.id = "brca_tcga_pub2015_3way_complete"
### BONUS EXERCISE: Level 2: count the total categories available for this set. 

# ok now using All Complete Tumors lets see if your sample is present. 

mysample = mycaselist[ mycaselist$case_list_id == case.list.id, ]$case_ids
mysample[ grepl("TCGA-OL-A66K-01", mysample)]
mysample[ grepl("TCGA-FAKE-A66K-01", mysample)]

### BONUS EXERCISE: Level 1: search for you sample and see if its there. 


### Level 3: BONUS EXERCISE: loop through each one and find the the samples that All tumor samples with methylation data AND either in ER- breast tumors OR Her2-positive breast tumors


# Now lets see what information is available for your study 
mygeneticprofile = getGeneticProfiles(mycgds,brca.study)
View ( mygeneticprofile )
# so looking at that we now know it contains several interesting modalities.  Lets try the mutation
mutation = mygeneticprofile[11, 1]

# ok now we can get the actual mutation data. However first lets download a set of genes that are known cancer genes
# grab this from cosmic 
cosmic = read.csv("https://www.dropbox.com/s/naheek0wicegf77/cancer.list.csv?dl=1")
# we dont need the annotations just the gene however lets take a look at this. 
View ( head (cosmic ))

cancer.gene = as.character ( unique ( cosmic$gene) )
length( cancer.gene)
# as you can see there are 719 genes. 
# lets see if BRCA is in here to make sure! 
cosmic[grepl("^BRCA|^ATM$|^BARD1$|^CDH1$|^CHEK2$|^NBN$|^NF1$|^PALB2$|^PTEN$", cosmic$gene), ]

# BONUS why did we add ^ in front and $ for only some genes and not others?

### BONUS can you check if your favorite gene is in here? 
# lets collect this through a loop so not to overwheml the system 

total = ceiling ( length( cancer.gene)/100  )  *100
mutations = data.frame ( stringsAsFactors = F )
e = 1
for(i in seq(from=200, to=total, by=200)){
   
    if ( i > length(cancer.gene)){
        i = length(cancer.gene)
    }
    
    print ( paste ( e, i ))
    temp = getMutationData(mycgds , brca.study, mutation, cancer.gene[e:i])
    e = i 
    mutations = rbind ( mutations, temp)
    
    Sys.sleep (2) # lets give the system a break
    
}

dim ( mutations )
colnames ( mutations )

# lets take a few minutes here to go over the different fields. 

# lets check if ALL your samples are availble 

samples = c("TCGA-C8-A3M7-01","TCGA-GM-A2D9-01","TCGA-BH-A1FL-01", "TCGA-E2-A14Z-01","TCGA-BH-A1FC-01","TCGA-AC-A2QJ-01", "TCGA-EW-A1P8-01","TCGA-E2-A1LE-01","TCGA-E2-A1LK-01","TCGA-AC-A2FE-01")

final =  mutations [ mutations$case_id %in% samples, ]


unique ( final$case_id)

# lets study your brca! 

# lets see what are the mutations types. 

unique ( mutations$mutation_type)
# lets talk about each one. 

# lets tabulate the types of mutations mostly seen with BRCA 

mutation.type = data.frame ( table ( mutations$mutation_type))

# BONUS compare this with another disease cohort. 



mutation.type = mutation.type[ order ( mutation.type$Freq), ]
mutation.type$Var1 = factor ( mutation.type$Var1, levels = mutation.type$Var1)

ggplot(mutation.type, aes(Var1,Freq), label=Freq ) +
    geom_bar(aes(fill = Freq), stat="identity", position = "dodge") +
    coord_flip() +
    scale_fill_distiller(palette = "RdBu") + xlab("") + ylab("") +
    theme(strip.text.y = element_text(angle = 0), legend.position="none") +
    geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=.4, hjust = .5, size=5) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), 
          text = element_text(size=16),  # size of label
          axis.text.x = element_text(angle=0, hjust=1) )

# lets figure out what type of mutations are BRCA
### how would you do this?  
#### .... 5 mins 

mutation.brca = data.frame ( table ( mutations[ grepl("BRCA", mutations$gene_symbol), ] $mutation_type))

mutation.brca = mutation.brca[ order ( mutation.brca$Freq), ]
mutation.brca$Var1 = factor ( mutation.brca$Var1, levels = mutation.brca$Var1)

ggplot(mutation.brca, aes(Var1,Freq), label=Freq ) +
    geom_bar(aes(fill = Freq), stat="identity", position = "dodge") +
    coord_flip() +
    scale_fill_distiller(palette = "RdBu") + xlab("") + ylab("") +
    theme(strip.text.y = element_text(angle = 0), legend.position="none") +
    geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=.4, hjust = .5, size=5) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), 
          text = element_text(size=16),  # size of label
          axis.text.x = element_text(angle=0, hjust=1) )

# Bonus can you figure out which samples have the BRCA samples and which don't?

# now lets see what are the main genes that are present in missense mutations. 

genes.missense = mutations[mutations$mutation_type == "Missense_Mutation", ]

genes.missense.genes = data.frame ( table ( genes.missense$gene_symbol))
genes.missense.genes = genes.missense.genes[ order ( -genes.missense.genes$Freq), ]

View ( head ( genes.missense.genes, 15))

