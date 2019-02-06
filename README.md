# workshop.mutation.2019
In this workshop you will learn about the fundamentals of mutation identification in cancer.  This will include the following. 

1. Basic bioinformatic workflow: from sample to identifying potentially targetable mutations. 
Types of sequencing, e.g. WES and WGS, and also DNA vs RNA
Standard pipelines to go from reads to actionable information
The size of modern data sets (TCGA) and their current applications, e.g. in healthcare
What does a standard genome look like? For example, 4-5 million variants in an average genome compared to the reference human genome. Then this can lead in to talking about germline vs somatic mutations, and how this is crucial for studying cancer 

2. Basic vocabulary and concepts. 
* Classes of somatic mutations

# Pt mutations
* Coding 
* Silent
* Missense
* Nonsense
* Noncoding ( UTR ) 
* Intronic
* Intergenic
* Splice site variants?

# Small regional mutations
* Insertion
* Deletions
* Duplications

3. Deciphering nomenclature of sequence variations.

What is his section? Is this deciphering single letter vs three letter mutatioion codes etc.?


4. Identifying functionally relevant mutations - passenger vs driver
How RNA and DNA sequencing data can be integrated to find functional variants?
Variant prediction tools – e.g. CADD scores, and other methods
Comparison to known cancer genes or even known cancer causing variants
Pathway analysis – how this leads in to the design/search of drugs
Structural biology for coding variants


Workshop ( this will constitute the bulk of the workshop )

* We will start off with basic data mining.  To begin with will learn how to directly download mutation data from R.  There are many sources and API's however here will be using  cbioportal. 
* Although this is not a course in R per se, but you will learn how to manipulate/wrangle a the mutation data.frame.  
* Subset type of mutations.  
* Aggregate by attributes such as types of mutations. 
* Query ( eg. for specific variants ) 
* Tabulate mutations ( eg. frequency tables ) . 
* How to identify what could be potentially be pathogenic and cross reference it with existing data. 
* How to take existing mutation data and predict possible actionable targets for either druggability, diagnostic or prognosis. 
* You will also learn a few ways to plot the data. 
* Basic plotting of your mutation table. 
* How to generate figures to look for total burden across different chromosome/regions.
* Potentially – how to analyse mutational signatures (if there is time)



