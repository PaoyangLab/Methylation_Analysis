# Loading packages
foo <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! suppressPackageStartupMessages(require( i , character.only = TRUE )) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      suppressPackageStartupMessages(require( i , character.only = TRUE , quietly = TRUE))
    }
  }
}
foo( c("argparser" , "tidyverse") )

p <- arg_parser("Conversion_rate")
p <- add_argument(p, "input", help="input CGmap", type="character")
argv <- parse_args(p)

cat("[",format(Sys.time(), "%X"),"]","Calculating bisulfite conversion rate","\n")

# Import CGmap file
CGmap <- read.table(argv$input, header = F)
colnames(CGmap) <- c("chr","nt","site","type","text","C/(C+T)","C","C+T")
# Calculate bisulfite conversion rate
conversion_rate <- ((sum(CGmap$C+T-CGmap$C)) / sum(CGmap$C+T)) * 100

# Print conversion rate
cat("[",format(Sys.time(), "%X"),"]","Bisulfite conversion rate:", conversion_rate, "%","\n")