#args<-commandArgs(T) 
#sv_data<-args[1]
library('argparse')
parser <- ArgumentParser(description="Calculation of the posterior pairing probability")
parser$add_argument('--input', type="character", help="input sv data, ")
args <- parser$parse_args()
sv_data<-args$input

source('code/classify/palimpsest_main.R')
sv_data<-read.csv(sv_data, colClasses=c('character','character','character','numeric','character','numeric'))
signature<-palimpsest_main(sv_data) 

