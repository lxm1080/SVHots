library('argparse')

parser <- ArgumentParser()
parser$add_argument('--a', type="double") #number of element in hot
parser$add_argument('--b', type="double") #hotlen
parser$add_argument('--c', type="double") #number of element in nonhot
parser$add_argument('--d', type="double") #nonhotlen
args <- parser$parse_args()

a<-args$a
b<-args$b
c<-args$c
d<-args$d


output<-poisson.test(x=a,T=b,r=c/d,alternative="two.sided")

p<-signif(output$p.value,3)

cat(p)
