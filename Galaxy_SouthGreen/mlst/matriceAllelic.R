matrix.allelic = function(comp, rac, name) 
{
  fm <- paste(rac,"-matrix.csv",sep="");
  matrix.1 <- read.table(fm, header=TRUE, row.names=1 , sep=";");
  pa <- vector("numeric", length = 0);
  for(j in 1:8) { 
    seq = read.dna(paste(rac,j,'.fas',sep="") ,format = "fasta", as.character = FALSE, as.matrix = FALSE);
    locus(comp,seq)
    if (pa[])
  }      
  msg <- name;
}
  
locus = function(comp,seq)  
{
  
    trouve=0;
    len = length(seq);
    for(i in 1:len) {
      st = as.character(seq[i]);
      if(as.alignment(st)$seq == tolower(comp))  trouve=i;
    }
    if(trouve >0) {
      pa[j] = trouve;
    } else { 
      pa[j] = max(matrix.1[,j])+1;
    }
    msg=paste(msg,pa[j], sep="\t");
  }
  
  return(pa);
}