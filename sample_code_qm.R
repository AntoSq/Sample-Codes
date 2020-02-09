### Portion of the Quantile Matching R code developed during the PhD programme 
### of Antonello A. Squintu, in cooperation between KNMI and Wageningen University
### Details in https://doi.org/10.1002/joc.5874
### This module takes two portion of a series, before and after a break
### (i.e. the cause of an inhomogeneity) and calculates their percentiles month
### by month, so that the differences in the temperature pdf 
### are identified. The same procedure is performed to some homogeneous series
### so that the climatic signal is detected.

library(gdata)

### DIRECTORY

setwd('~/Desktop')

### FILES AND PARAMETERS

file.ori=list.files(pattern=glob2rx("T*ori.txt"))       # File with the inhomogeneous series
files.ref <- list.files(pattern=glob2rx("T*ref.txt"))   # Files with the reference series

chp_y=c(1967, 1989)                                     # Year and month of the inhomogeneities
chp_m=c(5,1)                                            # generally stored in a .RData file
endmonth=c(31,28,31,30,31,30,31,31,30,31,30,31)
thre_y=5              #minimum years of overlapping data needed to calculate the percentiles
ampl_bin=5            #amplitude bins of percentiles
fq=ampl_bin
lq=100-fq
fb=(ampl_bin*3/2)          #first boundary
lb=100-fb                       #last boundary
nb=length(seq(fb,lb,ampl_bin))+2  #number of boundaries                    
nq=nb-1 #number of intervals
if (nq%%2==0)
{
  fb=fb-ampl_bin/2
  lb=lb+ampl_bin/2  
}


### FUNCTIONS (usually in a separate file, pasted here for simplicity)

#source(paste0(sourcepath,'/hom_functions.R'))

### application of the adjustment to the series
adjust <- function(dat,bndrs_ds,nb,corrs)
{
  m=as.numeric(dat['month'])
  value=as.numeric(dat['donat'])
  perc=c(NA,NA)
  add=c(NA,NA)
  exst=NA
  check0=sapply(bndrs_ds, function(x)  if (x[10,m]<(-99)) {return(0)} else {return(1)}) ### check if there are missing data in the boundaries in the months
  check0=sum(check0) ### how many months have non-missing data?
  if (check0>=3) ### at least 3 references?
  {
    for (j in 1: length(refers))
    {
      vec_b=bndrs_ds[[j]][,m]
      if (!is.na(vec_b[2]) && vec_b[2]>-90)
      {
        if(value<vec_b[2]) ### inspect in which quantile bin is the value that has to be adjusted
        {
          perc[j]=5
        }else if (value>vec_b[nb-1])
        {
          perc[j]=95
        }else
        {
          low_b=max(as.numeric(labels(which(vec_b<=value))))
          high_b=min(as.numeric(labels(which(vec_b>=value))))
          if (low_b < high_b & high_b-low_b==5) # check that two consecutive boundaries have been selected (e.g. 17.5 and 22.5)
          {
            perc[j]=(low_b+high_b)/2 ### determine the name of the quantile bin ((17.5,25.5) -> 20)
          } else if(low_b==high_b) ### this can be the case when a quantile sequence is "flat"
          {
            perc[j]=(abs(low_b-50)-2.5)*sign(low_b-50)+50 ### if the boundaries is the same, the percentile is chosen to be the one closer to the median (conservative approach)
            ### (eg. (22.5,22.5) -> 25, (62.5,62.5) -> 60 )
          } else  ### if the two boundaries are too far from eachother (eg. 32.5 and 42.5)
          {
            if(low_b<50 & high_b<50) ### if both are below the median take the one closest to the median  (e.g. (32.5,42.5) -> 45)
            {
              perc[j]=high_b-2.5
            } else if(low_b>50 & high_b>50) ### if both are above the median take the one closest to the median  (e.g. (82.5,92.5) -> 85)
            {
              perc[j]=low_b+2.5
            } else ### if they range over the median, assignt eh value to the centre of the distribution
            {
              perc[j]=50
            }
          }
        }
        
        add[j]=corrs[[j]][which(rownames(corrs[[j]])==perc[j]),m] ### what is the adjustment corresponding the quantile and the month?
        exst[j]=value+add[j] ### calculate the exstimation
      }
    }
  }
  if (length(unlist(exst)[!is.na(unlist(exst))])<3) ### if there are less than 3 exstimations, change into missing value 
  {
    hom=-999.9
  } else ### otherwise take a "wise median" of the values
  {
    hom=wise_median(unlist(exst))
    if(value>hom)
    {
      hom=round(hom+0.01,digits=1) ### in case the median gives a 0.05 resolution and the value (e.g. 16.4) was larger than the hom, this rounds it up (e.g. 17.05 -> 17.1) (conservative approach)
    }else if(value<hom)
    {
      hom=round(hom-0.01,digits=1) ### in case the median gives a 0.05 resolution and the value (e.g. 17.7) was lower than the hom, this rounds it down (e.g. 17.05 -> 17) (conservative approach)
    }
  }
  return(hom)
}

# calculation of the boundaries of the quantile bins for an array of values. 
# Example: for quantile 25 calculates values of quantile 22.5 and 27.5 
bndrs_calc <- function(vec,nq,ampl_bin,thre_y)
{
  if (length(vec[!is.na(vec)])>=thre_y*90)
  {
    ll=length(vec[!is.na(vec)])/(100/ampl_bin)                    #number of data in each interval
    #to calculate values of the bndrs you need to know the center of each interval
    card_qnt1=round(seq(0,nq+1)*ll)
    vec1_table=data.frame(value=sort(vec[!is.na(vec)]),cat=NA)
    for (b in (1:(nq+1)))
    {
      vec1_table[((card_qnt1[b]+1):card_qnt1[b+1]),'cat']=b
    }
    val_bnd1=tapply(vec1_table$value,vec1_table$cat,wise_median)
    val_bnd1[1]=-999.9
    val_bnd1[nq+1]=999.9
  }else
  {
    val_bnd1=rep(-999.9,nq+1)
  }
  val_bnd1=round(val_bnd1,digits = 1)
  return(val_bnd1)
}

# calculation of boundaries: takes as input a time table and calls bndrs_calc.
bndrs_d_3m <- function(overlap,nq,fb,lb,ampl_bin,thre_y1)
{
  bndrs=list()
  if (colnames(overlap)[4]=='donat')
  {
    aux_qnt=list()
    for (m in 1:12)
    {
      if (m==1)
      {
        vec=overlap[overlap$month==12 | overlap$month==1 | overlap$month==2, 'donat']
      }else if (m==12)
      {
        vec=overlap[overlap$month==11 | overlap$month==12 | overlap$month==1, 'donat']
      }else
      {
        vec=overlap[overlap$month==m-1 | overlap$month==m | overlap$month==m+1, 'donat']
      }
      aux_qnt[[m]]=bndrs_calc(vec,nq,ampl_bin,thre_y)
      for (b in 1:(nq+1))
      {
        if(b<=10)
        {
          aux_qnt[[m]][b]=round(aux_qnt[[m]][b]-0.01,digits=1)
        } else
        {
          aux_qnt[[m]][b]=round(aux_qnt[[m]][b]+0.01,digits=1)
        }
      }
    }
    bndrs=sapply(aux_qnt, unlist)
    colnames(bndrs)=c('01','02','03','04','05','06','07','08','09','10','11','12')
    row.names(bndrs)=c(0,seq(fb,lb,ampl_bin),100)
  }
  return(bndrs)
}

### check of the negative slope in a quantile sequence, see Squintu et al.,2019 for the theory
check_neg_slope <- function(corr,qntl)
{
  for (k in 11:19)
  {
    if(corr[k]-corr[k-1]<qntl[k-1]-qntl[k])
    {
      corr[k]=corr[k-1]+qntl[k-1]-qntl[k]
    }
  }
  for (k in 9:1)
  {
    if(corr[k]-corr[k+1]>qntl[k+1]-qntl[k])
    {
      corr[k]=corr[k+1]+qntl[k+1]-qntl[k]
    }
  }
  return(corr)
}

# calculates the matrix of quantile sequences (month x quantile) for basis/donat and for reference, then takes the difference
corr_qntl_ref_3m <- function(overlap,nb,nq,fq,lq,ampl_bin,thre_y1)
{
  corr=list()
  if (colnames(overlap)[4]=='basis')
  {
    aux_qnt=list()
    for (m in 1:12)
    {
      if (m==1)
      {vec=overlap[overlap$month %in% c(1,2,12) & !is.na(overlap$basis),'basis']
      }else if (m==12)
      {vec=overlap[overlap$month %in% c(11,12,1) & !is.na(overlap$basis),'basis']
      }else
      {vec=overlap[overlap$month %in% seq(m-1,m+1) & !is.na(overlap$basis),'basis']}
      aux_qnt[[m]]=qntls_3m(vec,nb,nq,ampl_bin,thre_y)
    }
  }else if (colnames(overlap)[4]=='donat')
  {
    aux_qnt=list()
    for (m in 1:12)
    {
      if (m==1)
      {vec=overlap[overlap$month %in% c(1,2,12) & !is.na(overlap$donat),'donat']
      }else if (m==12)
      {vec=overlap[overlap$month %in% c(11,12,1) & !is.na(overlap$donat),'donat']
      }else
      {vec=overlap[overlap$month %in% seq(m-1,m+1) & !is.na(overlap$donat),'donat']}
      aux_qnt[[m]]=qntls_3m(vec,nb,nq,ampl_bin,thre_y=5)
    }
  }
  aux_r_qnt=list()
  for (m in 1:12)
  {
    if (m==1)
    {
      vec=overlap[overlap$month==12 | overlap$month==1 | overlap$month==2, 'refer']
    }else if (m==12)
    {
      vec=overlap[overlap$month==11 | overlap$month==12 | overlap$month==1, 'refer']
    }else
    {
      vec=overlap[overlap$month==m-1 | overlap$month==m | overlap$month==m+1, 'refer']
    }
    aux_r_qnt[[m]]=qntls_3m(vec,nb,nq,ampl_bin,thre_y)
  }
  qnt=sapply(aux_qnt, unlist)
  colnames(qnt)=c('01','02','03','04','05','06','07','08','09','10','11','12')
  row.names(qnt)=seq(fq,lq,ampl_bin)
  r_qnt=sapply(aux_r_qnt, unlist)
  colnames(r_qnt)=c('01','02','03','04','05','06','07','08','09','10','11','12')
  row.names(r_qnt)=seq(fq,lq,ampl_bin)
  corr=qnt-r_qnt
  return(corr)
}

#calculates overlap between basis/donat and refer
overlapize_ref <- function(ref,ser,thre_y=5) 
{
  overlap=merge(ser,ref)
  overlap=overlap[order(overlap$day),]
  overlap=overlap[order(overlap$month),]
  overlap=overlap[order(overlap$year),]
  if (colnames(overlap)[4]=='basis')
  {
    append=data.frame(year=rep(9999,12),month=seq(1,12,1),day=rep(45,12),basis=NA,refer=NA)  ### makes the process robust
  } else if (colnames(overlap)[4]=='donat')
  {
    append=data.frame(year=rep(9999,12),month=seq(1,12,1),day=rep(45,12),donat=NA,refer=NA)
  }
  overlap=rbind(overlap,append)
  overlap[overlap<=(-50)]<-NA
  if (nrow(overlap[!is.na(overlap[,4])&!is.na(overlap[,5]),])<=thre_y*365)
  {overlap[,4:5]=NA}
  return(overlap)
}

#takes the bndrs of the intervals as input to calculate (with the median) the values of the qntls
qntls_3m <- function(vec,nb,nq,ampl_bin,thre_y=5) 
{
  if (length(vec[!is.na(vec)])>=thre_y*80)
  {
    ll=length(vec[!is.na(vec)])/(100/ampl_bin)                    #number of data in each interval
    #to calculate the values of the quantile you need to know which are the values of the boundaries of each interval
    #card_bnd1 is the array of the indices of the bndrs needed, in this case 2.5, 7.5, 12.5, ... ,97.5 percentiles
    card_bnd=round(seq(0.5,nb-0.5)*ll)
    vec1_table=data.frame(value=sort(vec[!is.na(vec)]),cat=NA)
    for (q in (1:nq)) #for each quantile
    {
      vec1_table[((card_bnd[q]+1):card_bnd[q+1]),'cat']=q #assign every data to an interval, bounded by the boundaries in ca
    }
    val_qnt1=tapply(vec1_table$value,vec1_table$cat,wise_median)
  }else
  {
    val_qnt1=rep(NA,nq)
  }
  val_qnt1=round(val_qnt1,digits = 1)
  return(val_qnt1)
}

# smooths out matrices with 12 columns 
smoothen_nqx12 <- function(matrice) 
{
  nq=dim(matrice)[1]
  smooth=matrix(NA,nq,12)
  for (m in 1:12)
  {
    if (m==1)
    {
      for (q in 1:nq)
      {
        if (q==1)
        {
          smooth[q,m]=(matrice[1,12]+2*matrice[1,1]+matrice[1,2]+matrice[2,1])/5
        }else if (q==nq)
        {
          smooth[q,m]=(matrice[nq,12]+2*matrice[nq,1]+matrice[nq,2]+matrice[nq-1,1])/5
        }else
        {
          smooth[q,m]=(matrice[q,12]+matrice[q,1]+matrice[q,2]+matrice[q-1,1]+matrice[q+1,1])/5
        }
      }
    }else if (m==12)
    {
      for (q in 1:nq)
      {
        if (q==1)
        {
          smooth[q,m]=(matrice[1,11]+2*matrice[1,12]+matrice[1,1]+matrice[2,12])/5
        }else if (q==nq)
        {
          smooth[q,m]=(matrice[nq,11]+2*matrice[nq,12]+matrice[nq,1]+matrice[nq-1,12])/5
        }else
        {
          smooth[q,m]=(matrice[q,11]+matrice[q,12]+matrice[q,1]+matrice[q-1,12]+matrice[q+1,12])/5
        }
      }
    }else
    {
      for (q in 1:nq)
      {
        if (q==1)
        {
          smooth[q,m]=(matrice[1,m-1]+2*matrice[1,m]+matrice[1,m+1]+matrice[2,m])/5
        }else if (q==nq)
        {
          smooth[q,m]=(matrice[nq,m-1]+2*matrice[nq,m]+matrice[nq,m+1]+matrice[nq-1,m])/5
        }else
        {
          smooth[q,m]=(matrice[q,m-1]+matrice[q,m]+matrice[q,m+1]+matrice[q-1,m]+matrice[q+1,m])/5
        }
      }
    }
  }
  smooth=round(smooth,digits=1)
  return(smooth)
}

#median, but, if the result is not with the set number of digits after the point, the value is rounded up or down in the direction of the mean
wise_median <- function(x,digits=1) 
{
  y=median(x,na.rm=TRUE)   
  if ((y*(10^digits)*10)%%10==5)
  {
    z=mean(x,na.rm=TRUE)
    if(z<y)
    {y=y-5*(10^(-digits-1))}
    else if(z>=y)
    {y=y+5*(10^(-digits-1))}
  }
  return(round(y,digits=digits))
}

### OPERATIONS

original=read.table(file.ori)
colnames(original)=c('year','month','day','value')
refers=lapply(files.ref,read.table)
names(refers)=files.ref
for (j in 1:length(refers)) {colnames(refers[[j]])=c('year','month','day','refer')}
donats=list()                                           # list with the inhomogeneous parts of the original
series_name=file.ori

ser_rest=original
i=0
for (k in 1:length(chp_y)) ### for each break, the part preceeding it is stored in the donating subseries
{
  i=i+1
  ser_fin=ser_rest[ser_rest$year<chp_y[k] | (ser_rest$year==chp_y[k] & ser_rest$month<chp_m[k]),]
  ser_rest=ser_rest[ser_rest$year>chp_y[k] | (ser_rest$year==chp_y[k] & ser_rest$month>=chp_m[k]),]
  donats[[i]]=ser_fin
  colnames(donats[[i]])=c('year','month','day','donat')
  names(donats)[i]=paste0(series_name,'-0',k)
}
basis=ser_rest ### the remaining part, after the last break, is the basis
colnames(basis)[4]='basis'

overlap_dr=list()
bndrs_ds=list()
adjs=list()
adjs02=list()
qntl_drs=list()

for (d in length(donats):1)
{
overlap_br=lapply(refers, function(x) overlapize_ref(x,basis));   ### overlap between basis and reference
qntl_brs=lapply(overlap_br, function(x) corr_qntl_ref_3m(x,nb,nq,fq,lq,ampl_bin,thre_y)) ### quantiles of basis and references

donat=donats[[d]]
overlap_dr[[d]]=list()
qntl_drs[[d]]=list()
bndrs_ds[[d]]=list()
adjs[[d]]=list()
adjs02[[d]]=list()
donat_name=names(donats)[d]
print(paste('Adjustment calculation of',names(donats[[1]])))

overlap=merge(basis,donat)  ### overlap of basis and donat
error=FALSE
if (dim(overlap)[1]!=0)
{
  error=1=TRUE
  print('There is an overlap! Error somewhere')
}
if (error==FALSE)
{
  {
    overlap_dr[[d]]=lapply(refers, function(x) overlapize_ref(x,donat));                ### overlap between donors and references
    qntl_drs[[d]]=lapply(overlap_dr[[d]], function(x) corr_qntl_ref_3m(x,nb,nq,fq,lq,ampl_bin,thre_y));  ### quantiles don vs ref
    bndrs_ds[[d]]=lapply(overlap_dr[[d]], function(x) bndrs_d_3m(x,nq,fb,lb,ampl_bin,thre_y)) ### boundaries of the quantiles
    adjs_rough=list()  
    adjs0=matrix(NA,nq,12)
    adjs02[[d]]=matrix(NA,nq,12)     
    adjs_cns=list()
    adjs_sm=list()

    for (j in 1:length(refers)) ### loop on the references
    {
      
      adjs_rough[[j]]=round(qntl_brs[[j]]-qntl_drs[[d]][[j]],digits=1)    ### adjustments without smoothing              
      if(length(adjs_rough[[j]][!is.na(adjs_rough[[j]])])>=nq)
      {
        adjs_sm[[j]]=round(smoothen_nqx12(adjs_rough[[j]]),digits=1) ### smooth the adjustments
        overlap=overlap_br[[j]] 
        adjs[[d]][[j]]=matrix(NA,19,12)
        for (m in 1:12) ### loop on the months to perform the check of the negative slopes
        {
          adj_sm=adjs_sm[[j]][,m]
          if (length(adj_sm[!is.na(adj_sm)])>0)
          {
            if (m==1)       {vec=overlap[overlap$month %in% c(1,2,12) & !is.na(overlap$basis),'basis']
            }else if (m==12){vec=overlap[overlap$month %in% c(11,12,1) & !is.na(overlap$basis),'basis']
            }else           {vec=overlap[overlap$month %in% seq(m-1,m+1) & !is.na(overlap$basis),'basis']}
            qntl=qntls_3m(vec,nb,nq,ampl_bin,thre_y) ### quantiles of this specific case are needed to know what the max neg slope can be
            adjs[[d]][[j]][,m]=check_neg_slope(adj_sm,qntl) ### checks if negative slopes are less than the slopes in quantiles (and smooths them)    
          }else
          {
            adjs[[d]][[j]][,m]=rep(NA,19)
          }
        }                      
        colnames(adjs[[d]][[j]])=seq(1,12)
        rownames(adjs[[d]][[j]])=seq(fq,lq,ampl_bin)
      }
    }
  }
}

### ADJUSTMENT OF THE SERIES
homdon=donat[c('year','month','day','donat')] ### initializing the homogenized version of the donating series    
homdon$hom=apply(donat,1,function(x) adjust(x,bndrs_ds[[d]],nb,adjs[[d]]))                       
colnames(homdon)[5]='basis' ### change of column name to allow the concatenation
basis=rbind(homdon[,c('year','month','day','basis')],basis)  ### the adjusted subseries is preappended to the basis
### the new longer basis is used for the backwards homogenization of the previous donating portion
}

homogenized=basis ### final homogenized series
write.fwf(homogenized,'hom_series.txt')
