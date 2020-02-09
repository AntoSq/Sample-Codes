### SCRIPT FOR THE CALCULATION OF QUANTILE MATCHING ADJUSTMENTS
# Candidate and reference series must be in separate files in the inpath with format 'yyyy mm dd value'.

# Name of the file is fundamental:
# paste0(ele,'_serid_',ser_id,'_ori.txt') for the candidate ('ori.txt' and 'parthom.txt' 
# _ori.txt indicates the original series, while _parthom.txt indicates the partially homogenized series after the 1st it.
# paste0(ele,'_serid_',ser_id,'_ref.txt') for the references ('parthomref.txt' and 'oriref.txt')
# where oriref are references in their original form and parthomref are reference partially homogenized
# where ele is 'TN', 'TX' or 'TG', and ser_id is a positive integer lower than 999999

# breaks can be the output of '../BreakDetection/BreakDetectionSoftware/Break_daily_w_arg.R'
# Otheriwse
# breaks must be written in the format of a list
# names of the elements have to be in the format xTX_serid_000000_ori.txt
# (for all the series, including references the extension has to be '_ori.txt') ('_ori.txt','_parthom.txt' for iter 2)
# the fields of the elementhave to include $range, $all2 and $month
# the list of breks must be stored in the inpath as an R workspace
# with the name 'n_breaks.RData', 'x_breaks.RData' or 'g_breaks.RData'

### READING ARGUMENTS

args <- commandArgs(trailingOnly = TRUE)

if(length(args)!=4){
  stop("Please provide basedir iteration element debugflag")
}else{
  basepath=args[1]
  iter=args[2]
  ele=args[3]
  debugflag=args[4]
}


### PACKAGES

suppressMessages(library(gdata))
suppressMessages(library(trend))

###DIRECTORIES

#basepath='/usr/people/squintu/Documents/ECAD/Package/'
inpath=paste0(basepath,'/WorkFolder')                                 #where are the series?
sourcepath='/home/besselaa/HomScripts/Functions'                        #where are the functions? 
outpath=paste0(basepath,'/Output_temp')                                #where to put the temporary results?
finalpath=outpath                                                   #where to put the final results?

##PARAMETERS
thre_y=5              #minimum years of overlapping data needed to calculate the percentiles
ampl_bin=5            #amplitude bins of percentiles
thre_corr=0.75        #threshold on correlation between candidate and reference series
min_num_ref=3         #below this number the correction is not performed
soft_min_num_ref=5    #below this number of references some inhomogeneous series are reintegrated 
max_num_ref=18        #above this number only the highest correlated references are selected

ELE=toupper(ele)
prefix=strsplit(ele,split='t')[[1]][2]

fq=ampl_bin
lq=100-fq
fb=(ampl_bin*3/2)          #first boundary
lb=100-fb                       #last boundary
nb=length(seq(fb,lb,ampl_bin))+2  #number of boundaries                    
endmonth=c(31,28,31,30,31,30,31,31,30,31,30,31)
length_ref=20
nq=nb-1 #number of intervals


if (nq%%2==0)
{
  fb=fb-ampl_bin/2
  lb=lb+ampl_bin/2  
}

### FUNCTIONS

source(paste0(sourcepath,'/hom_functions.R'))


##GET THE SERIES

setwd(inpath)


files.ori <- list.files(pattern=glob2rx("T*ori.txt")) #list of candidates: ORIginal version (withouth bdscore==FALSE portions)
files.parthom <- list.files(pattern=glob2rx("T*parthom.txt")) # list of candidates: PARTially HOMogenized during 1st iter
                                                         # though the parts that were removed because non homogenizable
                                                         # have been reintegrated for a 2ns attempt during 2nd iter
                                                         # withouth portions having bdscore==FALSE
files.cnd =c(files.ori, files.parthom) #list of all the candidates 
files.ref <- list.files(pattern=glob2rx("T*ref.txt"))  #list of files belonging to nearby stations
files=c(files.ori,files.parthom,files.ref)            #general list of files


series.cnd=lapply(files.cnd,read.table) 
names(series.cnd)=files.cnd

if (iter==1)
{
  load(paste0('../../BreakDetection/Output_',iter,'it_',ele,'/',prefix,'_breaks.RData'))
} 
if (iter==2)
{
  load(paste0('../../BreakDetection/Output_',iter,'it_',ele,'/',prefix,'_breaks.RData'))
}

if (length(files.ref)<3) ### the number of references is not sufficient to perform the adjustments
{
  for (s in 1:length(series.cnd)) #loop on the series of the same station
  {
    series=series.cnd[[s]]
    colnames(series)=c('year','month','day','series')
    series_name=names(series.cnd[s])
    ser_id=extract_ser_id(series_name)
    aux=substr(series_name,1,16)
    if (iter==1)  
    {
      output_mysql_s=paste0(aux,'mysqlhom1.txt')
      output_simple=paste0(aux,'hom1.txt')
    }
    if (iter==2) {
      output_mysql_s=paste0(aux,'mysqlhom2.txt')
      output_simple=paste0(aux,'hom2.txt')
    }
    chp=read.breaks(ser_id,prefix,iter,breaks)
    chp_y=chp$year
    chp_m=chp$month
    if (length(chp_y)!=0) ### if there are breaks but not enough references
    {
      if(debugflag>=1){print(paste0('Series ',ser_id,' is not homogeneous and cannot be homogenized before ', chp_y[length(chp_y)] ,' due to lack of references'))}
      donats=list() ### initialize list of donating subseries
      i=0
      ser_rest=series ### list of all the series
      for (k in 1:length(chp_y)) ### for each break, the part preceeding it is stored in the donating subseries
      {
        i=i+1
        ser_rest[ser_rest$year<chp_y[k] | (ser_rest$year==chp_y[k] & ser_rest$month<chp_m[k]),'series']=-999.9
      }
      basis=ser_rest ### the remaining part, after the last break, is the basis      
      firstday=as.Date(paste(series[1,1],series[1,2],series[1,3],sep='-'))
      lastday=as.Date(paste(series[dim(series)[1],1],series[dim(series)[1],2],series[dim(series)[1],3],sep='-'))
      final=basis
      colnames(final)=c('year','month','day','hom')
      final=complete(final,firstday,lastday) ### inserts all the missing days
      final$qc=8  ### qc temporarily set to 8, in a following step each day inherits the qc flag of its original version
      final$qca=8
      final[final$hom<(-90),c('qc','qca')]=9 ### set qc of missing data to 9

    }else ### the series is homogeneous and no references where downloaded
    {
      if(debugflag==1){print(paste0('Series ',ser_id,' is already homogeneous'))}
      firstday=as.Date(paste(series[1,1],series[1,2],series[1,3],sep='-'))
      lastday=as.Date(paste(series[dim(series)[1],1],series[dim(series)[1],2],series[dim(series)[1],3],sep='-'))
      final=series
      colnames(final)=c('year','month','day','hom')
      final=complete(final,firstday,lastday) ### inserts all the missing days
      final$qc=8  ### qc temporarily set to 8, in a following step each day inherits the qc flag of its original version
      final$qca=8
      final[final$hom<(-90),c('qc','qca')]=9 ### set qc of missing data to 9
     }
  }      
  out4mysql=data.frame(ser_id=ser_id,ser_date=final$year*10000+final$month*100+final$day,value=final$hom*10,
                       qc=final$qc,qca=final$qca,qcm=final$qca)    
  write.fwf(out4mysql,paste0(outpath,'/',output_mysql_s),colnames=FALSE,sep=', ',eol='\n',append=FALSE)
  write.fwf(final[c('year','month','day','hom')], paste0(outpath,'/',output_simple),colnames=FALSE,sep=' ',eol='\n',append=FALSE)
}else ### the number of references is sufficient to perform the adjustments
{
  series.ref=lapply(files.ref,read.table) ### read all the references
  names(series.ref)=files.ref
  
  ## DIVIDE REFERENCE SERIES INTO HOMOGENEOUS SUBSERIES
  
  list_rm_ref=NULL ### initialize the list of the reference that have been removed from the list
  
  for (ref_name in names(series.ref))
  {
    ref_id=extract_ser_id(ref_name)
    colnames(series.ref[[ref_name]])=c('year','month','day','refer')
    if (dim(series.ref[[ref_name]])[1]<thre_y*365) ### if the reference is too short
    {
      series.ref[[ref_name]]=NULL
      list_rm_ref=c(list_rm_ref,ref_name)
    }
  }
  
  ref_names=names(series.ref) ### list of references that are long enough to be considered
  refers=list() ### initialize the list of the references that will be used for a particular break
  i=0
  for (j in 1:length(series.ref))
  {
    ref_ser=series.ref[[j]]
    colnames(ref_ser)=c('year','month','day','refer')
    ref_name=names(series.ref[j])
    ref_id=extract_ser_id(ref_name)
    chp_ref=read.breaks(ref_id,prefix,iter,breaks) ### read the breaks of the considered reference series
    chp_y_ref=chp_ref$year
    chp_m_ref=chp_ref$month   
    if (length(chp_y_ref)==0) ### if there are no breaks the reference is stored as it is
    { 
      i=i+1
      refers[[i]]=ref_ser
      names(refers)[i]=paste0(ref_name,'-00')
    } else ### if there are breaks the reference series is divide into homogeneous subseries
    {
      ref_rest=ref_ser
      colnames(ref_rest)=c('year','month','day','refer')  
      for (k in 1:length(chp_y_ref))  ### loop on the breaks
      {
        ref_fin=ref_rest[ref_rest$year<chp_y_ref[k] | (ref_rest$year==chp_y_ref[k] & ref_rest$month<chp_m_ref[k]),] ### homogeneous subseries between two breaks
        ref_rest=ref_rest[ref_rest$year>chp_y_ref[k] | (ref_rest$year==chp_y_ref[k] & ref_rest$month>=chp_m_ref[k]),]  ### remaning part
        if (nrow(ref_fin)>0) ### if the subseries is not empty
        {
          i=i+1
          refers[[i]]=ref_fin ### stored in the list of reference ready to be used
          names(refers)[i]=paste0(ref_name,'-0',k) ### append a progressive code (-01,-02,-03, etc.) to the name of the reference
        }
      }
      if (nrow(ref_rest)>0) ### if the remaining part is not empty
      {
        i=i+1
        refers[[i]]=ref_rest
        names(refers)[i]=paste0(ref_name,'-0',k+1)
      }
    }
  }
  refers_all=refers ### all the references after the splitting
  
  for (s in 1:length(series.cnd)) #loop on the series of the same station
  {
    series=series.cnd[[s]] 
    colnames(series)=c('year','month','day','series')
    series_name=names(series.cnd[s])
    ser_id=extract_ser_id(series_name)
    if (ser_id==11847||ser_id==11848) ### these two series have very weird values between 1969 and 1972
    {      
      series=series[series$year<1969 |
                      series$year>1972 |
                      (series$year==1969 & series$month<=2) |
                      (series$year==1972 & series$month>=11),] 
    }
    
    first_day=as.Date(paste(series[1,1],series[1,2],series[1,3],sep='-'))
    last_day=as.Date(paste(series[dim(series)[1],1],series[dim(series)[1],2],series[dim(series)[1],3],sep='-'))
    
    #refers_ori=series.ref ### original references before the splitting
    refers_ori=lapply(refers, function(x) {colnames(x)=c('year','month','day','refer'); return(x)})
    aux=substr(series_name,1,16)
    if (iter==1)
    {
      output_mysql_a=paste0(aux,'mysqladj1.txt')
      output_mysql_s=paste0(aux,'mysqlhom1.txt')
      output_simple=paste0(aux,'hom1.txt')
    }
    if (iter==2)
    {
      output_mysql_a=paste0(aux,'mysqladj2.txt')
      output_mysql_s=paste0(aux,'mysqlhom2.txt')
      output_simple=paste0(aux,'hom2.txt') 
    }
    
    ## DIVIDE CANDIDATE SERIES INTO HOMOGENEOUS SUBSERIES
    
    fys=series[1,'year'] ### first year of the candidate
    lys=series[dim(series)[1],'year'] ### last year of the candidate
    chp=read.breaks(ser_id,prefix,iter,breaks)
    chp_y=chp$year
    chp_m=chp$month
    
    chp_y_old=chp_y
    if (debugflag==1){print('Breaks:'); print(chp_y)}
    chp_m=chp_m[which(chp_y>fys & chp_y<lys)] ### select only those breaks that are after the start of the series
    chp_y=chp_y[which(chp_y>fys & chp_y<lys)] ### 
    
    
    if (length(chp_y)==0) ### if there are no breaks
    {
      if (debugflag==1){print(paste('No breaks have been found in',ser_id))} 
      final=series 
      final$qc=8
      final$qca=8
      final[final$series<(-90),c('qc','qca')]=9
      colnames(final)=c("year", "month", "day", "hom","qc","qca")
      out4mysql=data.frame(ser_id=ser_id,ser_date=final$year*10000+final$month*100+final$day,value=final$hom*10,qc=final$qc,qca=final$qca,qcm=9)      
      
    }else ### if there is at least one break
    {
      donats=list() ### initialize list of donating subseries
      i=0
      ser_rest=series ### list of all the series
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
      if (dim(basis)[1]<thre_y*365)  ### if the basis is too short the last donating subseries is attached to it 
      {
        colnames(donats[[length(chp_y)]])=c('year','month','day','series')
        basis=rbind(donats[[length(chp_y)]],basis)
        donats[[length(chp_y)]]=NULL
      }
      basis_name=paste0(series_name,'-0',k+1)
      colnames(basis)=c('year','month','day','basis')
      basis_ser_id=extract_ser_id(series_name)
      benchmark=basis
      
      ### initializing tables that will be needed during the adjustment calculation
      overlap_dr=list() ### overlaps between donating and references
      qntl_drs=list()   ### quantiles of donating series in each overlap with each reference series
      bndrs_ds=list()   ### boundaries of the above-defined quantile bins
      adjs=list()      
      adjs01=list()
      adjs02=list()
      basis=benchmark   
      
      
      
      # intro=paste0("#",TNoTX," 2 iter adjustments for original series: ",ser_id)
      # hdr=paste0(format('#ser_id,',width=12,justify='right'),
      #            format('break_start,',width=12,justify='right'),
      #            format('break_stop,',width=12,justify='right'),
      #            format('ref_id,',width=12, justify='right'),
      #            format('month,',width=12, justify='right'),
      #            format('percentile,',width=12, justify='right'),
      #            format('lower_bnd,',width=12, justify='right'),
      #            format('upper_bnd,',width=12, justify='right'),
      #            format('adj,',width=12, justify='right'))
      
      
      if (length(donats)<1) ### if there are no donating (i.e. the only donating was attached to the basis)
      {
        if(debugflag==1){print(paste('Adjs calculation of',basis_ser_id,'not possible, since the latest segment is too short'))}        
        final=basis ### the basis is the final series
      }else
      {        
        for (d in length(donats):1) #loop (backwards) on donats
        {
          donat=donats[[d]]
          donat_ser_id=extract_ser_id(names(donats)[d])
          if (dim(donat)[1]<thre_y*365) ### if the donating subseries is not long enough
          {
            if (dim(donat)[1]>0)
            {
              fyd=donat[1,'year']
              lyd=donat[dim(donat)[1],'year']  
              if (debugflag==1){print(paste('Corrs calculation of',donat_ser_id,'from',fyd,'to',lyd,'is impossible due to the short period'))}
              aux=donat[c('year','month','day','donat')] ### auxiliary table
              colnames(aux)[4]='basis' 
              basis=rbind(aux[,c('year','month','day','basis')],basis) ### the uncorrected donating subseries is pre-appended to the basis without corrections
            }
          }else ### if the donating subseries is long enough
          {
            ##CHOOSE A SET OF REFERENCE SERIES REMOVING NON OVERLAPPING SERIES
            
            refers=refers_all                                        ### for each donating subseries the whole set of references is considered again
            fyb=basis[1,'year']                                      ### first year of the basis
            fmb=basis[1,'month']                                     ### first month of the basis
            lyb=basis[dim(basis)[1],1]                               ### last year of the basis
            fyd=donat[1,'year']                                      ### first year of the donor
            fmd=donat[1,'month']                                     ### first month of the donor
            lyd=donat[dim(donat)[1],'year']                          ### last year of the donor  
            lmd=donat[dim(donat)[1],'month']                         ### last month of the donor
            
            ### TO DO: condense refers, nyr, fyr and lyr in a single list
            
            nyr=lapply(refers, function(x) {x[dim(x)[1],1]-x[1,1]})  ### number of years in the references
            fyr=lapply(refers, function(x) {x[1,1]})                 ### fyrst year of the references
            lyr=lapply(refers, function(x) {x[dim(x)[1],1]})         ### last year of the references   
            i=1 
            while (i<=length(refers)) ### preliminary check whether the references have overlaps with both basis and donats
            {
              if (fyr[[i]]>lyb-thre_y+1 | lyr[[i]]<fyb+thre_y-1) ### if it doesn't have a long enough overlap
              {                                                 
                refers[[i]]=NULL                                 ### remove from the list of references
                fyr[[i]]=NULL
                lyr[[i]]=NULL
                nyr[[i]]=NULL
              }else                                              ### if it has a long enough overlap with the basis
              {
                if (fyr[[i]]>lyd-thre_y+1 | lyr[[i]]<fyd+thre_y-1)### then check the overlap with the donating (same criteria)
                { 
                  refers[[i]]=NULL                               ### remove from the list of refernces
                  fyr[[i]]=NULL
                  lyr[[i]]=NULL
                  nyr[[i]]=NULL
                }else 
                {
                  i=i+1                                          ### overlapping periods are long enough, proceed to next ref.
                }
              }
            } 
            
            refers=lapply(refers,function(x) crop.series(series=x,y=fyb,m=fmb,ny=length_ref)) ### crop the references (only 20 years of overlap are considered)
            nyr=lapply(refers, function(x) {x[dim(x)[1],1]-x[1,1]}) ### new lenghth in years
            fyr=lapply(refers, function(x) {x[1,1]})                ### new start
            lyr=lapply(refers, function(x) {x[dim(x)[1],1]})        ### new stop
            overlaps=lapply(refers,function(x) overlapize_ref(x,basis)) ### write down the overlapping periods and the values
            correlations=sapply(overlaps,function(x) cor(x$refer,x$basis,use='pairwise.complete.obs')) ### calculate correlations
            
            if (debugflag==1) ### write down correlations of the used references
            {
              value=sprintf("%06d",ser_id)
              refdetfile=paste0(outpath,'/Reference_Details/details_ref_',toupper(ele),'_',value,'_it',iter,'_',chp_y[d],'.txt')
              if (length(correlations)>0)
              {
                write(paste(names(correlations)[1],fyr[[1]],lyr[[1]],round(correlations[1],digits = 4)),refdetfile,append=FALSE)
                if (length(correlations)>1)
                {
                  for (i in 2:length(correlations))
                  {
                    write(paste(names(correlations)[i],fyr[[i]],lyr[[i]],round(correlations[i],digits = 4)),refdetfile,append=TRUE)
                  }
                }
              }else
              {write('No references',refdetfile,append=FALSE)}
            }
            
            ### REMOVAL OF REFERENCE WITH NA CORRELATION
            
            NAcorr=correlations[is.na(correlations)] ### NA correlations, due to insufficient data or other issues
            NAcorr0=NAcorr ### auxiliary table
            while(length(NAcorr)>0) ### remove references with NA correlations from the list of references
            {
              inac=which(is.na(correlations))[1]
              correlations=correlations[-inac]
              refers[[inac]]=NULL
              fyr[[inac]]=NULL
              lyr[[inac]]=NULL
              nyr[[inac]]=NULL
              NAcorr=correlations[is.na(correlations)]
            }
            
            ### REMOVAL OF REFERENCES WITH LOW CORRELATION
            
            correls=correlations 
            i=1
            j=1
            low_correl_ref=list()
            fylc=list()
            lylc=list()
            nylc=list()
            low_correls=NULL
            
            while (i<=length(refers)) ### loop on the references
            {
              if (correls[i]<thre_corr) ### if the correlation is too low
              { 
                low_correl_ref[[j]]=refers[[i]] ### add to the list of low correlated references
                low_correls[j]=correls[i]       ### add corresponding correlation to the list of low correlations
                names(low_correl_ref)[j]=names(refers)[i]
                names(low_correls)[j]=names(refers)[i]
                fylc[[j]]=fyr[[i]]
                names(fylc)[j]=names(refers)[i]
                lylc[[j]]=lyr[[i]]
                names(lylc)[j]=names(refers)[i]
                nylc[[j]]=nyr[[i]]
                names(nylc)[j]=names(refers)[i]
                refers[[i]]=NULL                ### remove it from the list of references
                fyr[[i]]=NULL
                lyr[[i]]=NULL
                nyr[[i]]=NULL
                correls=correls[-i]             ### remove from the list of correlations
                j=j+1
              }else
              {
                i=i+1
              }
            }
            
            i=1
            if (length(refers)>max_num_ref) ### if there are more than max_num_ref references, select the ones which are most correlated with the basis
            { 
              index_cr=order(correls,decreasing=TRUE)  ### sort the correlations
              mincr=correls[[index_cr[max_num_ref+1]]] ### first excluded one -> threshold
              if (debugflag==1)
              {
                write(paste(''),refdetfile,append=TRUE)
                write(paste('threshold corr      ','    ','    ',round(mincr,digits = 4)),refdetfile,append=TRUE)
              }
              
              i=1
              while (i<=length(refers)) ### check series and remove those with correlation below the threshold
              {
                if (correls[i]<=mincr)
                {
                  refers[[i]]=NULL
                  fyr[[i]]=NULL
                  lyr[[i]]=NULL
                  nyr[[i]]=NULL
                  correls=correls[-i]
                }else
                {
                  i=i+1
                }
              }
            }
            
            
            ### RECOVERY OF LOW CORRELATED SERIES (COMMENTED OUT)
            {
              # if (length(refers)<=5) ### add low correlated series if the subseries are not enough
              # {
              #   write( 'Not enough high correlated reference series. Some low correlated series will be used as references.',rprt,append=TRUE)
              #   print( 'Not enough high correlated reference series. Some low correlated series will be used as references.')
              #   num_ref_rec_need=5-length(refers)
              #   if (length(low_correls)>num_ref_rec_need)
              #   {
              #     index_lcr=order(low_correls,decreasing=TRUE)
              #     minlcr=low_correls[[index_lcr[num_ref_rec_need]]]
              #     i=1
              #     while (i<=length(low_correl_ref))
              #     {
              #       if (low_correls[i]<=minlcr-0.00001)
              #       {
              #         low_correl_ref[[i]]=NULL
              #         fylc[[i]]=NULL
              #         lylc[[i]]=NULL
              #         nylc[[i]]=NULL
              #         low_correls=low_correls[-i]
              #       }else
              #       {
              #         i=i+1
              #       }
              #     }
              #     refers=c(refers,low_correl_ref)
              #   }
              # }
            } #recovery of low correlated series (commented out)
            
            
            ### RECOVERY OF INHOMOGENEOUS REFERENCE SERIES (IF NEEDED)
            
            ref_ser_id=sapply(names(series.ref),extract_ser_id)
            
            if (length(refers)<soft_min_num_ref) 
            {
              if (debugflag==1){print('Not enough homogeneous reference series. Some non homogeneous series will be selected.')}
              
              ref_ser_id_ori=sapply(names(series.ref),extract_ser_id)
              ref_recover=series.ref[c(!ref_ser_id_ori%in%ref_ser_id)] ### the ref to be inspected are those that are not already selected
              nyrr=lapply(ref_recover, function(x) {x[dim(x)[1],1]-x[1,1]}) ### number of years in the recovery references
              fyrr=lapply(ref_recover, function(x) {x[1,1]})                ### first year
              lyrr=lapply(ref_recover, function(x) {x[dim(x)[1],1]})        ### last year
              i=1
              while (i<=length(ref_recover)) ### lop of the recovery references (same checks as above)
              {
                if (fyrr[[i]]>lyb-thre_y | lyrr[[i]]<fyb+thre_y)
                { 
                  ref_recover[[i]]=NULL
                  fyrr[[i]]=NULL
                  lyrr[[i]]=NULL
                  nyrr[[i]]=NULL
                }else
                {
                  if (fyrr[[i]]>lyd-thre_y | lyrr[[i]]<fyd+thre_y)
                  { 
                    ref_recover[[i]]=NULL
                    fyrr[[i]]=NULL
                    lyrr[[i]]=NULL
                    nyrr[[i]]=NULL
                  }else
                  {
                    i=i+1
                  }
                }
              }
              
              ref_recover=lapply(ref_recover,function(x) crop.series(series=x,y=fyb,m=fmb,ny=length_ref)) ### crop series (max 20 years are considered)
              nyrr=lapply(ref_recover, function(x) {x[dim(x)[1],1]-x[1,1]})
              fyrr=lapply(ref_recover, function(x) {x[1,1]})
              lyrr=lapply(ref_recover, function(x) {x[dim(x)[1],1]})
              
              i=1
              while (i<=length(ref_recover))
              {
                if (is.integer(fyrr) && is.integer(lyrr))
                {
                  if (!is.na(fyrr) && !is.na(lyrr))
                  {
                    if (fyrr[[i]]>lyb-thre_y | lyrr[[i]]<fyb+thre_y)
                    { 
                      ref_recover[[i]]=NULL
                      fyrr[[i]]=NULL
                      lyrr[[i]]=NULL
                      nyrr[[i]]=NULL
                    }else
                    {
                      if (fyrr[[i]]>lyd-thre_y | lyrr[[i]]<fyd+thre_y)
                      { 
                        ref_recover[[i]]=NULL
                        fyrr[[i]]=NULL
                        lyrr[[i]]=NULL
                        nyrr[[i]]=NULL
                      }else
                      {
                        i=i+1
                      }
                    }
                  }
                  else
                  {
                    ref_recover[[i]]=NULL
                    fyrr[[i]]=NULL
                    lyrr[[i]]=NULL
                    nyrr[[i]]=NULL
                  }
                }
                else
                {
                  ref_recover[[i]]=NULL
                  fyrr[[i]]=NULL
                  lyrr[[i]]=NULL
                  nyrr[[i]]=NULL
                }
              }
              
              overlaps_rec=lapply(ref_recover,function(x) overlapize_ref(x,basis))
              i=1
              while (i<=length(ref_recover))
              {
                if (nrow(overlaps_rec[[i]])<=thre_y*365)
                { 
                  overlaps_rec[[i]]=NULL
                  ref_recover[[i]]=NULL
                  fyrr[[i]]=NULL
                  lyrr[[i]]=NULL
                  nyrr[[i]]=NULL
                }else
                {
                  i=i+1
                }
              }
              correlations_rec=sapply(overlaps_rec,function(x) cor(x$refer,x$basis,use='pairwise.complete.obs'))
              if (debugflag==1)
              {          
                value=sprintf("%06d",ser_id)
                if (length(correlations_rec)>0)
                {
                  write(paste(''),ref_det_file,append=TRUE)
                  write(paste('Recovery of non-homogenized stations:'),refdetfile,append=TRUE)
                  write(paste(names(correlations_rec)[1],fyrr[[1]],lyrr[[1]],round(correlations_rec[1],digits = 4)),ref_det_file,append=TRUE)
                  if (length(correlations_rec)>1)
                  {
                    for (i in 2:length(correlations_rec))
                    {
                      write(paste(names(correlations_rec)[i],fyrr[[i]],lyrr[[i]],round(correlations_rec[i],digits = 4)),ref_det_file,append=TRUE)
                    }
                  }
                }else
                {write('No more references series to recover',refdetfile,append=TRUE)}
              }
              
              NAcorr=correlations_rec[is.na(correlations)]  ### remove of NA correlations
              NAcorr0=NAcorr
              while(length(NAcorr)>0)
              {
                inac=which(is.na(correlations_rec))[1]
                correlations_rec=correlations_rec[-inac]
                ref_recover[[inac]]=NULL
                fyrr[[inac]]=NULL
                lyrr[[inac]]=NULL
                nyrr[[inac]]=NULL
                NAcorr=correlations_rec[is.na(correlations_rec)]
              }
              
              correls=correlations_rec
              i=1
              j=1
              low_correl_ref=list()
              fylc=list()
              lylc=list()
              nylc=list()
              low_correls=NULL
              
              while (i<=length(ref_recover)) ### removal of low correlations
              {
                if (correls[i]<thre_corr)
                { 
                  low_correl_ref[[j]]=ref_recover[[i]]
                  low_correls[j]=correls[i]
                  names(low_correl_ref)[j]=names(ref_recover)[i]
                  names(low_correls)[j]=names(ref_recover)[i]
                  fylc[[j]]=fyrr[[i]]
                  names(fylc)[j]=names(ref_recover)[i]
                  lylc[[j]]=lyrr[[i]]
                  names(lylc)[j]=names(ref_recover)[i]
                  nylc[[j]]=nyrr[[i]]
                  names(nylc)[j]=names(ref_recover)[i]
                  ref_recover[[i]]=NULL
                  fyrr[[i]]=NULL
                  lyrr[[i]]=NULL
                  nyrr[[i]]=NULL
                  correls=correls[-i]
                  j=j+1
                }else
                {
                  i=i+1
                }
              }
              
              
              i=1
              num_ref_rec_need=soft_min_num_ref-length(refers) ### how many recovered references are needed to reach the soft limit?
              if (length(ref_recover)>num_ref_rec_need)
              {
                index_cr=order(correls,decreasing=TRUE) ### sort the recovered references by correlation
                mincr=correls[[index_cr[num_ref_rec_need+1]]]
                if (debugflag==1)
                {
                  write(paste(''),ref_det_file,append=TRUE)
                  write(paste('thresh corr rec  ','    ','    ',round(mincr,digits = 4)),ref_det_file,append=TRUE)
                }
                i=1
                while (i<=length(ref_recover)) ### remove those references that are not going to be recovered
                {
                  if (correls[i]<=mincr)
                  {
                    ref_recover[[i]]=NULL
                    fyrr[[i]]=NULL
                    lyrr[[i]]=NULL
                    nyrr[[i]]=NULL
                    correls=correls[-i]
                  }else
                  {
                    i=i+1
                  }
                }
              }              
              refers=c(refers,ref_recover) ### gather homogeneous references and inhomogeneous recovered references
              ref_ser_id=sapply(names(refers),extract_ser_id)
            }
            
            if (length(refers)<min_num_ref) ### if even after the recovery the number of refs is below the hard limit
            {
              if (debugflag==2){print(paste('Adjustment calculation for series',ser_id,' from ',fyd,' to ',lyd,' not possible due to lack of references'))}
              if (debugflag==1){print(paste0('Corrs calculation of ',donat_ser_id,' from ',fyd,' to ',lyd,' has not enough references'))}
              aux=donat[c('year','month','day','donat')]
              colnames(aux)[4]='basis'
              aux$basis=-999.9 ###non-homogenizable values are changed to missing data
              basis=rbind(aux[,c('year','month','day','basis')],basis)   ### the donating series is pre-appended to the basis
            }else                            ### if there are enough references
            {              
              if (debugflag==1)
              {
                if(debugflag==1){print('List of references:')}
                if(debugflag==1){print(names(refers))}
              }
            
            
            
            ###IDENTIFY OVERLAPPING PERIODS AND CALCULATE QUANTILES
            
            
           
            
         
              overlap_br=lapply(refers, function(x) overlapize_ref(x,basis));   ### overlap between basis and reference
              qntl_brs=lapply(overlap_br, function(x) corr_qntl_ref_3m(x,nb,nq,fq,lq,ampl_bin,thre_y)) ### quantiles of basis and references
                   
            donat_ser_ids=ser_id
            overlap_dr[[d]]=list()
            qntl_drs[[d]]=list()
            bndrs_ds[[d]]=list()
            adjs[[d]]=list()
            donat_name=names(donats)[d]
            donat_ser_id=ser_id
            if (debugflag>=1){print(paste('Adjustment calculation of',donat_ser_id,'from',fyd,'to',lyd))}
            
            overlap=merge(basis,donat)  ### overlap of basis and donat
            error=0
            if (dim(overlap)[1]!=0)
            {
              error=1
              print('There is an overlap! Error somewhere')
            }
            if (error==0)
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
                adjs[[d]]=list()
                
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
                    
                    
                    ### WRITE IN THE STATION FILE
                    
                    if (lyd%%4==0 & (lyd%%100!=0 | lyd%%400==0)) {endmonth[2]=29} ### leap years
                    ## PREPARING ADJUSTMENTS FILE
                    ### the below for loops prepares the lines with the adjustment details for the mysql table perc_adj_tx
                    ### for quantile bins on the left half: the upper boundary is decreased by 0.1 to avoid overlapping with the lower boundary of the next quantile
                    ### for quantile bins on the right half: the lower boundary is increased by 0.1
                    ### quantile bin of the median is left unchanged
                    for (m in 1:12)
                    {
                      for (q in 1:nq)
                      {
                        if(q<10) ### q spans between 1 and 19, so this is for quantiles below the median
                        {
                          if(!is.na(adjs[[d]][[j]][q,m])) ### if the adjustment is not NA
                          {
                            row=paste(format(ser_id,width=12),                                                         ### ser_id
                                      format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),          ### start
                                      format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'), ### stop
                                      format(ref_ser_id[[j]],width=12),                                                ### ref_id
                                      format(m,width=12),                                                              ### month
                                      format(as.numeric(rownames(adjs[[d]][[j]])[q]),width=12),                        ### quantile
                                      format(as.numeric(bndrs_ds[[d]][[j]][q,m]*10),width=12),                         ### lower boundary
                                      format((as.numeric(bndrs_ds[[d]][[j]][q+1,m]-0.1)*10),width=12),                 ### upper boundary (decreased by 0.1)
                                      paste0(format(adjs[[d]][[j]][q,m]*10,width=12)),sep=',')                         ### adjustment  
                          }else
                          {
                            if(abs(bndrs_ds[[d]][[j]][q,m])>100 & abs(bndrs_ds[[d]][[j]][q+1,m])>100) ### if the boundaries are  -999.9
                            {
                              row=paste(format(ser_id,width=12),
                                        format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                        format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                        format(ref_ser_id[[j]],width=12),
                                        format(m,width=12),
                                        format(as.numeric(rownames(adjs[[d]][[j]])[q]),width=12),
                                        format(-9999,width=12),
                                        format(-9999,width=12),
                                        paste0(format(-9999,width=12)),sep=',')
                            }else ### if the adjustment is NA/-999.9 but the boundaries are normal values
                            {
                              row=paste(format(ser_id,width=12),
                                        format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                        format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                        format(ref_ser_id[[j]],width=12),
                                        format(m,width=12),
                                        format(as.numeric(rownames(adjs[[d]][[j]])[q]),width=12),
                                        format(as.numeric(bndrs_ds[[d]][[j]][q,m])*10,width=12),
                                        format((as.numeric(bndrs_ds[[d]][[j]][q+1,m])-0.1)*10,width=12),
                                        paste0(format(-9999,width=12)),sep=',')
                            }
                          }
                        }
                        if(q==10) ### if we are working with the median
                        {
                          if(!is.na(adjs[[d]][[j]][q,m]))
                          {
                            row=paste(format(ser_id,width=12),
                                      format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                      format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                      format(ref_ser_id[[j]],width=12),
                                      format(m,width=12),
                                      format(as.numeric(rownames(adjs[[d]][[j]])[q]),width=12),
                                      format((as.numeric(bndrs_ds[[d]][[j]][q,m]))*10,width=12), 
                                      format((as.numeric(bndrs_ds[[d]][[j]][q+1,m]))*10,width=12), ###  upper quantile increased
                                      paste0(format(adjs[[d]][[j]][q,m]*10,width=12)),sep=',')
                          }else
                          {
                            if(abs(bndrs_ds[[d]][[j]][q,m])>100 & abs(bndrs_ds[[d]][[j]][q+1,m])>100)
                            {
                              row=paste(format(ser_id,width=12),
                                        format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                        format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                        format(ref_ser_id[[j]],width=12),
                                        format(m,width=12),
                                        format(as.numeric(rownames(adjs[[d]][[j]])[q]),width=12),
                                        format(-9999,width=12),
                                        format(-9999,width=12),
                                        paste0(format(-9999,width=12)),sep=',')
                            }else
                            {
                              row=paste(format(ser_id,width=12),
                                        format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                        format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                        format(ref_ser_id[[j]],width=12),
                                        format(m,width=12),
                                        format(as.numeric(rownames(adjs[[d]][[j]])[q]),width=12),
                                        format((as.numeric(bndrs_ds[[d]][[j]][q,m]))*10,width=12),
                                        format((as.numeric(bndrs_ds[[d]][[j]][q+1,m]))*10,width=12),
                                        paste0(format(-9999,width=12)),sep=',')
                            }
                          }
                        }
                        if(q>10)  ### for the right half of the sequence, same procedure as above
                        {
                          if(!is.na(adjs[[d]][[j]][q,m]))
                          {
                            row=paste(format(ser_id,width=12),
                                      format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                      format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                      format(ref_ser_id[[j]],width=12),
                                      format(m,width=12),
                                      format(as.numeric(rownames(adjs[[d]][[j]])[q]),width=12),
                                      format((as.numeric(bndrs_ds[[d]][[j]][q,m])+0.1)*10,width=12),
                                      format(as.numeric(bndrs_ds[[d]][[j]][q+1,m])*10,width=12),
                                      paste0(format(adjs[[d]][[j]][q,m]*10,width=12)),sep=',')
                          }else
                          {
                            if(abs(bndrs_ds[[d]][[j]][q,m])>100 & abs(bndrs_ds[[d]][[j]][q+1,m])>100)
                            {
                              row=paste(format(ser_id,width=12),
                                        format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                         format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                          format(ref_ser_id[[j]],width=12),
                                          format(m,width=12),
                                          format(as.numeric(rownames(adjs[[d]][[j]])[q]),width=12),
                                          format(-9999,width=12),
                                          format(-9999,width=12),
                                          paste0(format(-9999,width=12)),sep=',')
                              }else
                              {
                                row=paste(format(ser_id,width=12),
                                          format(paste0(fyd,sprintf("%02d", fmd),'01'),width=12,justify='right'),
                                          format(paste0(lyd,sprintf("%02d", lmd),endmonth[lmd]),width=12,justify='right'),
                                          format(ref_ser_id[[j]],width=12),
                                          format(m,width=12),
                                          format(as.numeric(rownames(adjs[[d]][[j]])[q]),width=12),
                                          format((as.numeric(bndrs_ds[[d]][[j]][q,m])+0.1)*10,width=12),
                                          format(as.numeric(bndrs_ds[[d]][[j]][q+1,m])*10,width=12),
                                          paste0(format(-9999,width=12)),sep=',')
                              }
                            }
                          }
                          write(row,paste0(outpath,'/',output_mysql_a),append=TRUE)
                        }
                      }
                    }     
                  }
                }
              }else #donating series is a subset of the base series
              {
                print(paste('Error in series',ser_id,'.Overlapping period between basis and donating!'))
              }
            
            
            
              ### ADJUSTMENT OF THE SERIES
              homdon=donat[c('year','month','day','donat')] ### initializing the homogenized version of the donating series    
              homdon$hom=apply(donat,1,function(x) adjust(x,bndrs_ds[[d]],nb,adjs[[d]]))                       
              homdon$diff=homdon$hom-homdon$donat
              colnames(homdon)[5]='basis'
	      basis=rbind(homdon[,c('year','month','day','basis')],basis)  ### the adjusted subseries is preappended to the basis
            								 ### the new longer basis is used for the homogenization of the previous portion
            } # end of "if there are enought reference series for this donating series"
          } #end of "if the donating series is long enough"
        }#end of the loop on the donating series
        final=basis 
      }   
      
      
      
      colnames(final)=c('year','month','day','hom')
      fs=as.Date(paste(final[1,'year'],final[1,'month'],final[1,'day'],sep='-'))
      ls=as.Date(paste(final[nrow(final),'year'],final[nrow(final),'month'],final[nrow(final),'day'],sep='-'))
      final=complete(final,fs,ls) ### complete the series with eventual missing days
      final$qc=8
      final$qca=8
      final[final$hom<(-90),c('qc','qca')]=9 
      out4mysql=data.frame(ser_id=ser_id,ser_date=final$year*10000+final$month*100+final$day,value=final$hom*10,qc=final$qc,
                           qca=final$qca,qcm=final$qc) ### mysql version
    }    
    
    
    out4storage=final[c('year','month','day','hom')]
    write.fwf(out4storage,paste0(outpath,'/',output_simple),colnames=FALSE,sep=' ',eol='\n',append=FALSE)
    write.fwf(out4mysql,paste0(outpath,'/',output_mysql_s),colnames=FALSE,sep=', ',eol='\n',append=FALSE)
    
    
  }#end of if(length(files.ref)<=3)
}
