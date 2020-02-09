### FUNCTIONS USED IN THE HOMOGENIZATION PROCESSESS
### THE FUNCTIONS ARE SORTED IN ALPHABETICAL ORDER






### calculation of the boundaries of the quantile bins for an array of values. Example: for quantile 25 calculates values of quantile 22.5 and 27.5 
bndrs_calc <- function(vec,nq1,ampl_bin1,thre_y1)
{
  if (length(vec[!is.na(vec)])>=thre_y1*90)
  {
    ll=length(vec[!is.na(vec)])/(100/ampl_bin1)                    #number of data in each interval
    #to calculate values of the bndrs you need to know the center of each interval
    card_qnt1=round(seq(0,nq1+1)*ll)
    vec1_table=data.frame(value=sort(vec[!is.na(vec)]),cat=NA)
    for (b in (1:(nq+1)))
    {
      vec1_table[((card_qnt1[b]+1):card_qnt1[b+1]),'cat']=b
    }
    val_bnd1=tapply(vec1_table$value,vec1_table$cat,wise_median)
    val_bnd1[1]=-999.9
    val_bnd1[nq1+1]=999.9
  }else
  {
    val_bnd1=rep(-999.9,nq1+1)
  }
  val_bnd1=round(val_bnd1,digits = 1)
  return(val_bnd1)
}

### Calculation of boundaries: takes as input a time table and calls bndrs_calc.
bndrs_d_3m <- function(overlap,nq1,fb1,lb1,ampl_bin1,thre_y1)
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
      aux_qnt[[m]]=bndrs_calc(vec,nq1,ampl_bin1,thre_y1)
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
    row.names(bndrs)=c(0,seq(fb1,lb1,ampl_bin1),100)
  }
  return(bndrs)
}

boxPlots <- function(table,index,time_int,printToFile='TRUE',lab_size=1,plot_title=TRUE,plotmin=-0.8,plotmax=1.2)
{
  nam=paste0('boxplot_trends_',TNoTX,'_',index,'_',time_int,'_',index,'.png')
  width=7
  if(printToFile==TRUE){  png(nam, width = width, height = width/2,res=600, units='in')}
  maintext=ifelse(plot_title,paste('Trends distribution (interval:',time_int,')',index),'')
  par(mar=c(3,3,4,3))
  boxplot(table[,c('ORI','HOM1','HOM2')],
          horizontal = TRUE,
          #ylab=c('Steps of homogenization'),
          ylab=' ',
          #xlab='trends[C/decade]',
          xlab=' ',
          yaxt='n',
          ylim=c(plotmin,plotmax),
          axes=FALSE,
          #xlim=c(min(table[,c(2,4)])-0.1,max(table[,c(2,4)])+0.1),
          main=maintext,
          col=c(rgb(0,0,0,0.2),rgb(0,0.2,1,0.2),rgb(1,0,0,0.2)))
  
  title(ylab=c('Steps of homogenization'),cex.lab=lab_size,line=2)
  title(xlab='trends[C/decade]',cex.lab=lab_size,line=2)
  axis(side=1,cex.axis=lab_size)
  axis(side=2,cex.axis=lab_size,at=c(1,2,3),labels=c('ORI','HOM1', 'HOM2'))
  abline(v=seq(-10,10,0.05),lwd=0.8,col='grey',lty=3)
  abline(v=seq(-10,10,0.1),lwd=1.5,col='darkgrey',lty=1)
  box()
  legend('topleft',
         cex=0.8,
         legend=c('hom2','hom1','ori'),
         inset=c(0.03,0.03,0.3),
         pch=c(15,15,15),
         col=c(rgb(1,0,0,0.2),rgb(0,0.2,1,0.2),rgb(0,0,0,0.2)))
  if(printToFile==TRUE){dev.off()}
}

boxPlots4Counts <- function(table,index,time_int,printToFile='TRUE',lab_size=1)
{
  nam=paste0('boxplot_trends_',TNoTX,'_',index,'_',time_int,'_',index,'.png',lab_size=1)
  if(printToFile==TRUE){  png(nam, width = 600, height = 300)}
  
  boxplot(table[,c('ORI','HOM1','HOM2')],
          horizontal = TRUE,
          #ylab=c('Steps of homogenization'),
          ylab=' ',
          #xlab='trends[C/decade]',
          xlab=' ',
          #ylim=c(-3,3),
          xact='n',
          yact='n',
          axes=FALSE,
          ylim=c(min(table[,c(2,4)],na.rm=TRUE)-1,max(table[,c(2,4)],na.rm=TRUE)+1),
          main=paste('Trends distribution (interval:',time_int,')',index, TNoTX ),
          col=c(rgb(0,0,0,0.2),rgb(0,0.2,1,0.2),rgb(1,0,0,0.2)))
  title(ylab=c('Steps of homogenization'),xlab='trends[C/decade]',cex.lab=lab_size)
  axis(side=1,cex.axis=lab_size)
  axis(side=2,cex.axis=lab_size)
  legend('topleft',
         legend=c('hom2','hom1','ori'),
         inset=c(0.03,0.03,0.3),
         pch=c(15,15,15),
         col=c(rgb(1,0,0,0.2),rgb(0,0.2,1,0.2),rgb(0,0,0,0.2)))
  if(printToFile==TRUE){  dev.off()}
}

check_av <- function(vec,thre_y1)
{
  if (length(vec)>thre_y1*30)
  {control=1}else
  {control=0}
  return(control)
}

# check_latlon_bm <- function(ele,name,breaks)
# {
#   name=paste0(ele,'_',extract_ser_id(name,base=3)+3000000,'.tst')
#   if (length(breaks[[name]])>0)
#   {
#   latlon=list()
#   latlon$lat=breaks[[name]]$latitude
#   latlon$lon=breaks[[name]]$longitude
#   }
#   name=paste0(ele,'_',extract_ser_id(name,base=3)+3000000,'.ser')
#   if (length(breaks[[name]])>0)
#   {
#     latlon=list()
#     latlon$lat=breaks[[name]]$latitude
#     latlon$lon=breaks[[name]]$longitude
#   }
#   return(latlon)
# }

check_month_avail <- function(overlap,thre_y1)
{
  control=0
  control_month=tapply(overlap$basis,overlap$month,function(x) check_av(x,thre_y1))
  if (sum(control_month)>0)
  {
    control=1
  }
  return(control)
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

clrsFun <- function(numclr=6)
{
  clrs=rep(c('red','blue','gold1','green','purple','darkorange1','seagreen','pink','cyan','magenta','brown','black','grey',
             'darkolivegreen','violetred','cyan4','saddlebrown','dodgerblue4','hotpink3','cornsilk3'),time=10)
  clrs=clrs[1:numclr]
  return(clrs)
}


### calculates the matrix of quantile sequences (month x quantile) for basis/donat and for reference, then takes the difference
corr_qntl_ref_3m <- function(overlap,nb1,nq1,fq1,lq1,ampl_bin1,thre_y1)
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
      aux_qnt[[m]]=qntls_3m(vec,nb1,nq1,ampl_bin1,thre_y1)
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
       aux_qnt[[m]]=qntls_3m(vec,nb1,nq1,ampl_bin1,thre_y1=5)
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
    aux_r_qnt[[m]]=qntls_3m(vec,nb1,nq1,ampl_bin1,thre_y1)
  }
  qnt=sapply(aux_qnt, unlist)
  colnames(qnt)=c('01','02','03','04','05','06','07','08','09','10','11','12')
  row.names(qnt)=seq(fq1,lq1,ampl_bin1)
  r_qnt=sapply(aux_r_qnt, unlist)
  colnames(r_qnt)=c('01','02','03','04','05','06','07','08','09','10','11','12')
  row.names(r_qnt)=seq(fq1,lq1,ampl_bin1)
  corr=qnt-r_qnt
  return(corr)
}

### application of the adjustment to the series
adjust <- function(dat,bndrs_ds1,nb1,corrs1)
{
  m=as.numeric(dat['month'])
  value=as.numeric(dat['donat'])
  perc=c(NA,NA)
  add=c(NA,NA)
  exst=NA
  check0=sapply(bndrs_ds1, function(x)  if (x[10,m]<(-99)) {return(0)} else {return(1)}) ### check if there are missing data in the boundaries in the months
  check0=sum(check0) ### how many months have non-missing data?
  if (check0>=3) ### at least 3 references?
  {
    for (j in 1: length(refers))
    {
      vec_b=bndrs_ds1[[j]][,m]
      if (!is.na(vec_b[2]) && vec_b[2]>-90)
      {
        if(value<vec_b[2]) ### inspect in which quantile bin is the value that has to be adjusted
        {
          perc[j]=5
        }else if (value>vec_b[nb1-1])
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
        
        add[j]=corrs1[[j]][which(rownames(corrs1[[j]])==perc[j]),m] ### what is the adjustment corresponding the quantile and the month?
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

adjust_blend <- function(dat,basis_ser_id,donat_ser_ids,bndrs_ds,nb,corrs) ### same as above but for blended series (it needs the blend_id as input)
{
  blend_id=as.numeric(dat['blend_id'])
  value=as.numeric(dat['blend'])
  if (blend_id != basis_ser_id)
  {  
    i=which(donat_ser_ids==blend_id)
    if (length(corrs[[i]])>=3)
    {
      m=as.numeric(dat['month'])
      perc=c(NA,NA)
      add=c(NA,NA)
      exst=NA
      check=sapply(bndrs_ds[[i]], function(x)  if (x[10,m]<(-99) | is.na(x[10,m])) {return(0)} else {return(1)})
      check=sum(check)
      if (check>=3)
      {
        for (j in 1: length(bndrs_ds[[i]]))
        {
          vec_b=bndrs_ds[[i]][[j]][,m]
          if (!is.na(vec_b[2]) && vec_b[2]>-90)
          {
            if(value<vec_b[2])
            {
              perc[j]=5
            }else if (value>vec_b[nb-1])
            {
              perc[j]=95
            }else
            {
              low_b=max(as.numeric(labels(which(vec_b<=value))))
              high_b=min(as.numeric(labels(which(vec_b>=value))))
              if (low_b < high_b & high_b-low_b==5)
              {
                perc[j]=(low_b+high_b)/2
              } else if(low_b==high_b)
              {
                perc[j]=(abs(low_b-50)-2.5)*sign(low_b-50)+50
              } else
              {
                if(low_b<50 & high_b<50)
                {
                  perc[j]=low_b-2.5
                } else if(low_b>50 & high_b>50)
                {
                  perc[j]=high_b+2.5
                } else
                {
                  perc[j]=50
                }
              }
            }
            add[j]=corrs[[i]][[j]][which(rownames(corrs[[i]][[j]])==perc[j]),m]
            exst[j]=value+add[j]
          }
        }
      }
      if (length(unlist(exst)[!is.na(unlist(exst))])<=3)
      {
        hom=-999.9
      } else
      {
        hom=wise_median(unlist(exst))
        if(value>hom)
        {
          hom=round(hom+0.01,digits=1)
        }else if(value<hom)
        {
          hom=round(hom-0.01,digits=1)
        }
      } 
    }else
    {
      hom=-999.9
    }
  }else
  {
    hom=value
  }
  
  return(hom)
}

complete <- function(series1,fs1,ls1,missing=-999.9) ### given series, first day and last day, it completes it adding missing values
{
  series1$idate=as.Date(paste(series1$year,series1$month,series1$day,sep='-'))
  dates=seq.Date(fs1,ls1,1)
  removed=dates[-which(dates%in%series1$idate)]
  if (length(series1$idate)<length(dates))
  {
    adding=data.frame(year=as.numeric(format.Date(removed,'%Y')),
                      month=as.numeric(format.Date(removed,'%m')),
                      day=as.numeric(format.Date(removed,'%d')),
                      hom=missing,idate=removed)
    series1=rbind(series1,adding)
    series1=series1[order(series1$idate),] 
  }
  series1$idate=NULL
  return(series1)
}

create_pdf <- function(vec)
{
  btm=floor(min(vec-0.5))+0.5
  top=ceiling(max(vec+0.5))-0.5
  sub_vec=data.frame(bin=seq(btm-4.5,top+4.5),num=NA)
  sub_vec[c(1,2,3,4,5,
            dim(sub_vec)[1]-4,dim(sub_vec)[1]-3,dim(sub_vec)[1]-2,dim(sub_vec)[1]-1,dim(sub_vec)[1]),'num']=0
  for (i in (btm+0.5):(top-0.5))
  {
    sub_vec[sub_vec$bin==i,'num']=length(vec[vec>=i-0.5 & vec<i+0.5])
  }
  return(sub_vec)
}

crop.series <- function(series,y,m,y2=0,m2=0,ny) ### takes a series and crops it withing a given range (ny, e.g. 20) around given year (y, e.g. 1968) and month (m and m2) (e.g. from 1948 to 1988 )
{                                                ### y2 and m2 are used during the blending when the two series don'thave any overlap (e.g. one stops in 1970 and the other starts in 1973, meaning series are cropped from 1950 to 1993)
  if (y2==0)       ### if there is only one break                              
  {
    mm=(y-2000)*12+m 
    fmm=(series[1,'year']-2000)*12+series[1,'month']             #number of the first month
    lmm=(series[dim(series)[1],'year']-200)*12+series[1,'month'] #number of the last month
    if(mm-fmm>ny*12)
    {
      series1=series[series$year>(y-ny) | (series$year==(y-ny) & series$month>=m),]
    }else
    {series1=series}
    if(lmm-mm>((ny*12)-1))
    {
      series2=series1[series1$year<(y+ny) | (series1$year==(y+ny) & series1$month<m),]
    }else
    {series2=series1}
  }else            ### if there are two breaks
  {
    y1=y
    fy=min(y1,y2)
    ly=max(y1,y2)
    mm_left=(fy-2000)*12
    mm_right=(ly-2000+1)*12
    fmm=(series[1,'year']-2000)*12+series[1,'month'] #number of the first month
    lmm=(series[dim(series)[1],'year']-2000)*12+series[dim(series)[1],'month'] #number of the last month
    if(mm_left-fmm>ny*12)
    {
      series1=series[series$year>=(fy-ny),]
    }else
    {series1=series}
    if(lmm-mm_right>ny*12)
    {
      series2=series1[series1$year<=(ly+ny),]
    }else
    {series2=series1}
  }
  
  return(series2)
}

dens_lin_bm <- function(vec,nhom,bin,xmin,clr,
                        lwd=5,bw=0,plotlim=c(-9999,9999))
{ 
  a=(nhom+1)/bin
  b=1+(nhom/2)-(nhom+1)*xmin/bin
  c=length(vec)*bin
  if(length(vec[!is.na(vec)])>3)
  {
    vec=vec[!is.na(vec)]
    if (bw==0){
      plotx=density(vec)$x
      ploty=density(vec)$y
            ploty=ploty[plotx>(plotlim[1]) & plotx<(plotlim[2])]
            plotx=plotx[plotx>(plotlim[1]) & plotx<(plotlim[2])]
    }else
    {
      plotx=density(vec,bw=bw)$x
      ploty=density(vec,bw=bw)$y
      ploty=ploty[plotx>(plotlim[1]) & plotx<(plotlim[2])]
     plotx=plotx[plotx>(plotlim[1]) & plotx<(plotlim[2])]
     }
    lines(plotx*a+b,ploty*c,col=clr,lwd=lwd)
  abline(v=plotlim[1]*a+b,lwd=2)
  abline(v=plotlim[2]*a+b,lwd=2)
      }
}

extract_ser_id <- function(string,base=0) ### extract series id from conventional filename (e.g. TG_serid_000435_ori.txt)
{
  aux=strsplit(string, '_')[[1]][3]
  aux=as.numeric(aux)
  ser_id=aux-base*1000000
  return(ser_id)
}

extract_ser_id_old <- function(string,base=1) ### extract series id from old conventional filename (e.g. TG_1000435.ori)
{
  aux=strsplit(string, '_')[[1]][2]
  aux=substr(aux,1,7)
  aux=as.numeric(aux)
  ser_id=aux-base*1000000
  return(ser_id)
}

index_benchmark <- function(table,rm=TRUE) ### homogenization indicator used in the comparison of homogenization methods
{
  categories=colnames(table)
  categories=categories[-which(categories=='tst' | categories=='ser' | categories=='year' | categories=='ori' | categories=='inh')]
  colnames(table)[which(colnames(table)=='tst')]='inh'
  colnames(table)[which(colnames(table)=='ser')]='ori'
  
  for(cat in categories)
  {
    if (is.na(strsplit(cat,'\\.')[[1]][2]))
    {
      table[paste0(cat,'.ind')]=NA
      for (j in 1:length(table[,cat]))
      {
        if (!(any(is.na(table[j,c('ori','inh', cat)]))))
        {
          if (abs(table[j,'inh']-table[j,'ori'])>=0.009)
          {
            table[j,paste0(cat,'.ind')]=1-(table[j,cat]-table[j,'ori'])/(table[j,'inh']-table[j,'ori'])
          }
        }
      }
      if (rm)
      {
      table[,paste0(cat,'.ind','.rm')]=running_mean_2(table[,c('year',paste0(cat,'.ind'))])
      }
    }
  }
  return(table)
}


### FOLLOWING FUNCTIONS CALCULATE INDICES FROM GIVEN SERIES

ind_annmean <- function (series,fy=1700,ly=2100,rrm=11,period=365)
{
  annmean=data.frame(year=fy:ly,annmean=NA,annmean.rm=NA)
  old_colname=colnames(series)[4]
  for  (y in fy:ly)
  {
    if (length(series[series$year==y & !is.na(series[,4]),4])>period*0.82)
    {
      annmean[annmean$year==y,'annmean']=mean(series[series$year==y,4])
    }
  }
  annmean[,'annmean.rm']=running_mean_2(annmean[,1:2],rrm)
  return(annmean)
}

ind_frstdays <- function (series,fy=1700,ly=2100,rrm=11)
{
  fd=data.frame(year=fy:ly,fd=NA)
  for  (y in fy:ly)
  {
    if (length(series[series$year==y & series$month%in%c(1,2,12)& !is.na(series[,4]),4])>72)
    {
      fd[fd$year==y,'fd']=dim(series[series$year==y & series[,4]<0 & !is.na(series[,4]),])[1]
    }
  }
  fd[,3]=running_mean_2(fd,rrm)
  return(fd[,2:3])
}

ind_smrdays <- function (series,fy=1700,ly=2100,rrm=11)
{
  sd=data.frame(year=fy:ly,sd=NA)
  for  (y in fy:ly)
  {
    if (length(series[series$year==y & series$month%in%c(6,7,8)& !is.na(series[,4]),4])>72)
    {
      sd[sd$year==y,'sd']=dim(series[series$year==y & series[,4]>24.99 & !is.na(series[,4]),])[1]
    }
  }
  sd[,3]=running_mean_2(sd,rrm)
  return(sd[,2:3])
}

ind_qntl <- function (series, qntl, fy=1700,ly=2100,rrm=11,season='ann',period=90) ### calculated the yearly quantile related to a certain series, can also be seasonal
{
  range=1
  month_range=seq(1,12,1)
  if (season=='djf')
  {
    month_range=c(1)
    series[series$month==12,'year']=series[series$month==12,'year']+1
  }
  if (season=='mam')
  {
    month_range=c(4)
  }
  if (season=='jja')
  {
    month_range=c(7)
  }
  if (season=='son')
  {
    month_range=c(10)
  }
  if (is.numeric(season) )
  {
    if (as.numeric(season) %in% seq(1,12))
    {
      month_range=season
    range=0
    }
  }
  
    
  qq=data.frame(year=fy:ly,qq=NA)
  for  (y in fy:ly)
  {
    check_m=0
    if (length(series[series$year==y & !is.na(series[,4]),4])>=period*4*0.8)
    {
      for (m in month_range)
      {
        if (length(series[series$year==y & series$month%in%c(one2twelve(m-range),m,one2twelve(m+range)) & !is.na(series[,4]),4])>=period*0.8)
        {
          #print(m)
          check_m=check_m+1
        }
      }
    }
    if (season=='ann' & check_m==12)
    {
      ll=length(series[series$year==y & !is.na(series[,4]),4])
      qq[qq$year==y,2]=wise_median(sort(series[series$year==y & !is.na(series[,4]),4])[round(ll/100*(qntl-2.5)):round(ll/100*(qntl+2.5))],digits=1)
    }
    if (season!='ann' & check_m==3)
    {
      ll=length(series[series$year==y & !is.na(series[,4]),4])
      qq[qq$year==y,2]=wise_median(sort(series[series$year==y & !is.na(series[,4]),4])[round(ll/100*(qntl-2.5)):round(ll/100*(qntl+2.5))],digits=1)
    }
    if (is.numeric(season) & check_m==1)
    {
      ll=length(series[series$year==y & !is.na(series[,4]),4])
      qq[qq$year==y,2]=wise_median(sort(series[series$year==y & !is.na(series[,4]),4])[round(ll/100*(qntl-2.5)):round(ll/100*(qntl+2.5))],digits=1)
    }
  }
  qq[,'qq_rm']=running_mean_2(qq,rrm=11)
  
  return(qq)
}

ind_ty05 <- function (series,fy=1700,ly=2100,rrm=11)
{
  tn05=data.frame(year=fy:ly,tnn=NA)
  for  (y in fy:ly)
  {
    if (length(series[series$year==y & series$month%in%c(1,2,12) & !is.na(series[,4]),4])>72)
    {
      ll=length(series[series$year==y & !is.na(series[,4]),4])
      tn05[tn05$year==y,2]=sort(series[series$year==y & !is.na(series[,4]),4])[ll/100*5]
    }
  }
  tn05[,3]=running_mean_2(tn05,rrm)
  return(tn05[,2:3])
}

ind_ty95 <- function (series,fy=1700,ly=2100,rrm=11)
{
  tn95=data.frame(year=fy:ly,tnn=NA)
  for  (y in fy:ly)
  {
    if (length(series[series$year==y & series$month%in%c(6,7,8) & !is.na(series[,4]),4])>72)
    {
      ll=length(series[series$year==y & !is.na(series[,4]),4])
      tn95[tn95$year==y,2]=sort(series[series$year==y & !is.na(series[,4]),4])[ll/100*95]
    }
  }
  tn95[,3]=running_mean_2(tn95,rrm)
  return(tn95[,2:3])
}

ind_tyn <- function (series,fy=1700,ly=2100,rrm=11)
{
  tnn=data.frame(year=fy:ly,tnn=NA)
  for  (y in fy:ly)
  {
    if (length(series[series$year==y & series$month%in%c(1,2,12)& !is.na(series[,4]),4])>72)
    {
      tnn[tnn$year==y,'tnn']=min(series[series$year==y & !is.na(series[,4]),4])
    }
  }
  tnn[,3]=running_mean_2(tnn,rrm)
  colnames(tnn)=c('year','tnn','tnn.rm')
  return(tnn[,2:3])
}

ind_tyx <- function (series,fy=1700,ly=2100,rrm=11)
{
  tnn=data.frame(year=fy:ly,tnn=NA)
  for  (y in fy:ly)
  {
    if (length(series[series$year==y & series$month%in%c(6,7,8)& !is.na(series[,4]),4])>72)
    {
      tnn[tnn$year==y,'tnn']=max(series[series$year==y & !is.na(series[,4]),4])
    }
  }
  tnn[,3]=running_mean_2(tnn,rrm)
  colnames(tnn)=c('year','tnn','tnn.rm')
  return(tnn[,2:3])
}

### SOME MAPPING FUNCTIONS

map_analyzer <- function(table,TNoTX,index,operation,time_int, coord=c(-10,35,35,71),
                         table_off=data.frame(ser_id=NA,lat=NA,lon=NA,value=NA,signif=NA),
                         fact=5,finder=0,zoom=0,numbers=0,legend='topleft',printToFile=TRUE)
{
  text4main=index
  if (operation=='ori')
  {text4main='trends of ori.'}
  if (operation=='hom')
  {text4main='trends of 1 it.'}
  if (operation=='hom2')
  {text4main='trends of 2 it.'}
  if (operation=='diff')
  {text4main='diff. between trends of 2 it. and orig.'}
  if (operation=='diff1')
  {text4main='diff. between trends of 1 it. and orig.'}
  if (operation=='diff2')
  {text4main='diff. between trends of 2 it. and 1 it.'}
  colnames(table)[4]='value'
  colnames(table)[5]='signif'
  colnames(table_off)[4]='value'
  colnames(table_off)[5]='signif'
  nam=paste0('map_',TNoTX,'_',time_int,'_',index,'_',operation,'.png')
  width=7
  if (printToFile) {png(nam, width = width, height = width, res=600, units='in')}
  newmap <- getMap(resolution = "low")
  plot(newmap, xlim = coord[1:2], ylim = coord[3:4], asp = 1,
       main=paste0(TNoTX,', ', index,', ', text4main,', period: ', time_int))
  #mapCountryData(xlim = c(-10, 35), ylim = c(35, 71),oceanCol = 'lightblue')
  table_pos=table[table$value>0,]
  table_null=table[table$value==0,]
  table_neg=table[table$value<0,] 
  if (finder==0 & numbers==0)
  {
    points(table_null$lon, table_null$lat,cex=0.4,col='grey',lwd=3)
    points(table_neg$lon, table_neg$lat,cex=-fact*(table_neg$value),col='blue',lwd=ifelse(table_neg$signif>0.95 & !is.na(table_neg$signif),2.5,0.6))
    points(table_pos$lon, table_pos$lat,cex=(fact*table_pos$value),col='red',lwd=ifelse(table_pos$signif>0.95 & !is.na(table_pos$signif),2.5,0.6))
    points(table_off$lon, table_off$lat,cex=1.5,col='black',lwd=1.5,pch=4)
  }else
  {
    if (finder==1 & numbers==0)
    {
      points(table_null$lon, table_null$lat,pch=16,cex=fact*1.5,col='grey',lwd=3)
       points(table_pos$lon, table_pos$lat,pch=16,cex=fact*1.5,col=rgb(1,0,0,ifelse(table_pos$signif>0.95 & !is.na(table_pos$signif),1,0.3)),lwd=3)
       points(table_neg$lon, table_neg$lat,pch=16,cex=fact*1.5,col=rgb(0,0,1,ifelse(table_neg$signif>0.95 & !is.na(table_neg$signif),1,0.3)),lwd=3)
       points(table_off$lon, table_off$lat,cex=fact*1.5,col='black',lwd=1.5,pch=4)
    }
    else if (numbers==1 & finder==0)
    {
      text(table_off$lon, table_off$lat,cex=fact*0.8,col=('darkgrey'),lwd=1.5,labels=as.character(table_off$value))
      points(table_null$lon, table_null$lat,pch="0",cex=fact*1.5,col='darkgrey',lwd=3)
      text(c(table_neg$lon,NA), c(table_neg$lat,NA),labels=as.character(c(table_neg$value,NA)),
           cex=fact*1,col=rgb(0,0,1,ifelse(table_neg$signif>0.95 & !is.na(table_neg$signif),1,0.4)),lwd=3)
      text(c(table_pos$lon,NA), c(table_pos$lat,NA),labels=as.character(c(table_pos$value,NA)),
           cex=fact*1,col=rgb(1,0,0,ifelse(table_pos$signif>0.95 & !is.na(table_pos$signif),1,0.4)),lwd=3)
    }
  }
  if (substr(operation,1,4)=='diff')
  {
    legend('topleft',
           legend=c('neg. diff.','pos. diff.','no diff','size prop. to diff'),
           col=c('blue', 'red', 'grey','black'),
           pch=c(1,1,1,NA),
           #horiz=TRUE,
           inset=c(0.03,0.03),
           cex=fact)
  } else
  {
    if (finder==0)
    { 
      if (zoom==0)
      {
        legend(legend,
               legend=c('neg. trend','pos. trend','null trend','size prop. to trend'),
               col=c('blue', 'red', 'grey',NA),
               pch=c(16,16,16,NA),
               #horiz=TRUE,
               inset=c(0.03,0.03))
      }else
      {
        legend(legend,
               cex=1.3,
               legend=c('neg. trend','pos. trend','null trend','other stations','size prop. to trend'),
               col=c('blue', 'red', 'grey','black',NA),
               pch=c(1,1,1,4,NA),
               #horiz=TRUE,
               inset=c(0.03,0.03))
      }
    }else
    {
      if (zoom==0)
      {
        legend(legend,
               cex=1.2,
               legend=c('neg. trend','pos. trend','null trend','col. prop. to sign.'),
               col=c('blue', 'red', 'grey',NA),
               pch=c(16,16,16,NA),
               #horiz=TRUE,
               inset=c(0.03,0.03))
      }else
      {
        legend(legend,
               cex=1.2,
               legend=c('neg. trend','pos. trend','null trend','other stations','values in °C/dec'),
               col=c('blue', 'red', 'grey','darkgrey',NA),
               pch=c(16,16,16,16,NA),
               #cex=1.5,
               #horiz=TRUE,
               inset=c(0.03,0.03))
      }
    }
  }
  if (printToFile) {dev.off()}
}

map_trends <- function(table,TNoTX,index,operation,time_int,sup=9999,printToFile=TRUE)
{
  fct=4
  text4main=index
  if (operation=='ori')
  {text4main='trends of ori.'}
  if (operation=='hom')
  {text4main='trends of 1 it.'}
  if (operation=='hom2')
  {text4main='trends of 2 it.'}
  if (operation=='diff')
  {text4main='diff. between trends of 2 it. and orig.'
  fct=6}
  if (operation=='diff1')
  {text4main='diff. between trends of 1 it. and orig.'}
  if (operation=='diff2')
  {text4main='diff. between trends of 2 it. and 1 it.'}
  
  colnames(table)[4]='value'
  colnames(table)[5]='signif'
  
  
  nam=paste0('map_',TNoTX,'_',time_int,'_',index,'_',operation,'.png')
  width=7
  if(printToFile){  png(nam, width = width, height = width, units='in', res=600)}
  newmap <- getMap(resolution = "low")
  plot(newmap, xlim = c(-10, 35), ylim = c(35, 71), asp = 1,
       main=paste0(TNoTX,', ', index,', ', text4main,', period: ', time_int),cex.main=1.3)
  #mapCountryData(xlim = c(-10, 35), ylim = c(35, 71),oceanCol = 'lightblue')
  if (sup==9999)
  {
    table_pos=table[table$value>0,]
  }else
  {
    table_sup=table[table$value>sup,]
    table_pos=table[table$value>0 & table$value<=sup,]
  }
  
  table_null=table[table$value==0,]
  table_neg=table[table$value<0,] 
  points(table_null$lon, table_null$lat,cex=0.4,col='grey',lwd=3)
  points(table_pos$lon, table_pos$lat,cex=(fct*table_pos$value),col='red',lwd=ifelse(table_pos$signif>0.95,2.2,0.6))
  if (sup!=9999)
  {
    points(table_sup$lon, table_sup$lat,cex=(fct*table_sup$value),col='red4',lwd=ifelse(table_pos$signif>0.95,2.2,0.6))
  }
  points(table_neg$lon, table_neg$lat,cex=(-fct*table_neg$value),col='blue',lwd=ifelse(table_neg$signif>0.95,2.2,0.6))
  
  if (substr(operation,1,4)=='diff')
  {
    legend('topleft',
           legend=c('neg. diff.','pos. diff.','no diff','size prop. to diff'),
           col=c('blue', 'red', 'grey','black'),
           cex=1,
           pch=c(1,1,1,NA),
           #horiz=TRUE,
           inset=c(0.03,0.03))
  } else
  {
    if (sup==0)
    {
      legend('topleft',
             cex=1,
             legend=c('neg. trend','pos. trend','null trend','size prop. to trend'),
             col=c('blue', 'red', 'grey','black'),
             pch=c(1,1,1,NA),
             #horiz=TRUE,
             inset=c(0.03,0.03))
    }else
    {
      legend('topleft',
             legend=c('neg. trend',paste0('pos. trend, < ',sup,'°C/dec'),paste0('trend > ',sup,'°C/dec'),
                      'null trend','size prop. to trend'),
             cex=1,
             col=c('blue', 'red','red4', 'grey','black'),
             pch=c(1,1,1,1,NA),
             #horiz=TRUE,
             inset=c(0.03,0.03))
    }
  }
  
  if(printToFile){dev.off()}
}

map_trends_fill <- function(table,sig_table,TNoTX,index=NA,operation,time_int,qntl=50,
                            minmaxranges=c(-0.45,0.45),n=11,pltt=rev(brewer.pal(n, 'Spectral')),
                            coord=c(-15,41,28,70),printToFile=TRUE,width=7,height=5,magn=1,
                            homtest=0,season='ann',textit=0)                            
{ 
  ranges=data.frame(min=c(-9999,seq(minmaxranges[1]-0.0001,minmaxranges[2]+0.001,length.out=n-1)),
                    max=c(seq(minmaxranges[1]-0.0001,minmaxranges[2]+0.001,length.out=n-1),9999),
                    color=pltt,legend=NA)
  ranges$legend=paste(round(ranges$min,digits=2),'<',round(ranges$max,digits=2))
  ranges[1,'legend']=paste('min <',round(ranges$max[1],digits=2))
  ranges[nrow(ranges),'legend']=paste(round(ranges$min[nrow(ranges)],digits=2),'< max')
  
  table=table[,c('sta_id','lat','lon',paste0('q',sprintf('%02d',qntl)))]
  colnames(table)=c('sta_id','lat','lon','value')
  table$color=NA
  table$lwd=NA
  for (i in 1:nrow(table))
  {
    table[i,'color']=as.character(ranges[which(table[i,'value']>=ranges$min & table[i,'value']<=ranges$max),'color'])
    table[i,'lwd']=ifelse(sig_table[sig_table$sta_id==table[i,'sta_id'][1],paste0('q',sprintf('%02d',qntl))]>0.949,1.8,0.5)
  }
  table=table[order(abs(table$value)),]
  #nam=paste0('map_trends_',time_int,'_',TNoTX,'_',season,'_',operation,'_q',sprintf('%02d',qntl),'.png')
  if (operation=='ori'){opename='NewBlend'}
  if (operation=='hom'){opename='NewHomBlend'}
  if (qntl==0)
  {
    if (homtest==0)
    {
      title=paste0(TNoTX, ' ',season, ' diff of trends of 90p and 10p, ',time_int)
      nam=paste0('map_diff_trends_90p-10p_',time_int,'_',tolower(TNoTX),'_',season,'_',operation,'.png')
    }
    if (homtest!=0)
    {
      title=paste0(TNoTX, ' ',season, ' ', operation, ' of trends of ',homtest,'p (NewHomBlend-NewBlend), ',time_int)
      nam=paste0('map_diff_trends_hom-ori_',time_int,'_',tolower(TNoTX),'_',season,'_',operation,'_p',homtest,'.png')
    }
  }  else
    {title=paste0(TNoTX, ' ',season, ' ', opename, ' trends of ', qntl,'p ',time_int)
  nam=paste0('map_trends_',time_int,'_',tolower(TNoTX),'_',season,'_',operation,'_',sprintf('%02d',qntl),'p.png')}
  
  if(printToFile){png(nam,width=width,height=height,res=600,units='in')}
  par(mar=c(1,1,4,1))
  newmap <- getMap(resolution = "low")
  plot(newmap, xlim = c(coord[1],coord[2]), ylim = c(coord[3],coord[4]), asp = 1,
       main=title)
  
  points(table[table$lwd<1,'lon'],table[table$lwd<1,'lat'],lwd=table[table$lwd<1,'lwd'],
         cex=0.4+magn*abs(table[table$lwd<1,'value']*5),col='black',bg=table[table$lwd<1,'color'],pch=21)
  points(table[table$lwd>=1,'lon'],table[table$lwd>=1,'lat'],lwd=table[table$lwd>=1,'lwd'],
         cex=0.4+magn*abs(table[table$lwd>=1,'value']*5),col='black',bg=table[table$lwd>=1,'color'],pch=21)
  if (textit==1)
  {
    text(table$lon,table$lat,labels=table$sta_id)
  }
  
  
  legend('left',legend=c('[C/dec]',ranges$legend,'sign','non sign'),pch=c(NA,rep(21,n),1,1),col='black',
         pt.bg=c(NA,as.character(ranges$color),NA,NA),pt.lwd=c(NA,rep(1,n),1.8,0.5),
         pt.cex=c(NA,0.3+pmax(pmin(abs(round(ranges$min,digits=2)+round(ranges$max,digits=2))/2*5*magn,
                              c(seq(2,1.2,length.out=(n-1)/2),1,seq(1.2,2,length.out=(n-1)/2))),0.8),1,1),cex=0.8)
  if(printToFile){dev.off()}
  
}

monotonize <- function(vec,nq1) ### functions that takes a vector and makes it monotonous (not used anymore, replaced by the check of negative slopes)
{
  
  if (length(vec[!is.na(vec)])==nq1)
  {
    sen_slopes=matrix(NA,nq1,nq1)
    for (r in 1:nq1-1)
    {
      for (s in (r+1):nq1)
      {
        sen_slopes[r,s]=round((vec[s]-vec[r])/(s-r),digits=2)
      }
    }
    #fine=Sys.time()
    sgn=sign(median(sen_slopes, na.rm=TRUE))
    if (sgn==0)
    {
      if(abs(round(mean(sen_slopes, na.rm=TRUE),digits=3))>=0.002)
      {
        sgn=sign(round(mean(sen_slopes, na.rm=TRUE),digits=3))
      }
    }
    
    if (sgn!=0)
    {
      vec=vec*sgn   #enables to always work with increasing vectors
      r=1
      while (r<nq1)
      {
        if (vec[r+1]>=vec[r])
        {
          r=r+1;
        } else
        {
          s=r+1
          while(s<nq1 & vec[s+1]<=vec[s])
          {
            s=s+1
          }
          vec_d=vec[r:s]
          if (vec_d[length(vec_d)]>=0) #if the decreasing segment is above zero
          {
            l=r
            while(vec[l]>=vec_d[length(vec_d)] & l>=1)
            {
              l=l-1
              if (l==0)
              {
                break
              }
            }
            l=l+1
            vec[l:s]=vec[s]
            r=s
          }else if (vec_d[1]<=0) #if the decreasing segment is below zero
          {
            h=s
            while(vec[h]<=vec_d[1] & h<=nq1)
            {
              h=h+1
              if (h==nq1+1)
              {
                break
              }
            }
            h=h-1
            vec[r:h]=vec[r]
            r=h
          }else #if the decreasing segment is around zero
          {
            l=r
            h=s
            while(vec[l]>=0 & l>=1)
            {
              l=l-1
              if (l==0)
              {
                break
              }
            }
            l=l+1
            while(vec[h]<=0 & h<=nq1)
            {
              h=h+1
              if (h==nq1+1)
              {
                break
              }
            }
            h=h-1
            vec[l:h]=0
            r=h
          }
        }
        
      }
      vec=vec*sgn
      vec1=vec
      #durata=fine-inizio
      #print(durata)
    } else #there is no significant trend for the correction of the selected month
    {
      vec1=rep(round(mean(vec, na.rm=TRUE),digits=1),nq1)
    }
    
  }else
  {
    vec1=rep(NA,nq1)
  }
  return(vec1)
}

multi_hist_bm <- function(table,region,ele,operation,
                          bin,loc_leg='topright',dens_lwd=5,
                          printToFile=TRUE,outpath,maxcoeff=1.15,
                      inh=TRUE,suppl=NULL,bw=0,xlims=c(NA,NA))    ### creation of the histogram plots for homogenization
{
  setwd(outpath)
  if (inh){myclr=rbind(inh_row,myclr)}
  if (inh){clrs_legg=c(NA,clrs_legg)}
  if (operation=='pd05') 
  {
    xlab = 'PD05 [%]'
    opemain = 'PD05'
    opeleg = 'PD05'
    opefile = 'pd05'
  }
  if (operation=='rmse') 
  {
    xlab = 'RMSE'
    opemain = 'RMSE'
    opeleg = 'RMSE'
    opefile = 'rmse'
  }
  if (operation=='ind') 
  {
    xlab = paste0('ind on ',suppl)
    opemain = paste0('ind ',suppl)
    opeleg = paste0('ind ',suppl)
    opefile = paste0('ind_',suppl)
  }
  if (operation=='nonadj') 
  {
    xlab = paste0('non-adjusted data [%]')
    opemain = 'non-adjusted data'
    opeleg = ' perc. of non-adj. data'
    opefile = 'nonadj'
  }
  if (printToFile){png(paste0(operation,'_',suppl,'_',ele,'_',region,'.png'),height=600,width=1500)}
  table=table[,colSums(!is.na(table)) != 0]
  versions=c('qm','dap','homhom','splidhom','homdap','splidhomdap')
   myclr=data.frame(categories=versions,
                   clrs_hist=c(rgb(1,0,0.1,0.75),rgb(0,0.2,0.8,0.7),rgb(1,1,0,0.75),rgb(0,0.7,0.3,1),rgb(0.7,0,0.8,0.6),rgb(1,0.6,0,0.8)),
                   clrs_dens=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5),rgb(1,1,0,0.75),rgb(0,1,0,0.5),rgb(0.75,0,1,0.5),rgb(1,0.4,0,0.5)),
                   stringsAsFactors=FALSE)
   
   if(inh){versions=c('inh',versions)}
  ctgs=versions[versions%in%colnames(table)]
  table=table[,colnames(table)%in%ctgs]

  if(inh)
  {
    inh_row=c('inh',NA,rgb(0,0,0,0.5))
    myclr=rbind(inh_row,myclr)
  }
  ctgs_dens=ctgs
  ctgs=ctgs[which(ctgs!='inh')]
  clrs=myclr[myclr$categories%in%ctgs,'clrs_hist']
  nhom=length(ctgs)
  
  xmin=(floor(min(table[,ctgs],na.rm=TRUE)/bin)-1)*bin
  xmax=(ceiling(max(table[,ctgs],na.rm=TRUE)/bin)+1)*bin
  if(length(xlims[!is.na(xlims)])>=1)
  {
    x1=(floor(max(xmin,xlims[1],na.rm=TRUE)/bin))*bin
    x2=(ceiling(min(xmax,xlims[2],na.rm=TRUE)/bin))*bin
    xbreaks=unique(c(xmin,seq(x1-bin/2,x2+bin/2,bin),xmax))
    effmin=x1-bin
    plotlim=c(x1-bin/2,x2+bin/2)
  }  else
  {
    x1=xmin
    x2=xmax
    xbreaks=unique(seq(xmin-bin/2,xmax+bin/2,bin))
    effmin=xmin
    plotlim=c(-9999,9999)
  }
  if (operation== 'failure' | operation=='pd05')
  {
    x1=0
    x2=100
    xbreaks=unique(seq(x1-bin/2,x2+bin/2,bin))
    effmin=x1
  }
  cercamax=array(0)
  i=0
  for (ctg in ctgs)
  {
    i=i+1
    cercamax[i]=max(hist(table[,ctg],breaks=seq(xmin-bin/2,xmax+bin/2,bin),plot=FALSE)$counts)
  }
  multhist(table[,ctgs], col=clrs ,
           xlab=NA,
           breaks=xbreaks,cex.axis=1,
           ylim=c(0,maxcoeff*max(cercamax)))
  title(main=paste0(opemain,'  hist, ',ele,', dataset:',region),cex.main=1.5,
        xlab=xlab,cex.lab=1.5)
  
  if (operation=='failure'){abline(v=1+(nhom/2)-(nhom+1)*effmin/bin+(nhom+1)/bin*(80-bin/2),col='grey',lwd=2,lty=2)}
   for (ctg in ctgs_dens)
  {
    vec=table[,ctg]
    dens_lin_bm(vec=table[,ctg],clr=myclr[myclr$categories==ctg,'clrs_dens'],lwd=dens_lwd,
            nhom=nhom,xmin=effmin,bin=bin,plotlim=plotlim,bw=bw)
  }
  
  legg=data.frame(ctgs=ctgs_dens,legg=NA)
  for (ctg in legg$ctgs)
  {
    vec=table[,ctg][!is.na(table[,ctg])]
    legg[legg$ctg==ctg,'legg']=paste(sprintf("%s",ctg),':',sprintf("%.2f",round(mean(vec),digits=2)))
  }
  clrs_legg=myclr[myclr$categories%in%ctgs_dens,'clrs_dens']
  
 
  legend(loc_leg,
         legend=c(paste0('method : mean of ', opeleg),legg[legg$ctg%in%ctgs_dens,'legg']),
         col=c(NA,clrs_legg),cex=1.5,
         lwd=5)
  if(printToFile)dev.off()
}

one2twelve <- function(n) ### brings all integers to the [1,12] interval , useful in particular monthly based analyses
{
  if (!(n%in%(1:12))){n=n%%12}
  if (n==0) {n=12}
  return(n)
}

overlapize_ref <- function(ref,ser,thre_y=5) #calculates overlap between basis/donat and refer
{
  overlap=merge(ser,ref)
  overlap=overlap[order(overlap$day),]
  overlap=overlap[order(overlap$month),]
  overlap=overlap[order(overlap$year),]
  if (colnames(overlap)[4]=='basis')
  {
    append=data.frame(year=rep(9999,12),month=seq(1,12,1),day=rep(45,12),basis=NA,refer=NA)
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

### PLOT FUNCTIONS

plot_adj_from_matrix <- function(list_mtrx,month,ser_id,break_year,add,mdn_yon=0,minmin=NA,maxmax=NA) ### plot adjustments sequences from matrix
{
  m=month
  corrs_list=list()
  for (i in 1:length(list_mtrx))
  {
    corrs_list[[i]]=list_mtrx[[i]][,m]
  }
  if (is.na(minmin))
  {
    mins=sapply(corrs_list,min)
    mins=mins[!is.na(mins)]
    minmin=min(mins)
  }
  if (is.na(maxmax))
  {
    maxs=sapply(corrs_list,max)
    maxs=maxs[!is.na(maxs)]
    maxmax=max(maxs)
  }
  clr=rep(c('blue','red','green','magenta','orange','cyan','maroon'),3)
  pch_list=rep(c(16,17,18),each=7)
  nam=paste0('adjustments_',m,'_',ser_id,'_',break_year,'_',add,'.png')
  png(nam, width = 800, height = 500)
  par(mar=c(4, 4.1, 4, 8),xpd=FALSE)
  plot(seq(5,95,5),corrs_list[[1]],type='l',col=clr[1],
       ylim=c(minmin,maxmax),
       xlab='quantiles',ylab='adjustmens[C]',
       main=paste0('Adj. est. ',add,', month ' ,m, ', ser id ',ser_id, ', break ', break_year))
  points(seq(5,95,5),corrs_list[[1]],col=clr[1],pch=pch_list[1])
  abline(h=0,col='grey',lwd=3,lty=3)
  abline(h=c(-3,-2.5,-2,-1.5,-1,-0.5,0.5,1),col='grey',lwd=1,lty=3)
  for (i in 2:length(corrs_list))
  {
    lines(seq(5,95,5),corrs_list[[i]],col=clr[i])
    points(seq(5,95,5),corrs_list[[i]],col=clr[i],pch=pch_list[i])
  }
  
  if(mdn_yon==1)
  {
    mdn=rep(NA,19)
    aux=rep(NA,length(corrs_list))
    for (q in 1:19)
    {
      for (i in 1:length(corrs_list))
      {
        aux[i]=corrs_list[[i]][q]
      }
      aux=unlist(aux)
      aux=aux[!is.na(aux)]
      mdn[q]=wise_median(aux)
    }
    lines(seq(5,95,5),mdn,col='black',lwd=5)
  }
  
  
  ref_ids=sapply(names(list_mtrx),extract_ser_id)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  if (mdn_yon==1)
  {
    legend('right',
           legend=c(ref_ids,'median'),
           col=c(clr[1:length(ref_ids)],'black'),
           pch=c(pch_list[1:length(ref_ids)],NA),
           lwd=c(rep(1,length(ref_ids)),5),
           lty=c(rep(1,length(ref_ids)),1),
           horiz=FALSE,
           inset=c(0.03,0.03))
  }else if (mdn_yon==0)
  {
    legend('right',
           legend=c(ref_ids),
           col=c(clr[1:length(ref_ids)]),
           pch=c(pch_list[1:length(ref_ids)]),
           lwd=c(rep(1,length(ref_ids))),
           lty=c(rep(1,length(ref_ids))),
           horiz=FALSE,
           inset=c(0.03,0.03))
  }
  dev.off()
}

plot_ind <- function(table,nox,ind,time='ann',ser,namecountry='Elsewhere',homplot=1,chp1,chp2=NA,
                     lab_size=1,plot_title=TRUE,printToFile=TRUE)
{
  if (nox=='n'){TNoTX='TN'}
  if (nox=='x'){TNoTX='TX'}
  if (nox=='g'){TNoTX='TG'}
  if (ind=='mean'){index='Mean'}
  if (ind=='sd'){index='Summer Days'}
  if (ind=='fd'){index='Frost Days'}
  if (ind=='tnn'){index='TNn'}
  if (ind=='tnx'){index='TNx'}
  if (ind=='txn'){index='TXn'}
  if (ind=='txx'){index='TXx'}
  if (ind=='tn05'){index='TN05'}
  if (ind=='tn95'){index='TN95'}
  if (ind=='tx05'){index='TX05'}
  if (ind=='tx95'){index='TX95'}
  
  print(paste('Plotting',ind,time))
  if (homplot[1]!=0)
  {
    ifelse(is.numeric(ser),
           nam<-paste0('plot_index_',TNoTX,'_',ser+1000000,'_',ind,'_',time,'_2it.png'),
           nam<-paste0('plot_index_',TNoTX,'_',as.character(ser),'_',ind,'_',time,'_hom.png'))
    maintext=paste0(index,' ',time, '(2it) t',nox,' ',ser, ' ', namecountry)
  } else
  {
    ifelse(is.numeric(ser),
           nam<-paste0('plot_index_',TNoTX,'_',ser+1000000,'_',ind,'_',time,'_ori.png'),
           nam<-paste0('plot_index_',TNoTX,'_',as.character(ser),'_',ind,'_',time,'_ori.png'))
    maintext=paste0(index,' ',time, 'ori t',nox,' ',ser, ' ', namecountry)
  }
  fct=1
  if (printToFile==TRUE)
  {
    width=7
    fct=1.5
    png(nam, width = width, height = width/10*6, units='in',res=600)
  }

  if (is.na(chp2[1])) 
  {
    table$hom2=NA
    table$hom2.rm=NA
  }
  par(mar=c(4.2, 2.6, 2.5, 1),xpd=FALSE)
  plot(table$year,table$ori,type='p',
       xlab='', ylab='',
       xaxt='n',yaxt='n',
       ylim=c(min(table$ori[!is.na(table$ori)],table$hom[!is.na(table$hom)],table$hom2[!is.na(table$hom2)]),
              max(table$ori[!is.na(table$ori)],table$hom[!is.na(table$hom)],table$hom2[!is.na(table$hom2)])),
       lwd=2/fct,cex=1/fct)
  if (plot_title==TRUE) {title(main=maintext,cex.main=lab_size)}
  title(xlab='Years',
        ylab=paste0(index,' t',nox,' [C]'),
        cex.lab=lab_size, mgp =c(1.3, 1, 1))
  axis(side=1, cex.axis=lab_size, mgp=c(1.8,0.5,0))
  axis(side=2, cex.axis=lab_size, mgp=c(1.8,0.5,0))
  lines(table$year,table$ori.rm,col='black',lwd=2/fct)
  if (homplot[1]!=0)
  {
    points(table$year,table$hom,col='blue',lwd=2/fct,cex=1/fct)
    lines(table$year,table$hom.rm,col='blue',lwd=2/fct)
    points(table$year,table$hom2,col='red',lwd=2/fct,cex=1/fct)
    lines(table$year,table$hom2.rm,col='red',lwd=2/fct)
  }
  abline(v=seq(1750,2020,5),lwd=1/fct,col='grey',lty=3)
  if (ind!='sd' & ind!='fd')
  {
    abline(h=seq(-50,50,1),lwd=1/fct,col='grey',lty=3)
  }
  abline(h=seq(-50,450,10),lwd=1/fct,col='grey',lty=3)
  
  abline(v=chp1,col='blue',lwd=3/fct)
  if (homplot[1]!=0)
  {
    abline(v=chp2,col='red',lwd=3/fct)
  }
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  if (homplot[1]!=0)
  {
    if (!is.na(chp2[1]))
    {
    legend('bottom',
           cex=lab_size*0.8,
           legend=c('ori','iter.1','iter.2','gaussian weighted running mean'),
           col=c('black',  'blue','red',NA),
           pch=c(1,1,1,NA),
           lwd=c(2,2,2,NA),
           horiz=TRUE,
           inset=c(1,0.03,1),
           x.intersp = 0)
    }else
    {
      legend('bottom',
             legend=c('ori','iter.1','iter. 2 not needed','gaussian weighted running mean'),
             col=c('black',  'blue','red',NA),
             pch=c(1,1,1,NA),
             lwd=c(2,2,2,NA),
             horiz=TRUE,
             inset=c(1,0.03,1),
             x.intersp = 0)
    }
  }else
  {
    legend('bottom',
           legend=c('original series','gaussian weighted running mean'),
           col=c('black',NA),
           pch=c(1,NA),
           lwd=c(2,NA),
           horiz=TRUE,
           inset=c(0.03,0.03))
  }
  if (printToFile==TRUE)
  {
    dev.off()
  }
}

plot_ind_bh <- function(table, nox, ind, time='ann', sta_id, blend_id_series,  namecountry='Elsewhere', printToFile=TRUE,homplot=1,
                        lab_size=1,season='ann')
{

  index=ind
  if (nox=='n'){TNoTX='TN'}
  if (nox=='x'){TNoTX='TX'}
  if (nox=='g'){TNoTX='TG'}
  if (ind=='annmean'){index='Annual Mean'}
  if (ind=='sd'){index='Summer Days'}
  if (ind=='fd'){index='Frost Days'}
  if (ind=='tnn'){index='TNn'}
  if (ind=='tnx'){index='TNx'}
  if (ind=='txn'){index='TXn'}
  if (ind=='txx'){index='TXx'}
  if (ind=='tn05'){index='TN05'}
  if (ind=='tn95'){index='TN95'}
  if (ind=='tx05'){index='TX05'}
  if (ind=='tx95'){index='TX95'}
  
  
  value=sta_id+2000000
  if (homplot[1]!=0)
  {
    nam=paste0('plot_index_bh_',TNoTX,'_',value,'_',ind,'_',season,'_hom.png')
  } else
  {
    nam=paste0('plot_index_bh_',TNoTX,'_',value,'_',ind,'_',season,'_ori.png')
  }
  
  nox=tolower(TNoTX)
  fct=1
  if(printToFile){width=7
  fct=1.5
    (png(nam, width = 7, height = 7/10*6,units='in',res=600))}
  par(mar=c(4.5, 2.6, 2.5, 1),xpd=FALSE)
  #par(mar=c(8, 4.1, 4, 2),xpd=FALSE)
  plot(table$year,table$ori,type='p',
       xlab='', ylab='',
       xaxt='n',yaxt='n',
       lwd=2/fct,cex=1/fct,
       ylim=c(floor(min(min(table$ori,na.rm=TRUE),min(table$hom,na.rm=TRUE))),
              ceiling(max(max(table$ori,na.rm=TRUE),max(table$hom,na.rm=TRUE)))))
  title(main=paste0(index,' ',TNoTX,' ',
                    #sta_id, ' ', 
                    namecountry,' ',season),cex.main=lab_size/fct)
  title(xlab='Years',
        ylab=paste(ind,TNoTX,'[°C]'),
        cex.lab=lab_size, mgp =c(1.3, 1, 1))
  #points(table$year,table$value,col='darkgrey',lwd=2/fct,cex=1/fct,pch=2)
  #lines(table$year,table$value.rrm,col='darkgrey',lwd=2/fct)
  lines(table$year,table$ori.rm,col='black',lwd=2/fct)
  axis(side=1, cex.axis=lab_size, mgp=c(1.8,0.5,0))
  axis(side=2, cex.axis=lab_size, mgp=c(1.8,0.5,0))
  abline(v=seq(1750,2020,5),lwd=1,col='grey',lty=3)
  abline(h=seq(-50,50,1),lwd=1,col='grey',lty=3)
  blend_id_series=blend_id_series[blend_id_series$year>=min(table$year) & blend_id_series$year<=max(table$year),]
  blend_id_series$yvalue=floor(min(min(table$ori,na.rm=TRUE),min(table$hom,na.rm=TRUE)))-0.05
  points((blend_id_series$year+(blend_id_series$month)/12+(blend_id_series$day)/365),blend_id_series$yvalue,
         col=blend_id_series$col,cex=1/fct,pch="|")
  blend_id_syn=blend_id_series[blend_id_series$blend_id>=900000,]
  # if(length(unique(blend_id_series$col))<=3)
  # {
  # points((blend_id_syn$year+(blend_id_syn$month)/12+(blend_id_syn$day)/365),blend_id_syn$yvalue,
  #        col='green',cex=1/fct,pch="|",lwd=1/fct*1.2)
  # }else
  # {
  points((blend_id_syn$year+(blend_id_syn$month)/12+(blend_id_syn$day)/365),blend_id_syn$yvalue,
         col=blend_id_syn$col,cex=1/fct,pch="|",lwd=1/fct*1.2)
  # }
  if (homplot[1]!=0)
  {
    points(table$year,table$hom,col='red',lwd=2/fct,cex=1/fct)
    lines(table$year,table$hom.rm,col='red',lwd=2/fct)
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend('bottomleft',
         legend=c('NewBlend','NewBlendHom'),
         col=c('black',  'red'),
         pch=c(1,1),
         lwd=c(2,2),
         horiz=TRUE,
         inset=c(0.03,0.03),
         cex=0.8*lab_size)
  colsb=unique(blend_id_series$col)
  blends=unique(blend_id_series$blend_id)
  #colsb[which(blends>900000)]='black'
  blends[blends>900000]=paste('WMO',sprintf('%05i',blends[blends>900000]-900000))
  legend('bottomright',
         legend=blends,
         col=colsb,
         pch=rep("|",length(blends)),
         cex=0.8*lab_size,
         pt.cex=1.5,
         pt.lwd=50,
         horiz=TRUE,
         inset=c(0.03,0.03))
  if(printToFile)dev.off()
  
}
  
plot_ind_bm <- function(table,nox,ind,time='ann',ser,namecountry='Elsewhere',num_plots=1,homplot=1,chp1=2000,chp2=NA,
                         lab_size=1,plot_title=TRUE,printToFile=TRUE,width=7,height=width/7*6,
                        leg_size=lab_size,ttl_size=lab_size)
{
  index=ind
  if (nox=='n'){TNoTX='TN'}
  if (nox=='x'){TNoTX='TX'}
  if (ind=='mean'){index='Mean'}
  if (ind=='qntl'){index='Quantile'}
  if (ind=='sd'){index='Summer Days'}
  if (ind=='fd'){index='Frost Days'}
  if (ind=='tnn'){index='TNn'}
  if (ind=='tnx'){index='TNx'}
  if (ind=='txn'){index='TXn'}
  if (ind=='txx'){index='TXx'}
  if (ind=='tn05'){index='TN05'}
  if (ind=='tn95'){index='TN95'}
  if (ind=='tx05'){index='TX05'}
  if (ind=='tx95'){index='TX95'}
  fy=table[1,1]
#  categories=colnames(table)[c(FALSE,TRUE)]
  categories=c('ori','inh','qm','dap','homhom','splidhom','homdap','splidhomdap')
  
    for (cat in categories)
  {
    if(length(table[!is.na(table[,cat]),cat])>=24)
    {
    assign(x=paste0('trend_',cat),value=round(sens.slope(table[!is.na(table[,cat]),cat])$estimates*10,digits=2))
    }
    else
    {
      assign(x=paste0('trend_',cat),NA)
    }
  }
  
  print(paste('Plotting',ind,time))
  if (homplot[1]!=0)
  {
    ifelse(is.numeric(ser),
           nam<-paste0('plot_index_',TNoTX,'_',ser+1000000,'_',ind,'_',time,'.png'),
           nam<-paste0('plot_index_',TNoTX,'_',as.character(ser),'_',ind,'_',time,'_hom.png'))
    maintext=paste0('T',toupper(nox),' ',index,' ',time, ' ',ser, ' ', namecountry)
  } else
  {
    ifelse(is.numeric(ser),
           nam<-paste0('plot_index_',TNoTX,'_',ser+1000000,'_',ind,'_',time,'_ori.png'),
           nam<-paste0('plot_index_',TNoTX,'_',as.character(ser),'_',ind,'_',time,'_ori.png'))
  }
  if (printToFile==TRUE)
  {
    png(nam, width = width, height = height, res=600, units='in')
  }
  
  layout(matrix(c(1,2),2,1), widths = 1, heights = c(3,2))
  par(mar = c(lab_size, 4.1, 4.1, 2.1))
  #par(mar=c(8, 4.1, 4, 2),xpd=FALSE)
  plot(table$year,table$inh,type='p',pch=6,lwd=4,
       xlab='', ylab='',
       xaxt='n',yaxt='n',
       ylim=c(min(table[,2:ncol(table)],na.rm=TRUE),
              max(table[,2:ncol(table)],na.rm=TRUE)),
       panel.first=c(lines(table$year,table$ori.rm,col='grey',lwd=8),
                     points(table$year,table$ori,col='grey',lwd=4,pch=2)))
  if (plot_title==TRUE) {title(main=maintext,cex.main=ttl_size)}
  title(ylab=paste0(index,' t',nox,' [°C]'),
        cex.lab=lab_size)
  axis(side=1, cex.axis=lab_size)
  axis(side=2, cex.axis=lab_size)
  lines(table$year,table$inh.rm,col='black',lwd=6)
  table=index_benchmark(table)
  
  if (length(chp1)==0){chp1=1979}
  if (is.na(max(chp1))){chp1=1979}
  if (min(chp1)<1970){chp1=1979}
  if (max(chp1)>2000){chp1=1979}
    
  if (homplot[1]!=0)
  {
#    table[,c('ori','ori.rm','inh','inh.rm')]=NULL
    numhom=length(categories[-which(categories=='ori'|categories=='inh')])
    clrs=clrsFun(numhom)
    if(any(categories=='qm')){clrs[which(categories=='qm')-2]='red'}
    if(any(categories=='dap')){clrs[which(categories=='dap')-2]='blue'}
    if(any(categories=='homhom')){clrs[which(categories=='homhom')-2]='gold1'}
    if(any(categories=='splidhom')){clrs[which(categories=='splidhom')-2]='green'}
    if(any(categories=='homdap')){clrs[which(categories=='homdap')-2]='purple'}
    if(any(categories=='splidhomdap')){clrs[which(categories=='splidhomdap')-2]='darkorange1'}
    if (numhom>=1)
    {
      for (j in c(1,3,4,5,6,2))
     
     {
        print(colnames(table)[2*j+4])
        lines(table$year,table[,(2*j)+1+4],col=clrs[j],lwd=3)
      }
      lines(table$year,table$inh.rm,col='black',lwd=4)
      points(table$year,table$ori,col='grey',lwd=4,pch=2)
      for (j in c(1,3,4,5,6,2))
        
      {
        points(table[table$year<max(chp1),'year'],
               table[table$year<max(chp1),(2*j+4)],col=clrs[j],lwd=2)
      }
    }
  } 
  # points(table$year,table$inh,col='black',lwd=2)
  #

  abline(v=seq(1800,2020,5),lwd=1,col='grey',lty=3)
  if (ind!='sd' & ind!='fd')
  {
    abline(h=seq(-50,50,1),lwd=1,col='grey',lty=3)
  }
  abline(h=seq(-50,450,10),lwd=1,col='grey',lty=3)
  if (numhom>=1)
  {  
    if (homplot[1]!=0)
  {
    abline(v=chp1,col='black',lwd=3)
  }
    par(mar = c(4.1, 4.1, lab_size-0.5, 2.1))
    #par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    
    plot(0, 0,
         #ylim=c(min(table[,((numhom*2+2):ncol(table))],-1,na.rm = TRUE),
         #        max(table[,((numhom*2+2):ncol(table))],3,na.rm = TRUE)),
         ylim=c(-1,3),
         xlim=c(table$year[1],table$year[nrow(table)]),
         xlab=NA,ylab=NA,xaxt='n',yaxt='n',
         axes=FALSE)
    title(ylab=paste0('hom index [°C]'),
          cex.lab=lab_size)
    axis(side=2, cex.axis=lab_size)
    axis(side=3, labels=FALSE)
  
    
  text(table$year[1]+5,2.25,'overcorrected',col='darkgrey')
  text(table$year[1]+5,-0.25,'wrongcorrected',col='darkgrey')
  lines(c(1750,max(chp1)),c(1,1),lwd=3,col='darkgreen',lty=1)
  lines(c(1750,max(chp1)),c(0,0),lwd=3,col='yellow',lty=1)
  lines(c(1750,max(chp1)),c(2,2),lwd=3,col='yellow',lty=1)
  lines(c(1750,max(chp1)),c(-1,-1),lwd=3,col='red',lty=1)
  lines(c(1750,max(chp1)),c(3,3),lwd=3,col='red',lty=1)
  for (ctg in categories[-which(categories=='ori'|categories=='inh')][c(1,3,4,5,6,2)])
  {
    j=which(categories==ctg)-2
    points(table[table$year<max(chp1),'year'],table[table$year<max(chp1),paste0(ctg,'.ind')],
           col=clrs[j],lwd=3,pch=2)
    lines(table[table$year<max(chp1),'year'],table[table$year<max(chp1),paste0(ctg,'.ind.rm')],
          col=clrs[j],lwd=3)
  }
  rect(1750, 2, max(chp1), 1000, border = NA, col = rgb(1,0,0.1,0.1))
  rect(1750, 1.5, max(chp1), 2, border = NA, col = rgb(1,1,0,0.1))
  rect(1750, 0.5, max(chp1), 1.5, border = NA, col = rgb(0,1,0.2,0.1))
  rect(1750, 0, max(chp1), 0.5, border = NA, col = rgb(1,1,0,0.1))
  rect(1750, -1000, max(chp1), 0, border = NA, col = rgb(1,0,0.1,0.1))
  abline(v=seq(1800,max(chp1),5),lwd=1,col='grey',lty=3)
 
  for (h in seq(-50,450,0.5)[-which(seq(-50,450,0.5)==1)])
  {
    lines(c(1750,max(chp1)),c(h,h),lwd=1,col='grey',lty=3)
  }
  
  #abline(v=chp1,col='black',lwd=1.5)
  legg=array(NA,length(categories))
  nams=categories
  nams[1:2]=c('benchm.','test')
  j=0
  for (i in 1:length(categories))
  {
    ctg=categories[i]
    if (length(table[!is.na(table[ctg]),ctg])>1)
    {
      if (ind!='sd' & ind!='fd')
      {
        j=j+1
        legg[j]<-paste0(nams[i],' (',sprintf('%2.2f',eval(parse(text=paste0('trend_',ctg)))),'°C/dec)')
      }else
      {
        legg[j]<-paste0(nams[i],' (',sprintf('%2.2f',eval(parse(text=paste0('trend_',ctg)))),' days/dec)')
      }
    }
  }
  
  # aux=legg[2]
  # legg[2]=legg[1]
  # legg[1]=aux
  
  numhom=length(legg[which(!is.na(legg[3:8]))])
  
  legend('right',
         legend=legg[which(!is.na(legg))],
         col=c('grey','black',clrs[which(!is.na(legg[3:8]))]),
         pch=c(2,2,rep(1,numhom)),
         pt.lwd=c(8,6,rep(2,numhom)),
         lwd=c(8,6,rep(3,numhom)),bty='n',
         #ncol=numhom%/%4+1,
         cex=leg_size)
  }
  if (printToFile==TRUE)
  {
    dev.off()
  }
}


qntls_3m <- function(vec,nb1,nq1,ampl_bin1,thre_y1=5) #OK.Takes the bndrs of the intervals as input to calculate  (with the median) the values of the qntls
{
  if (length(vec[!is.na(vec)])>=thre_y1*80)
  {
    ll=length(vec[!is.na(vec)])/(100/ampl_bin1)                    #number of data in each interval
    #to calculate the values of the quantile you need to know which are the values of the boundaries of each interval
    #card_bnd1 is the array of the indices of the bndrs needed, in this case 2.5, 7.5, 12.5, ... ,97.5 percentiles
    card_bnd1=round(seq(0.5,nb1-0.5)*ll)
    vec1_table=data.frame(value=sort(vec[!is.na(vec)]),cat=NA)
    for (q in (1:nq1)) #for each quantile
    {
      vec1_table[((card_bnd1[q]+1):card_bnd1[q+1]),'cat']=q #assign every data to an interval, bounded by the boundaries in ca
    }
    val_qnt1=tapply(vec1_table$value,vec1_table$cat,wise_median)
  }else
  {
    val_qnt1=rep(NA,nq1)
  }
  val_qnt1=round(val_qnt1,digits = 1)
  return(val_qnt1)
}

qntl_corrs_ref <- function(overlap)
{
  control=check_month_avail(overlap,thre_y)
  if (control==1)
  {
    if (colnames(overlap)[4]=='basis')
    {
      qnt=as.data.frame.list(tapply(overlap$basis,overlap$month,function (x) quantilize(x,nb,nq,ampl_bin,thre_y)))
    }else if (colnames(overlap)[4]=='donat')
    {
      qnt=as.data.frame.list(tapply(overlap$donat,overlap$month,function (x) quantilize(x,nb,nq,ampl_bin,thre_y)))
    }
    
    if (length(colnames(qnt))==12)
    {
      colnames(qnt)=c('01','02','03','04','05','06','07','08','09','10','11','12')
    }else
    {
      printf('Error: number of month columns of quantile matrix is not 12')
      break
    }
    row.names(qnt)=seq(fb,lb,ampl_bin)
    r_qnt=as.data.frame.list(tapply(overlap$refer,overlap$month,function (x) quantilize(x,nb,nq,ampl_bin,thre_y)))
    colnames(r_qnt)=c('01','02','03','04','05','06','07','08','09','10','11','12')
    row.names(r_qnt)=seq(fb,lb,ampl_bin)
    corr=round(qnt-r_qnt, digits=1)
  }  else
  {
    corr=data.frame('01'=rep(NA,nb),'02'=NA,'03'=NA,'04'=NA,'05'=NA,'06'=NA,'07'=NA,'08'=NA,'09'=NA,'10'=NA,'11'=NA,'12'=NA)
    colnames(corr)=c('01','02','03','04','05','06','07','08','09','10','11','12')
  }
  return(corr)
}

quality_check <- function(QC1,couples1,TNoTX1,series1,ser_id1)
{
  if(TNoTX1=='TN')
  {
    ser_id2=couples[couples$tn_id==ser_id1,'tx_id']
    subQC=QC1[QC1$ser_id==ser_id2 & QC$element=='TMIN',]
    subQCI=QC1[QC1$ser_id==ser_id2 & QC$element=='TMIN' & QC$flag=='I',]
  } else if (TNoTX1=='TX')
  {
    ser_id2=ser_id1
    subQC=QC1[QC1$ser_id==ser_id2 & QC1$element=='TMAX',]
    subQCI=QC1[QC1$ser_id==ser_id2 & QC$element=='TMAX' & QC$flag=='I',]
  }
  series1$idate=series1$year*10000+series1$month*100+series1$day
  subQC$idate=subQC$year*10000+subQC$month*100+subQC$day
  if (ser_id1==75025)
  {
    subQC=subQC[subQC$year!=1940 | subQC$month!=1,]
  }
  if (dim(subQC)[1]>0 & dim(series1[which(series1$idate%in%subQC$idate),])[1]>0)
  {
    series2=series1[-which(series1$idate%in%subQC$idate),]
  }else
  {
    series2=series1
  }
  
  years=unique(subQCI$year)
  
  for (y in years)
  {
    if (dim(subQCI[subQCI$year==y,])[1]>dim(series1[series1$year==y,])[1]*0.25 & dim(series1[series1$year==y,])[1]>0)
    {
      series2=series2[(which(series2$year!=y)),]
      print (paste0('Flagged I remotion of series ',ser_id1,'. Year:',y))
    }
  }
  
  ### REPEATED VALUES
  mem=colnames(series2)[4]
    colnames(series2)[4]='series'
  i=0
  while (i<(nrow(series2)))
  {
    i=i+1
    checkd=abs(series2[i+1,'series']-series2[i,'series'])
    d=2
    while (checkd<0.1 & d<((nrow(series2))-i))
    {
      checkd=abs(series2[i+d,'series']-series2[i,'series'])
      d=d+1
    }
    d=d-2
    if (d>=5) 
    {
      print(paste0(ser_id1,': ',d,' repeated values'))
      series2[i:(i+d),'series']=-999.9
      i=i+d
    }
  }
  series2=series2[series2$series>(-50),]
  colnames(series2)[4]=mem
  series2$idate=NULL
  
  return(series2)
  
}

quantilize <- function(vec,nb1,nq1,ampl_bin1,thre_y1) ### simple function for the creation of a quantile sequence
{
  if (length(vec)>=thre_y1*29)
  {
    ll=length(vec[!is.na(vec)])/(100/ampl_bin1)  #number of data in each bin
    card_bnd1=round(seq(0.5,nb1-0.5)*ll)       #position in the sorted data corresponding to the boundaries of the bins
    vec1_table=data.frame(value=sort(vec[!is.na(vec)]),cat=NA)
    for (q in 1:nq1)
    {
      vec1_table[((card_bnd1[q]+1):card_bnd1[q+1]),'cat']=q
    }
    val_qnt1=tapply(vec1_table$value,vec1_table$cat,wise_median)
  }else
  {
    val_qnt1=rep(NA,nq1)
  }  
  return(val_qnt1)
}

read.breaks <- function(ser_id,prefix,iter,breaks)
{
  if (prefix=='n') ELE='TN'
  if (prefix=='x') ELE='TX'
  if (prefix=='g') ELE='TG'
  value=sprintf("%06d",ser_id)
  if (iter==1)
  {
    bd_name=paste0(prefix,ELE,'_serid_',value,'_ori.txt')
    range=breaks[bd_name][[1]]$range
  }
  if (iter==2)
  {
    bd_name=paste0(prefix,ELE,'_serid_',value,'_parthom.txt') ### is there a partially homogenized version of this series?
    range=breaks[bd_name][[1]]$range
    if (length(range)==0)
    {
      bd_name=paste0(prefix,ELE,'_serid_',value,'_ori.txt') ### is there a ori version of this series?
      range=breaks[bd_name][[1]]$range
      if (length(range)==0)
      {
        print(paste('Error: the series',ref_id,'is neither hom or ori.'))
      }
    }
  }   
  if (length(range)==0)
  {
   print(paste0('Break Detection has failed on series ',ser_id))
   chp=list()
   chp$year=NULL
   chp$month=NULL
   chp$fail=TRUE
   return(chp)
  }
   
  chp_y=breaks[bd_name][[1]]$all2
  chp_m=breaks[bd_name][[1]]$m
  chp_fail=FALSE
  chp=list(chp_y,chp_m,chp_fail)
  names(chp)=c('year','month','fail')  
  return(chp)
}

read.breaks.bm <- function(ser_id,nox,iter,breaks,break_directory='/usr/people/squintu/Documents/ECAD/Benchmarking/Test_Series',
                           base=3,code='ser')   ### read breaks for benchmarking
{
  setwd(break_directory)
  if (nox=='n') TNoTX='TN'
  if (nox=='x') TNoTX='TX'
  if (nox=='g') TNoTX='TG'
  value=base*1000000+ser_id
  yuri_name=paste0(TNoTX,'_',value,'.',code)
  chp_y=breaks[yuri_name][[1]]$year
  chp_m=breaks[yuri_name][[1]]$month
  chp=list(chp_y,chp_m)
  names(chp)=c('year','month')
  return(chp)
}

read.chp.y <- function (filename,TNoTX,bd_num)
{
  ser=extract_ser_id(filename)
  nox=ifelse(TNoTX=='TN','n','x')
  setwd(workfolder)
  if (TNoTX=='TN')
  {
    load("n_breaks.RData")
    breaks1=breaks
    if (bd_num==2)
    {
      load("n_breaks_2.RData")
      breaks2=breaks
    }
    nox='n'
  }  else if (TNoTX=='TX')
  {
    load("x_breaks.RData")
    breaks1=breaks
    if (bd_num==2)
    {
      load("x_breaks_2.RData")
      breaks2=breaks
    }
    nox='x'
  }
  setwd(inpath)
  yuri_value=100000+ser
  yuri_name=paste0(nox,TNoTX,'_SOUID',yuri_value,'.txt')
  chp1=breaks1[yuri_name][[1]]$all2
  chp=chp1
  chp_2=NULL
  if (bd_num==2)
  {
    yuri_value=1000000+ser
    yuri_name=paste0(nox,TNoTX,'_',yuri_value,'.hom')
    if(length(breaks2[yuri_name][[1]]$range)==0)
    {
      yuri_name=paste0(nox,TNoTX,'_',yuri_value,'.ori')
    }
    chp2=breaks2[yuri_name][[1]]$all2
    chp=list(chp1=chp1,chp2=chp2)
  }
  
}

read.csv.4.ind <- function (filename,TNoTX,bd_num=2,inpath, workfolder= "/usr/people/squintu/Documents/ECAD/Homog/WorkFolder")
{
  series=read.table(filename)
  if (ncol(series)==6)
  {
    series_raw=read.csv(filename,header=FALSE)
    colnames(series_raw)=c('ser_id','date','value','qc0','qc1','qc2')
    series_raw[series_raw == "-9999"] <- NA
    series=data.frame(year=NA,month=NA,day=NA,value=series_raw$value/10)
    series$year=floor(series_raw$date/10000)
    series$month=floor((series_raw$date-floor(series_raw$date/10000)*10000)/100)
    series$day=series_raw$date-series$year*10000-series$month*100
  }
  colnames(series)=c('year','month','day','temp')
  first=series[1,1:3] #first day of the series
  last=series[dim(series)[1],1:3] #last day of the series
  fy0=as.matrix(first[1])
  ly0=as.matrix(last[1])
  days=seq(as.Date(paste(first[1],first[2],first[3],sep="/")), as.Date(paste(last[1],last[2],last[3],sep="/")),"day")
  for (y in fy0:ly0)
  {
    for (m in 1:12)
    {
      if (dim(series[series$year==y & series$month==m,])[1]==0)
      {
        append=data.frame(year=y,month=m,day=32,temp=NA)
        series=rbind(series,append)
      }
    }
  }
  output=list(series=series,fy=fy0,ly=ly0)
  return(output)
}

running_mean <- function(table,rrm,CRU=1) ### running mean for the plots
{
  #colnames(table)[2]=c('value')
  thre=ceiling(rrm/2)
  fy=table[1,1]
  ly=table[dim(table)[1],1]
  if (CRU==0)
  {
    for (y in fy+rrm:ly-rrm)
    {
      
      if (length(table[table$year>=y-rrm & table$year<=y+rrm & !is.na(table$ori),'ori'])>thre)
      {
        table[table$year==y,'ori.rm']=mean(table[table$year>=y-rrm & table$year<=y+rrm,'ori'],na.rm=TRUE)
      }
      if (length(table[table$year>=y-rrm & table$year<=y+rrm & !is.na(table$hom),'hom'])>thre)
      {
        table[table$year==y,'hom.rm']=mean(table[table$year>=y-rrm & table$year<=y+rrm,'hom'],na.rm=TRUE)
      }
      if (length(table[table$year>=y-rrm & table$year<=y+rrm & !is.na(table$hom2),'hom2'])>thre)
      {
        table[table$year==y,'hom2.rm']=mean(table[table$year>=y-rrm & table$year<=y+rrm,'hom2'],na.rm=TRUE)
      }
      if (length(table[table$year>=y-rrm & table$year<=y+rrm & !is.na(table$ben),'ben'])>thre)
      {
        table[table$year==y,'ben.rm']=mean(table[table$year>=y-rrm & table$year<=y+rrm,'ben'],na.rm=TRUE)
      }
    }
  }else
  {
    if (ncol(table)==2)
    {
      table$annmean.rm=filter.cru(rrm,table$annmean)$tslow
    }else
    {
      table$ori.rm=filter.cru(rrm,table$ori)$tslow
      table$hom.rm=filter.cru(rrm,table$hom)$tslow
      if(length(table$hom2)==length(table$ori)){table$hom2.rm=filter.cru(rrm,table$hom2)$tslow}
      if(length(table$ben)==length(table$ori)){table$ben.rm=filter.cru(rrm,table$ben)$tslow}
    }
  }
  return(table)
}

running_mean_2 <- function(table,rrm=5)
{
  thre=ceiling(rrm/2)
  if (ncol(table)==2)
  {
    fy=table[1,1]
    ly=table[dim(table)[1],1]
    table[,3]=filter.cru(rrm,table[,2])$tslow
  }
  return(table[,3])
}



select_season <- function(series,season)
{
  if (season=='djf')
  {
    series=series[series$month%in%c(12,1,2),]
  }
  if (season=='mam')
  {
    series=series[series$month%in%c(3,4,5),]
  }
  if (season=='jja')
  {
    series=series[series$month%in%c(6,7,8),]
  }
  if (season=='son')
  {
    series=series[series$month%in%c(9,10,11),]
  }
  return(series)
}

sen_slope <- function(vec)
{
  ll=length(vec)
  sen_slopes=matrix(NA,ll,ll)
  for (r in 1:ll-1)
  {
    for (s in (r+1):ll)
    {
      sen_slopes[r,s]=round((vec[s]-vec[r])/(s-r),digits=3)
    }
  }
  trend=wise_median(sen_slopes)
  return(trend*10)
}

series_2_sta <- function(table) {
  aux=data.frame(matrix(integer(0),0,ncol(table)))
  colnames(aux)=colnames(table)
  for (sta_id in unique(table$sta_id))
  {
    row=table[table$sta_id==sta_id,]
    row=row[which.max(row$length),]
    row=row[which.max(row$stop),][1,]
    aux=rbind(aux,row)
  }
  return (aux)
}

smoothen_nqx12 <- function(matrice) ### smoothens out matrices with 12 columns 
{
  nq1=dim(matrice)[1]
  smooth=matrix(NA,nq1,12)
  for (m in 1:12)
  {
    if (m==1)
    {
      for (q in 1:nq1)
      {
        if (q==1)
        {
          smooth[q,m]=(matrice[1,12]+2*matrice[1,1]+matrice[1,2]+matrice[2,1])/5
        }else if (q==nq1)
        {
          smooth[q,m]=(matrice[nq1,12]+2*matrice[nq1,1]+matrice[nq1,2]+matrice[nq1-1,1])/5
        }else
        {
          smooth[q,m]=(matrice[q,12]+matrice[q,1]+matrice[q,2]+matrice[q-1,1]+matrice[q+1,1])/5
        }
      }
    }else if (m==12)
    {
      for (q in 1:nq1)
      {
        if (q==1)
        {
          smooth[q,m]=(matrice[1,11]+2*matrice[1,12]+matrice[1,1]+matrice[2,12])/5
        }else if (q==nq1)
        {
          smooth[q,m]=(matrice[nq1,11]+2*matrice[nq1,12]+matrice[nq1,1]+matrice[nq1-1,12])/5
        }else
        {
          smooth[q,m]=(matrice[q,11]+matrice[q,12]+matrice[q,1]+matrice[q-1,12]+matrice[q+1,12])/5
        }
      }
    }else
    {
      for (q in 1:nq1)
      {
        if (q==1)
        {
          smooth[q,m]=(matrice[1,m-1]+2*matrice[1,m]+matrice[1,m+1]+matrice[2,m])/5
        }else if (q==nq1)
        {
          smooth[q,m]=(matrice[nq1,m-1]+2*matrice[nq1,m]+matrice[nq1,m+1]+matrice[nq1-1,m])/5
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

step095 <- function (value)
{
  if (value>=0.95) {out=2.5}
  else  {out=0.6}
  return(out)
}

wise_median <- function(x,digits=1) ### median, but in case there are too many decimals after the point, the value is rounded up or down in the direction of the mean
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

