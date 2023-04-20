print("Starting PopOff")
#setwd("/root")
if(!require(RSQLite)){
  install.packages("RSQLite")
}
library(RSQLite)
if(!require(shiny)){
  install.packages("shiny")
}
library(shiny)
if(!require(DT)){
  install.packages("DT")
}
library(DT)
if(!require(purrr)){
  install.packages("purrr")
}
library(purrr)
if(!require(stringi)){
  install.packages("stringi")
}
library(stringi)
if(!require(ggplot2)){
  install.packages("ggplot2")
}
library(ggplot2)
if(!require(scales)){
  install.packages("scales")
}
library(scales)
if(!require(ggpubr)){
  #  install.packages("pbkrtest")
  install.packages("ggpubr")
}
library(ggpubr)
if(!require(dplyr)){
  install.packages("dplyr")
}
library(dplyr)
if(!require(tidyr)){
  install.packages("tidyr")
}
library(tidyr)
#if(!require(plyr)){
#  install.packages("plyr")
#}
#library(plyr)
if(!require(grid)){
  install.packages("grid")
}
library(grid)
if(!require(gridExtra)){
  install.packages("gridExtra")
}
library(gridExtra)

callback_js = JS(
  "table.column(1).nodes().to$().css({cursor: 'pointer'});",
  "",
  "// make the table header of the nested table",
  "var format = function(d, childId){",
  "  if(d != null){",
  "    var html = ",
  "      '<table class=\"display_compact_hover\" id=\"' + childId + '\"><thead><tr>';",
  "    for (var key in d[d.length-1][0]) {",
  "      html += '<td>' + key + '</td>';",
  "    }",
  "    html += '</tr></thead></table>'",
  "    return html;",
  "  } else {",
  "    return '';",
  "  }",
  "};",
  "",
  "// row callback to style the rows of the child tables",
  "var rowCallback = function(row, dat, displayNum, index){",
  "  if($(row).hasClass('odd')){",
  "    $(row).css('background-color', 'white');",
  "    $(row).hover(function(){",
  "      $(this).css('background-color', 'white');",
  "    }, function() {",
  "      $(this).css('background-color', 'white');",
  "    });",
  "  } else {",
  "    $(row).css('background-color', 'white');",
  "    $(row).hover(function(){",
  "      $(this).css('background-color', 'white');",
  "    }, function() {",
  "      $(this).css('background-color', 'white');",
  "    });",
  "  }",
  "};",
  "",
  "// header callback to style the header of the child tables",
  "var headerCallback = function(thead, data, start, end, display){",
  "  $('th', thead).css({",
  "    'border-top': '3px solid indigo',",
  "    'color': 'black',",
  "    'background-color': 'white'",
  "  });",
  "};",
  "",
  "// make the datatable",
  "var format_datatable = function(d, childId){",
  "  var dataset = [];",
  "  var n = d.length - 1;",
  "  for(var i = 0; i < d[n].length; i++){",
  "    var datarow = $.map(d[n][i], function (value, index) {",
  "      return [value];",
  "    });",
  "    dataset.push(datarow);",
  "  }",
  "  var id = 'table#' + childId;",
  "    var subtable = $(id).DataTable({",
  "            'data': dataset,",
  "            'autoWidth': true,",
  "            'deferRender': true,",
  "            'info': false,",
  "            'lengthChange': false,",
  "            'ordering': d[n].length > 1,",
  "            'order': [],",
  "            'paging': false,",
  "            'scrollX': false,",
  "            'scrollY': false,",
  "            'searching': false,",
  "            'sortClasses': false,",
  "            'rowCallback': rowCallback,",
  "            'headerCallback': headerCallback,",
  "            'columnDefs': [",
  "              {targets: -1, visible: false},",
  "              {targets: 0, orderable: false, className: 'details-control'},",
  "              {targets: '_all', className: 'dt-center'}",
  "             ]",
  "          }).column(0).nodes().to$().css({cursor: 'pointer'});",
  "  }",
  "",
  "// display the child table on click",
  "table.on('click', 'td.details-control', function(){",
  "  var tbl = $(this).closest('table'),",
  "      tblId = tbl.attr('id'),",
  "      td = $(this),",
  "      row = $(tbl).DataTable().row(td.closest('tr')),",
  "      rowIdx = row.index();",
  "  if(row.child.isShown()){",
  "    row.child.hide();",
  "    td.html('&oplus;');",
  "  } else {",
  "    var childId = tblId + '-child-' + rowIdx;",
  "    row.child(format(row.data(), childId)).show();",
  "    td.html('&CircleMinus;');",
  "    format_datatable(row.data(), childId);",
  "  }",
  "});")
rev_com <- function(string){
  return(stri_reverse(chartr("ATGCatgcUuNn","TACGtacgAaNn",string)))
}
find_overlaps <- function(line,PAM,ALL_vars){
  line <- data.frame(t(line))
  #print(line)
  if (line$Sense == '-'){
    line$gRNA <- paste(rev_com(PAM),rev_com(line$gRNA),sep = "")
    line$Ref_seq <- paste(rev_com(PAM),rev_com(line$Ref_seq),sep = "")
    #line$Start <- as.numeric(line$Start)-nchar(PAM)
  }else{
    line$gRNA <- paste(line$gRNA,PAM,sep="")
    line$Ref_seq <- paste(line$Ref_seq,PAM,sep="")
    #line$Stop <- as.numeric(line$Stop)+nchar(PAM)
  }
  #subset <- as.numeric(ALL_vars$Loc) >= as.numeric(as.character(line$Start)) & as.numeric(ALL_vars$LOC) <= as.numeric(as.character(line$Stop))
  vars <- ALL_vars[(as.character(ALL_vars$Chr) == as.character(line$Chr)) & (as.numeric(as.character(ALL_vars$Loc)) >= as.numeric(as.character(line$Start))) & (as.numeric(as.character(ALL_vars$Loc)) <= as.numeric(as.character(line$Stop))),]
  #vars <- ALL_vars[(as.numeric(as.character(ALL_vars$Loc)) >= 1002387) & (as.numeric(as.character(ALL_vars$LOC)) <= 1002414),]
  #print(vars)
  #print('done')
  vars$x <- as.numeric(vars$Loc)-as.numeric(line$Start)+1
  if (nrow(vars) == 0){
    return(NA)
  }else if(nrow(vars) >= 1){
    Effect <- apply(vars,1,find_effect, start = line$Start, gRNA = line$gRNA, ref = line$Ref_seq, ref_pam = line$Ref_pam, sense = line$Sense, PAM = PAM)
    
    subsub <- apply(subset(vars, select = c('XX_AF','XY_AF','AF_African.African_American','AF_Amish','AF_Ashkenazi_Jewish','AF_East_Asian','AF_European_Finnish','AF_European_Non.Finnish','AF_Latino.Admixed_American','AF_Middle_Eastern','AF_Other','AF_South_Asian')),1,sub_sub)
    simple_var <- apply(vars,1,simplify_var)
    var_loc <- apply(vars,1,visulise_var_location, seq = line$Ref_seq, sense = line$Sense, PAM = PAM)
    vars <- cbind(" " = "&oplus;", Effect = Effect, Variant = simple_var, 'Variant Location' = var_loc , subset(vars, select = c('AF','Popmax')), details = I(subsub))
    #vars <- cbind(" " = "&oplus;", Effect = Effect, Variant = simple_var, 'Variant Location' = var_loc , subset(vars, select = c('AF','Popmax')), details = "NA")
    #print("done")
    return(I(vars))
  }
}
visulise_var_location <- function(line,seq,sense,PAM){
  line <- data.frame(t(line))
  if (sense == "-"){
    seq <- chartr("ATGCatgcUuNn","TACGtacgAaNn",seq)
  }
  tmp = as.list(strsplit(seq,"")[[1]])
  tmp[as.numeric(line$x)]= paste('<b>',line$Alt,'</b>',sep="")
  if (sense == "-"){
    tmp <- rev(tmp)
  }
  tmp[length(tmp)-nchar(PAM)+1] = paste('<u><i>',tmp[length(tmp)-nchar(PAM)+1],sep="")
  tmp[length(tmp)] = paste(tmp[length(tmp)],'</u></i>',sep="")
  return(paste(tmp,collapse=""))
}
simplify_var <- function(line){
  line <- data.frame(t(line))
  return(paste(line$Chr,line$Loc,paste(line$Ref,line$Alt,sep = ">"),sep = ":"))
}
match <- function(a,b){
  count <- 0
  score_mat <- c('A' = 'ANRMWDHV',
                 'T' = 'TNYKWBDH',
                 'G' = 'GNRKSADV',
                 'C' = 'CNYMSBHV',
                 'N' = 'N',
                 'R' = 'RAGN',
                 'Y' = 'YCTN',
                 'K' = 'KGTN',
                 'M' = 'MACN',
                 'S' = 'SCGN',
                 'W' = 'WATN',
                 'B' = 'BCGTYKSN',
                 'D' = 'DAGTRKWN',
                 'H' = 'HACTYMWN',
                 'V' = 'VACGRMSN'
  )
  a <- unlist(strsplit(as.character(a),NULL))
  b <- unlist(strsplit(as.character(b),NULL))
  #print(a)
  #print(b)
  for (x in 1:length(b)){
    if (b[[x]] %in% unlist(strsplit(score_mat[a[[x]]],NULL))){
    }else{
      count <- count+1
    }
  }
  return(count)
}
find_effect <- function(vars,start,gRNA,ref,ref_pam,sense,PAM){
  vars <- lapply(data.frame(t(vars)), as.character)
  #print(vars)
  if (nchar(vars$Ref) != nchar(vars$Alt)){
    if ((sense == '+' & as.numeric(vars$x) > (nchar(gRNA)-nchar(PAM))) | (sense == '-' & as.numeric(vars$x) <= nchar(PAM))){
      old <- match(ref_pam,PAM)
      if (old == 0){
        return('PAM Indel')
      }else{
        return('Non-PAM Indel')
      }
    }else{
      return('Indel')
    }
  }
  x <- as.numeric(vars$Loc)-as.numeric(start)
  sub_gRNA <- substr(gRNA, as.numeric(vars$x), as.numeric(vars$x)+nchar(vars$Ref)-1)
  #sub_ref <- substr(ref, as.numeric(vars$x), as.numeric(vars$x)+nchar(vars$Ref)-1)
  if ((sense == '+' & as.numeric(vars$x) > (nchar(gRNA)-nchar(PAM))) | (sense == '-' & as.numeric(vars$x) <= nchar(PAM))){
    new_pam <- unlist(strsplit(as.character(ref_pam), NULL))
    if (sense == '+'){
      x <- 1+nchar(PAM)-(nchar(gRNA)-x)
      new_pam[[x]] <- vars$Alt
    }else{
      x <- nchar(PAM)-x
      new_pam[[x]] <- rev_com(vars$Alt)
    }
    old <- match(ref_pam,PAM)
    new <- match(paste(new_pam,sep=''),PAM)
    if (old == new){
      if (old == 0){
        return("PAM Neutral")
      }else{
        return("Non-PAM Neutral Guide Substitution")
      }
    }
    if (old > new){
      if (new == 0){
        return("Novel PAM")
      }else{
        return("Non-PAM Enhancing Guide Substitution")
      }
    }else{
      if (old == 0){
        return("PAM Loss")
      }else{
        return("Non-PAM Ablative Guide Substitution")
      }
    }
  }else{
    if (sub_gRNA == vars$Ref & sub_gRNA != vars$Alt){
      return("Ablative Guide Substitution")
    }
    if (sub_gRNA != vars$Ref & sub_gRNA == vars$Alt){
      return("Enhancing Guide Substitution")
    }
    if (sub_gRNA != vars$Ref & sub_gRNA != vars$Alt){
      return("Neutral Guide Substitution")
    }else{
      return(NA)
    }
    return("ERROR")
  }
  #return("ERROR")
}

convert_table = function(line){
  deets <- lapply(line[length(line)],convert_subdat)
  return(paste(c(unlist(line[3:length(line)-1]),deets), collapse = '\t'))
}
convert_subdat = function(line){
  if (is.na(line)){
    return('NA\tNA') #<======================= Need to update this later!
  }else{
    tmp <- line[[1]]
    pers <- lapply(tmp[length(tmp)],convert_subsubdat)
    
    deets <- cbind(as.data.frame(tmp[3:length(tmp)-1]),pers)
    ret <- c()
    for(x in 1:length(deets)){
      ret <- c(ret,as.character(lapply(deets[x],paste,collapse=',')[1]))
    }
    return(paste(unname(ret),collapse = '\t'))
  }
}
convert_subsubdat = function(line){
  return(paste(unlist(line$AF), collapse = '|'))
}
count_vars = function(line){
  line <- data.frame(t(line))
  if (anyNA(line$details)){
    return(0)
  }else{
    return(nrow(data.frame(line$details)))
  }
}
add_symbol = function(line){
  line <- data.frame(t(line))
  if (anyNA(line$details)){
    return(' ')
  }else{
    return('&oplus;')
  }
}
sub_sub = function(line){
  df <- data.frame(Population = as.vector(c('XX Allele Frequency','XY Allele Frequency','African/African American Allele Frequency','Amish Allele Frequency','Ashkenazi Jewish Allele Frequency','East Asian Allele Frequency','European Finnish Allele Frequency','European Non-Finnish Allele Frequency','Latino/Admixed American Allele Frequency','Middle Eastern Allele Frequency','Other Populations Allele Frequency','South Asian Allele Frequency')), 'Allele Frequency' = as.vector(unname(unlist(line))), stringsAsFactors = FALSE)
  df$details <- I(NA)
  return(df)
}
translcate_location = function(line){
  line <- data.frame(t(line))
  if (line$Sense == '+'){
    return(paste(line$Chr,paste(line$Start,line$Stop,sep = '-'),sep=":"))
  }else{
    return(paste(line$Chr,paste(line$Stop,line$Start,sep = '-'),sep=":"))
  }
}
make_graphs = function(data){
  ret = c()
  ret$palete <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                           "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#C3D7A4","#293352")
                           
  #ret$theme_no_legend_no_lines <- theme_pubr() + 
  #  theme(plot.title = element_text(size = 20),
  #        legend.position = 'none',
  #        plot.caption = element_text(hjust = 0, face= "italic"),
  #        plot.title.position = "plot",
  #        line = element_blank(),
  #        text = element_blank(),
  #        title = element_blank()
  #        #Here
  #  )
  ret$theme_no_legend_angled <- theme_void() + theme(plot.title = element_text(size = 20),
                                                     legend.position = 'none',
                                                     plot.caption = element_text(hjust = 0, face= "italic"),
                                                     plot.title.position = "plot",
                                                     axis.text.x = element_text(size = 15, angle = -45, hjust = 0.5, vjust = -0.4)
  )
  #top_bar_dat <- data.frame(names = c('No Variant','Overlappig\nVariant(s)'), 
  #                          values = c(length(data$'Number of Population Variants')-length(data[data$'Number of Population Variants' != 0,]$'Number of Population Variants'),length(data[data$'Number of Population Variants' != 0,]$'Number of Population Variants'))
  #                          )
  #top_bar_dat$per <- label_percent(accuracy = 0.01)(top_bar_dat$values/sum(top_bar_dat$values))
  #top_bar_dat$labels <- paste(top_bar_dat$per, top_bar_dat$names, sep="\n")
  #ret$top_bar <- ggplot(top_bar_dat, aes(x="", y=values, fill=names, label = labels)) +
  #  geom_bar(stat="identity", width = 1) +
  #  geom_text(size = 5, position = position_nudge()) +
  #  ret$theme_no_legend_no_lines +
  #  scale_fill_manual(values = ret$palete) +
  #  coord_polar("y", start=0) + coord_flip()
  
  if (nrow(data %>% filter(!is.na(details))) >= 1){
    bottom_bar_dat <- data.frame(t(table(unlist(apply((data %>%filter(!is.na(details))),1,list_effects)))))
    bottom_bar_dat$per <- label_percent(accuracy = 0.01)(bottom_bar_dat$Freq/sum(bottom_bar_dat$Freq))
    bottom_bar_dat <- bottom_bar_dat %>% arrange(desc(as.numeric(Freq)))
    bottom_bar_dat$Var2 <- factor(as.character(bottom_bar_dat$Var2), levels=as.character(bottom_bar_dat$Var2))
    ret$bottom_bar <- ggplot(bottom_bar_dat, aes(x=Var2,y=Freq/sum(Freq)*100,fill=Var2,label=per)) +
      geom_bar(stat = "identity") +
      xlab('') +
      geom_text(size = 7, vjust=-1, colour = "black", position = position_dodge(width = 1)) +
      ret$theme_no_legend_angled +
      ylab("Percentage") + 
      scale_fill_manual(values = ret$palete) + 
      ylim(0,max(bottom_bar_dat$Freq/sum(bottom_bar_dat$Freq)*100)+5) +
      ggtitle("Variant Effects")
    #scale_x_discrete(position = "top") 
    
    
    hist_a <- bin_freq(unlist(apply(data[data$details != 'NA',],1,list_AF)))
    ret$hist_a_plot <- ggplot(hist_a, aes(x=Bin,y=Val)) + geom_bar(stat="identity") + 
      xlab("") + ylab("Count") + 
      ggtitle("Allele Frequency") +
      scale_fill_manual(values = ret$palete) +
      ret$theme_no_legend_angled
    
    hist_b <- bin_freq(unlist(apply(data[data$details != 'NA',],1,list_PMAX)))
    ret$hist_b_plot <- ggplot(hist_b, aes(x=Bin,y=Val)) + geom_bar(stat="identity") + 
      xlab("") + ylab("") + ylab("Count") + 
      ggtitle("Maximum Population Frequency") +
      scale_fill_manual(values = ret$palete) +
      ret$theme_no_legend_angled
    
  }else{
    ret$bottom_bar <- geom_blank()
    ret$hist_a_plot <- geom_blank()
    ret$hist_b_plot <- geom_blank()
  }
  return(ret)
}
list_effects <- function(line){
  #line <- data.frame(t(line))
  tmp <- line$details
  return(tmp$Effect)
}
list_AF <- function(line){
  #line <- data.frame(t(line))
  tmp <- line$details
  return(tmp$AF)
}
list_PMAX <- function(line){
  #line <- data.frame(t(line))
  tmp <- line$details
  return(tmp$Popmax)
}
bin_freq <- function(list){
  return(data.frame(Bin = c('>1%','1-10%','11-20%','21-30%','31-40%','41-50%','51-60%','61-70%','71-80%','81-90%','91-100%'),
                    Val = c(length(list[list<0.01]),length(list[list>=0.01 & list<0.1]),length(list[list>=0.1 & list<0.2]),length(list[list>=0.2 & list<0.3]),length(list[list>=0.3 & list<0.4]),length(list[list>=0.4 & list<0.5]),length(list[list>=0.5 & list<0.6]),length(list[list>=0.6 & list<0.7]),length(list[list>=0.7 & list<0.8]),length(list[list>=0.8 & list<0.9]),length(list[list>=0.9]))        
  ))
}
foramt_input <- function(data,form,PAM){
  if (form == 'Cas-OFFinder_Portable' || form == 'Cas-OFFinder'){
    if(ncol(data) > 6){
      ret <- data.frame(Chr = sub("^chr","",data[,4]))
      ret$Sense <- data[,6]
      ret$Start <- data[,5] +1
      ret$Stop <- data[,5] + nchar(data[,2])
      ret$gRNA <- substr(data[,2],1,(nchar(data[,2])-nchar(PAM)))
      ret$Ref_seq <- substr(data[,3],1,(nchar(data[,3])-nchar(PAM)))
      ret$Ref_pam <- substr(data[,3],(nchar(data[,3])-nchar(PAM)+1),nchar(data[,3]))
      ret$Mismatches <- data[,7]
    }
    else{
      ret <- data.frame(Chr = sub("^chr","",data[,2]))
      ret$Sense <- data[,5]
      ret$Start <- data[,3] + 1
      ret$Stop <- data[,3] + nchar(data[,1])
      ret$gRNA <- substr(data[,1],1,(nchar(data[,1])-nchar(PAM)))
      ret$Ref_seq <- substr(data[,4],1,(nchar(data[,4])-nchar(PAM)))
      ret$Ref_pam <- substr(data[,4],(nchar(data[,4])-nchar(PAM)+1),nchar(data[,4]))
      ret$Mismatches <- data[,6]
    }
    
  }
  if (form == 'Off-Spotter'){
    ret <- data.frame(Chr = as.character(data[,1]))
    ret$Sense <- data[,2]
    ret$Start <- data[,3]
    ret$Stop <- data[,4]
    ret$gRNA <- data[,5]
    ret$Ref_seq <- substr(data[,6],1,(nchar(data[,6])-(nchar(PAM)+1)))
    ret$Ref_pam <- substr(data[,6],(nchar(data[,6])-(nchar(PAM)-1)),nchar(data[,6]))
    ret$Mismatches <- data[,7]
  }
  if (form == 'CRISPOR'){
    ret <- data.frame(Chr = sub("^chr","",data$chrom))
    ret$Sense <- data$strand
    ret$Start <- data$start
    ret$Stop <- data$end
    ret$gRNA <- substr(data$guideSeq,1,(nchar(data$guideSeq)-nchar(PAM)))
    data$offtargetSeq <- substr(data$offtargetSeq,1,(nchar(data$offtargetSeq)-nchar(PAM)))
    ret$Ref_seq <- apply(data,1,format_Mismatch)
    ret$Ref_pam <- PAM
    ret$Mismatches <- data$mismatchCount
  }
  if (form == 'tsv' || form == 'csv'){
    ret <- data.frame(Chr = data[,1])
    ret$Sense <- data[,2]
    ret$Start <- data[,3]
    ret$Stop <- data[,4] # can remove?
    ret$gRNA <- data[,5]
    ret$Ref_seq <- data[,6]
    ret$Ref_pam <- data[,7]
    #ret$Mismatches <- length(Reduce(setdiff, strsplit(c(toupper(data[,5]),toupper(data[,6])),split = "")))
    ret$Mismatches <- apply(data,1,function(x){
      tmp <- data.frame(strsplit(toupper(c(x[6], x[5])),split=""))
      return(nrow(tmp[tmp[,1] != tmp[,2],]))
    })
    
  }
  return(ret)
}
format_Mismatch <- function(line){
  line <- data.frame(t(line))
  tmp <- data.frame(char=strsplit(line$offtargetSeq,split = NULL), val = strsplit(line$mismatchPos,split = NULL))
  colnames(tmp) <- c('char','val')
  ret <- paste(apply(tmp,1,replace_char, tar= '*'),collapse='')
  return(ret)
}
replace_char <- function(line, tar){
  line <- data.frame(t(line))
  if (line$val == tar){
    ret <- chartr("ATGCUN","atgcun",line$char)
  }else{
    ret <- line$char
  }
  return(ret)
}
convert_details <- function(line){
  line <- data.frame(t(line))
  if (is.na(line$details)){
    return('NA,NA,NA,NA,NA')
  }
  else{
    det <- (as.list(line$details)[[1]])
    Effect <- paste(det$Effect,collapse=";")
    Variant <- paste(det$Variant,collapse=";")
    AF <- paste(det$AF,collapse=";")
    Popmax <- paste(det$Popmax,collapse=";")
    All_Fq <- list()
    for (x in 1:nrow(det)){
      All_Fq <- append(All_Fq,paste(det$details[[x]]$Allele.Frequency,collapse="|"))
    }
    return(paste(Effect,Variant,AF,Popmax,paste(All_Fq,collapse = ';'),sep=','))
  }
}

suppressWarnings({ ui <- basicPage(
  tags$head(tags$style("#err{color: red;
                       font-style: bold;}")),
  tabsetPanel(id = "MAIN",
              tabPanel("Upload Data",
                       fluidRow(HTML("<center><h1>PopOff: A Population Variant Annotation Tool for Off-Target Sites</h1></center>"),
                                HTML("<center>Welcome to PopOff a tool for annotating CRISPR off-target (or on-target) sites with overlaping population variants. To get started upload off target sites following the insturctions in the <b>Home</b> tab.</center>"),
                                hr()
                       ),
                       sidebarPanel(selectInput(inputId = "PAM",
                                                label = "PAM Sequence",
                                                list('NGG','NGRRT','NGRRN','TTTV','TTN','OTHER')),
                                    conditionalPanel(condition = "input.PAM == 'OTHER'",
                                                     textInput("PAM_OTH", HTML("<i>please enter PAM</i>"))),
                                    tags$hr(),
                                    selectInput(inputId = "IN_FOR",
                                                label = "Input Source",
                                                list('Cas-OFFinder','CRISPOR','Off-Spotter','tsv','csv')), #Add blast fasta?
                                    conditionalPanel(condition = "input.IN_FOR == 'Cas-OFFinder'",
                                                     HTML("<a href=\"http://www.rgenome.net/cas-offinder/\" target=\"_blank\">Looking for the tool? You can find it here!</a>")),
                                    conditionalPanel(condition = "input.IN_FOR == 'CRISPOR'",
                                                     HTML("<a href=\"http://crispor.tefor.net/\" target=\"_blank\">Looking for the tool? You can find it here!</a>")),
                                    conditionalPanel(condition = "input.IN_FOR == 'Off-Spotter'",
                                                     HTML("<a href=\"https://cm.jefferson.edu/Off-Spotter/\" target=\"_blank\">Looking for the tool? You can find it here!</a>")),
                                    conditionalPanel(condition = "input.IN_FOR == 'tsv' || input.IN_FOR == 'csv'",
                                                     HTML("Please provide in input file with the following columns:<ul><li>1: Chromosome</li><li>2: Sense ('+' or '-')</li><li>3: Locus Start Location</li><li>4: Locus End Location</li><li>5: gRNA sequence</li><li>6: Locus sequence (excluding PAM)</li><li>7: Locus PAM sequence</li></ul>")                     ),
                                    tags$hr(),
                                    fileInput("IN", "Choose Input File",
                                              multiple = FALSE,
                                              accept = c("text/csv/tsv","text/comma-separated-values,text/plain,text/tab-separated-values",".csv,.txt,.tsv")),
                                    fluidRow(column(12,textOutput('err'),
                                                    actionButton("upload", "Upload Data"),
                                                    align = "center")
                                    )
                                    
                                    
                       ),
                       mainPanel(
                         tabsetPanel(id = "UPLOAD",
                                     tabPanel("Home",
                                              p(),
                                              HTML("Nearly 1 in 5 off-target loci overlap with a population variant. In addition, these variants can increase the similarity between guides and off-target loci. Using this tool, you can identify population variants that overlap with off-target loci predicted by tools such as CRISPOR and Cas-OFFinder, determine their impact, and see how common they are in different population groups."),
                                              HTML("For more information on PopOff and why you should be assessing population variation at off-target loci please refer to our publication: "),
                                              HTML("<a href=\"https://www.biorxiv.org/\" target=\"_blank\">Detecting CRISPR Off-Target Loci Altered by Population Variants with PopOff</a>."),
                                              p(),
                                              HTML("To use this tool you first need to identify off-target loci. PopOff can use the direct outputs of these tools or you can provide them as a csv or tsv formatted file."),
                                              p(),
                                              HTML("If you are looking for a link to a particular tool it will appear under the dropdown box to the left when you select it."),
                                              p(),
                                              HTML("You will be able to view your data before running the tool under the <b>Data</b> tab which will show up when you upload some data."),
                                              p(),
                                              HTML("Once you have uploaded your data, click the <b>Run PopOff</b> button on the bottom left to run the tool. You will be automatically directed to your results once it is finished running where you can save these results for later use."),
                                              p(),
                                              HTML("Please be aware that submitting more than 1000 off-target loci can result in quite run times exceeding 10 minutes. For large queries we suggest downloading the local version of the tool, available "),
                                              HTML("<a href=\"https://github.com/\" target=\"_blank\">here</a>"),
                                              p(),
                                              #HTML("Tips and Tricks"),
                                              #p(),
                                              #HTML("You can look for novel off-target loci generated by population variants by using the local version of Cas-OFFinder but replacing nucleotides in the PAM sequence with \"N\" and doing the same in PopOff"),
                                              #p(),
                                              HTML("Have a feature you want to request? Let us know at <a href=\"https://github.com/\" target=\"_blank\">here</a>")
                                     ),
                                     tabPanel("Data",
                                              fluidRow(column(12,
                                                              p(),
                                                              HTML("Please check your input. If it is correct click <b>Run PopOff</b> to run the tool. If needed you can re-upload your input data and this table will update."),
                                                              hr()
                                              )),
                                              fluidRow(column(12,actionButton("button", "Run PopOff"),
                                                              align = "center")
                                              ),
                                              fluidRow(column(12,
                                                              DTOutput("input_table")
                                              )
                                              )
                                     )
                         ),
                         fluidRow(hr()
                         )
                       )
                       
                       
              ),
              tabPanel("Results",
                       tabsetPanel(id = "Results_SUB",
                                   tabPanel("Overview",
                                            fluidRow(column(12,
                                                            h2(textOutput('PER'))),
                                                     align = "center"),
                                            hr(),
                                            fluidRow(column(12,plotOutput("graph"))),
                                            hr(),
                                            h3("Effect Descriptions"),
                                            tabsetPanel(id = "Effect_Explainer",
                                                        tabPanel("Ablative Guide Substitution", #Mismatch Gain
                                                                 HTML("<h4>Population variant generates a new missmatch between the guide and target locus</h4>"),
                                                                 p(),
                                                                 HTML("<body><h4>
                                             <table style=\"width: 50%\">
                                             <tr>
                                             <td style=\"text-align: right\">Guide:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;margin-left: 10px;\">GCTAGCTAGCGCTAGCTGAC-NGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Reference Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAG<b>C</b>GCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-CGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Variant Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAG<b><FONT COLOR=\"RED\">T</b></FONT>GCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-CGG</td>
                                             </tr>
                                             </table>
                                             </h4></body>
                                             ")
                                                        ),
                                                        tabPanel("Enhancing Guide Substitution", #Mismatch Loss
                                                                 HTML("<h4>Population variant removes a missmatch between the guide and target locus</h4>"),
                                                                 p(),
                                                                 HTML("<body><h4>
                                             <table style=\"width: 50%\">
                                             <tr>
                                             <td style=\"text-align: right\">Guide:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;margin-left: 10px;\">GCTAGCTAGCGCTAGCTGAC-NGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Reference Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<b><FONT COLOR=\"RED\">C</FONT></b>CTGAC-CGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Variant Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<b>G</b>CTGAC-CGG</td>
                                             </tr>
                                             </table>
                                             </h4></body>
                                             ")
                                                        ),
                                                        tabPanel("Neutral Guide Substitution", #Mismatch Silent
                                                                 HTML("<h4>Population variant substitues a missmatch between the guide and target locus for another</h4>"),
                                                                 p(),
                                                                 HTML("<body><h4>
                                             <table style=\"width: 50%\">
                                             <tr>
                                             <td style=\"text-align: right\">Guide:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;margin-left: 10px;\">GCTAGCTAGCGCTAGCTGAC-NGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Reference Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<b><FONT COLOR=\"RED\">C</FONT></b>CTGAC-CGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Variant Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<b><FONT COLOR=\"RED\">T</FONT></b>CTGAC-CGG</td>
                                             </tr>
                                             </table>
                                             </h4></body>
                                             ")
                                                        ),
                                                        tabPanel("Indel",
                                                                 HTML("<h4>Population variant removes or adds bases in the guide region at the target locus</h4>"),
                                                                 p(),
                                                                 HTML("<body><h4>
                                             <table style=\"width: 50%\">
                                             <tr>
                                             <td style=\"text-align: right\">Guide:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;margin-left: 10px;\">GCTAGCTAGCGCTAGCTGAC-NGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Reference Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAG<b>CGC</b>TA<FONT COLOR=\"RED\">C</FONT>CTGAC-CGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Variant Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAG<b>---</b>TA<FONT COLOR=\"RED\">T</FONT>CTGAC-CGG</td>
                                             </tr>
                                             </table>
                                             </h4></body>
                                             ")
                                                        ),tabPanel("PAM Loss",
                                                                   HTML("<h4>Population variant removes the PAM site of the target locus</h4>"),
                                                                   p(),
                                                                   HTML("<body><h4>
                                             <table style=\"width: 50%\">
                                             <tr>
                                             <td style=\"text-align: right\">Guide:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;margin-left: 10px;\">GCTAGCTAGCGCTAGCTGAC-NGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Reference Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-C<b>G</b>G</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Variant Allele:&ensp;</td>
                                              <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-C<FONT COLOR=\"RED\"><b>T</b></FONT>G</td>
                                             </tr>
                                             </table>
                                             </h4></body>
                                             ")
                                                        ),tabPanel("PAM Neutral",
                                                                   HTML("<h4>Population variant occurs in the PAM site of the target locus but does not impact it</h4>"),
                                                                   p(),
                                                                   HTML("<body><h4>
                                             <table style=\"width: 50%\">
                                             <tr>
                                             <td style=\"text-align: right\">Guide:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;margin-left: 10px;\">GCTAGCTAGCGCTAGCTGAC-NGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Reference Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-<b>C</b>GG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Variant Allele:&ensp;</td>
                                              <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-<FONT COLOR=\"RED\"><b>A</b></FONT>GG</td>
                                             </tr>
                                             </table>
                                             </h4></body>
                                             ")
                                                        ),tabPanel("PAM Indel",
                                                                   HTML("<h4>Population variant removes or adds bases in the PAM region at the locus</h4>"),
                                                                   p(),
                                                                   HTML("<body><h4>
                                             <table style=\"width: 50%\">
                                             <tr>
                                             <td style=\"text-align: right\">Guide:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;margin-left: 10px;\">GCTAGCTAGCGCTAGCTGAC-N--GG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Reference Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-<b>C--</b>GG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Variant Allele:&ensp;</td>
                                              <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-<FONT COLOR=\"RED\"><b>CAT</b></FONT>GG</td>
                                             </tr>
                                             </table>
                                             </h4></body>
                                             ")
                                                        ),tabPanel("Novel PAM",
                                                                   HTML("<h4>Population variant generates PAM site of the locus</h4>"),
                                                                   p(),
                                                                   HTML("<body><h4>
                                             <table style=\"width: 50%\">
                                             <tr>
                                             <td style=\"text-align: right\">Guide:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;margin-left: 10px;\">GCTAGCTAGCGCTAGCTGAC-NGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Reference Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-C<FONT COLOR=\"RED\"><b>T</b></FONT>G</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Variant Allele:&ensp;</td>
                                              <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-C<b>G</b>G</td>
                                             </tr>
                                             </table>
                                             </h4></body>
                                             ")
                                                        ),tabPanel("Non-PAM Enhancing Guide Substitution",
                                                                   HTML("<h4>Population variant causes the prospective PAM site to be more similar to the required PAM sequence but does not generate a novel PAM site</h4>"),
                                                                   p(),
                                                                   HTML("<body><h4>
                                             <table style=\"width: 50%\">
                                             <tr>
                                             <td style=\"text-align: right\">Guide:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;margin-left: 10px;\">GCTAGCTAGCGCTAGCTGAC-NGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Reference Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-C<FONT COLOR=\"RED\"><b>T</b>A</FONT></td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Variant Allele:&ensp;</td>
                                              <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-C<b>G</b><FONT COLOR=\"RED\">A</FONT></td>
                                             </tr>
                                             </table>
                                             </h4></body>
                                             ")
                                                        ),tabPanel("Non-PAM Ablative Guide Substitution",
                                                                   HTML("<h4>Population variant causes the prospective PAM site to be less similar to the required PAM sequence for a locus that does not have a PAM site</h4>"),
                                                                   p(),
                                                                   HTML("<body><h4>
                                             <table style=\"width: 50%\">
                                             <tr>
                                             <td style=\"text-align: right\">Guide:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;margin-left: 10px;\">GCTAGCTAGCGCTAGCTGAC-NGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Reference Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-C<b>G</b><FONT COLOR=\"RED\">A</FONT></td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Variant Allele:&ensp;</td>
                                              <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-C<FONT COLOR=\"RED\"><b>T</b>A</FONT></td>
                                             </tr>
                                             </table>
                                             </h4></body>
                                             ")
                                                        ),tabPanel("Non-PAM Neutral Guide Substitution",
                                                                   HTML("<h4>Population variant occurs in prospective PAM for a locus that does not have a PAM site but does not alter it to be more or less like the required PAM sequence</h4>"),
                                                                   p(),
                                                                   HTML("<body><h4>
                                             <table style=\"width: 50%\">
                                             <tr>
                                             <td style=\"text-align: right\">Guide:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;margin-left: 10px;\">GCTAGCTAGCGCTAGCTGAC-NGG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Reference Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-<b>C</b><FONT COLOR=\"RED\">T</FONT>G</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Variant Allele:&ensp;</td>
                                              <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-<b>A</b><FONT COLOR=\"RED\">T</FONT>G</td>
                                             </tr>
                                             </table>
                                             </h4></body>
                                             ")
                                                        ),tabPanel("Non-PAM Indel",
                                                                   HTML("<h4>Population variant removes or adds bases in the prospective PAM for a locus that does not have a PAM site</h4>"),
                                                                   p(),
                                                                   HTML("<body><h4>
                                             <table style=\"width: 50%\">
                                             <tr>
                                             <td style=\"text-align: right\">Guide:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;margin-left: 10px;\">GCTAGCTAGCGCTAGCTGAC-N--GG</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Reference Allele:&ensp;</td>
                                             <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-<b>C--</b><FONT COLOR=\"RED\">A</FONT>G</td>
                                             </tr>
                                             <tr>
                                             <td style=\"text-align: right\">Variant Allele:&ensp;</td>
                                              <td style=\"font-family: Consolas, monaco, monospace;\">GC<FONT COLOR=\"RED\">G</FONT>AGCTAGCGCTA<FONT COLOR=\"RED\">C</FONT>CTGAC-<FONT COLOR=\"RED\"><b>CAT</b>A</FONT>G</td>
                                             </tr>
                                             </table>
                                             </h4></body>
                                             ")
                                                        )
                                            )
                                   ),
                                   tabPanel("Loci of Concern",
                                            HTML("<center><h3>This table contains the results for all loci overlaping with a 'Enhancing Guide Substitution' or 'Novel PAM' variant</h3></center>"),
                                            fluidRow(column(12, offest = 1, DTOutput("concern_table")))
                                   ),
                                   tabPanel("All Loci",
                                            HTML("<center><h3>This table contains the results for all loci</h3></center>"),
                                            fluidRow(column(12, offest = 1, DTOutput("table")))
                                   ),
                                   
                                   fluidRow(column(12,downloadButton("download_button", label = "Download All Loci as Table"),downloadButton("download_graph", label = "Download Graphs")))
                                   
                       )
              )
  )
  
)
})
server <- function(input, output,session) {
  DATA <- reactiveValues()
  hideTab(inputId = "MAIN", target = "Results")
  hideTab(inputId = "UPLOAD", target = "Data")
  #Add code to check for input
  observeEvent(input$upload,{
    req(input$IN)
    if (is.null(tryCatch(
      {
        if (input$IN_FOR %in% c('Cas-OFFinder','Off-Spotter','Cas-OFFinder_Portable')){
          DATA$in_table <- foramt_input(read.table(input$IN$datapath, sep = "\t", stringsAsFactors = FALSE),input$IN_FOR,chartr("Uu","Tt",input$PAM))
        }else if (input$IN_FOR == 'csv'){
          DATA$in_table <- foramt_input(read.table(input$IN$datapath, sep = ',', stringsAsFactors = FALSE, header = FALSE),input$IN_FOR,chartr("Uu","Tt",input$PAM))
        }else if (input$IN_FOR == 'tsv'){
          DATA$in_table <- foramt_input(read.table(input$IN$datapath, sep = '\t', stringsAsFactors = FALSE, header = FALSE),input$IN_FOR,chartr("Uu","Tt",input$PAM))
        }else{
          DATA$in_table <- foramt_input(read.table(input$IN$datapath, sep = '\t', stringsAsFactors = FALSE, header = TRUE),input$IN_FOR,chartr("Uu","Tt",input$PAM))
        }
      },
      error = function(e) {
        return(NULL)
      },
      warning = function(w) {
        return(NULL)
      },
      finally = {
        #print('X')
      }
    ))){
      output$err <- renderText({"ERROR: could not load input. Please check that it is properly formatted."})
    }else{
      output$err <- NULL
      output$input_table <- renderDT({
        datatable(DATA$in_table, selection = "none", rownames = FALSE
        )
      })
      showTab(inputId = "UPLOAD", target = "Data")
      updateTabsetPanel(session, "UPLOAD", selected = "Data")
    }
  })
  
  observeEvent(input$button,{
    showModal(modalDialog("Processing input....",footer="This can take 5-10 minuites"))
    db <- dbConnect(SQLite(), dbname = "POPOFF_DB")
    DATA$PROSS <- DATA$in_table
    DATA$PROSS$Chr <- sub("^[Cc]hr","",DATA$PROSS$Chr)
    DATA$PROSS$Range <- apply(DATA$PROSS,1,make_range)
    vars <- data.frame()
    print("Finding overlaps")
    for (x in unique(DATA$PROSS$Chr)){
      vars <- rbind(vars,dbGetQuery(conn = db, paste("SELECT * FROM '",as.character(x),"' WHERE IDX IN (",paste(unique(unlist(unnest(DATA$PROSS[DATA$PROSS$Chr == x,],Range)$Range)),collapse = ','),")",sep='')))
    }
    vars$IDX <- NULL
    vars <- vars %>% distinct()
    vars$Chr <- sub("^chr","",vars$Chr)
    print("Annoating overlaps")
    DATA$PROSS$details <- apply(DATA$PROSS,1,find_overlaps, PAM = chartr("Uu","Tt",input$PAM), ALL_vars = vars)
    print("Did Details")
    rm(vars)
    DATA$PROSS$Range <- NULL
    dbDisconnect(db)
    DATA$PROSS$Chr <- paste('chr',DATA$PROSS$Chr, sep = '')
    DATA$Table <- cbind(' ' = apply(DATA$PROSS,1,add_symbol), 'Number of Population Variants' = apply(DATA$PROSS,1,count_vars) , Location = apply(DATA$PROSS,1,translcate_location), 'Sequence in Reference' = DATA$PROSS$Ref_seq, 'Number of Mismatches' = DATA$PROSS$Mismatches, subset(DATA$PROSS, select = c('gRNA','details')))
    print("Did Table")
    DATA$Graphs <- make_graphs(DATA$Table)
    print("Did Graphs")
    showTab(inputId = "MAIN", target = "Results")
    updateTabsetPanel(session, "MAIN", selected = "Results")
    removeModal()
  })
  
  output$download_button <- downloadHandler(
    filename = "PopOff_results.txt",
    content = function(file) {
      #writeLines(paste(apply(DATA$data,1, convert_table), collapse = "\n"), file)
      Export <- DATA$Table
      Export$`Effects,Variants,AF,Popmax,Poplation Frequencies (XX|XY|African/African American|Amish|Ashkenazi Jewish|East Asian|European Finnish|European Non-Finnish|Latino/Admixed American|Middle Eastern|Other Populations|South Asian)` <- apply(Export,1,convert_details)
      Export$` ` <- NULL
      Export$details <- NULL
      write.csv(Export, file, quote = FALSE, row.names = FALSE)
      #write.table(paste(DATA$data,collapse=", "), file,col.names=FALSE)
    }
  )
  
  output$download_graph <- downloadHandler(
    filename = "PopOff_graphs.pdf",
    content = function(file) {
      #graphs <- grid.arrange(DATA$Graphs$bottom_bar,DATA$Graphs$hist_a_plot,DATA$Graphs$hist_b_plot,
      #                       layout_matrix = rbind(
      #                         c(1,1,1,1,2,2),
      #                         c(1,1,1,1,2,2),
      #                         c(1,1,1,1,3,3),
      #                         c(1,1,1,1,3,3)
      #                       )
      #)
      ggexport(DATA$Graphs$bottom_bar,DATA$Graphs$hist_a_plot,DATA$Graphs$hist_b_plot, filename = file)
      #ggsave(file, plot = DATA$Graphs$pie_a_plot)
    }
  )
  output$graph <- renderPlot({
    req(input$IN)
    #ggarrange(Graphs$pie_a_plot,Graphs$pie_b_plot,ggarrange(Graphs$hist_a_plot,Graphs$hist_b_plot,nrow = 2),ncol = 3)
    #ggarrange(DATA$Graphs$top_bar,DATA$Graphs$bottom_bar,ggarrange(DATA$Graphs$hist_a_plot,DATA$Graphs$hist_b_plot,nrow = 2),ncol = 3)
    grid.arrange(DATA$Graphs$bottom_bar,DATA$Graphs$hist_a_plot,DATA$Graphs$hist_b_plot,
                 layout_matrix = rbind(
                   c(1,1,1,1,2,2),
                   c(1,1,1,1,2,2),
                   c(1,1,1,1,3,3),
                   c(1,1,1,1,3,3)
                 )
    )
    #ggarrange(DATA$Graphs$pie_a_plot,DATA$Graphs$pie_b_plot,ggarrange(DATA$Graphs$hist_a_plot,DATA$Graphs$hist_b_plot,nrow = 2),ncol = 3)
  })
  output$table <- renderDT({
    req(input$IN)
    DT::datatable(DATA$Table, callback = callback_js, escape = -2, editable = FALSE, selection = 'none', #rownames = FALSE,
                  options = list(
                    columnDefs = list(
                      list(visible = FALSE, targets = ncol(DATA$Table)),
                      list(orderable = FALSE, className = 'details-control', targets = 1),
                      list(className = "dt-center", targets = "_all")
                    )
                  )
    )
  })
  output$concern_table <- renderDT({
    req(input$IN)
    DT::datatable(DATA$Table[apply(DATA$Table,1,function(x){
      ifelse(is.na(x['details']),
             FALSE,
             ifelse('Enhancing Guide Substitution' %in% sapply(x['details'],"[",'Effect')[[1]] || 'Novel PAM' %in% sapply(x['details'],"[",'Effect')[[1]],TRUE,FALSE)
      )}),],
      #print(sapply(test[!is.na(test$details),]$details,"[[",'Effect')) #[[1]]$Effect)
      #DATA$Table['Enhancing Guide Substitution' %in% DATA$Table$details || 'Novel PAM' %in% DATA$Table$details,],
      callback = callback_js, escape = -2, editable = FALSE, selection = 'none', #rownames = FALSE,
      options = list(
        columnDefs = list(
          list(visible = FALSE, targets = ncol(DATA$Table)),
          list(orderable = FALSE, className = 'details-control', targets = 1),
          list(className = "dt-center", targets = "_all")
        )
      )
    )
  })
  output$PER <- renderText({
    req(input$IN)
    HTML(paste("Of your loci,",
               as.character(label_percent(accuracy = 0.01)((length(DATA$Table$'Number of Population Variants')-length(DATA$Table[DATA$Table$'Number of Population Variants' == 0,]$'Number of Population Variants'))/length(DATA$Table$'Number of Population Variants'))),
               "overlap with population variants",
               sep = " ")
    )
  })
  
  
}
make_range <- function(line){
  line <- data.frame(t(line))
  #print(as.numeric(line[1,]$Start):as.numeric(line[1,]$Stop))
  return(paste(as.numeric(as.character(line$Start)):as.numeric(as.character(line$Stop)),sep = ',',collapse=','))
}

#setwd("C:/Users/csam488/Downloads")
shinyApp(ui = ui, server = server)
