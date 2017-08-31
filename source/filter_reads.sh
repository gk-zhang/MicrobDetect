#/bin/bash

#### the function is to filter out the mapped reads from the alignment result file(*.sam) of bwa
#### the reads are all considered as single reads
#### only usefule fields of each line of the mapped reads will be copied in the output file
#### to run the program, use bash instead of sh

function filter_reads {

  inputfile=$1 #input *.sam file, the results from bwa mapping
  mismatch_allowed=$2 #allowed mismatch, e.g.=0, means no mismatch is allowed, =2, means maximal 2 bases mismatch are allowed. =1000,which is much larger than the read length, means any mismatch is allowed
  keep_sq=$3 #=n, means @SQ lines will not be presented in the output; =y, means @SQ lines will be presented in the output
  
  # to filter the mapped reads
  command_m="awk -v mis=\"\$mismatch_allowed\" -v OFS='\t' '{if (\$1==\"@SQ\") {print \$0};if (\$6 ~ /M/) {split(\$6,a,\"[0-9]*M\");sum=0;for (i=1;i<=length(a);i++) {if (a[i] ~ /[0~9]*/) sum=sum+substr(a[i],0,length(a[i])-1)};if (sum<=mis) { if (\$12 ~ /NM/) {split(\$12,b,\":\");if (b[3]<=mis) print \$0}}}}'"

  if [[ $keep_sq -eq y ]]; then
    cat $inputfile |eval $command_m
  else
    if [[ $keep_sq -eq n ]]; then
      cat $inputfile |grep -v "@SQ" |eval $command_m
    fi
  fi    
}

filter_reads $1 $2
