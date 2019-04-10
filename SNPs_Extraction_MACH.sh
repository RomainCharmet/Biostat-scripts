#This script is used to extract specific SNPs from MACH imputation files
grep "SNP\|SNP1\|...\|SNPn" file.info
#Lists lines with the chosen SNPs
awk -v ORS='\\|' '{print $colonneSNP}'file
#Lists SNPs to be extracted with a grep-compatible format (SNP1\|...\|SNPn)
grep -n "SNP1\|...\|SNPn" file.info | awk -F ':' '{print $1}'
#Lists line number in info files for SNPs, the corresponding column in dose file will be line+1
var="$(grep -n "SNP\|SNP1\|...\|SNPn" file.info | awk -F ':' '{print $1}')"
#Same as above but stored in a variable
a=( $var )
for i in {0..1}
do
echo $(expr ${a[$i]} + 1)
done | awk -v ORS=',\$' '{print $1}'
#Lists SNPs column numbers with a ,$ split
awk -F '\t' '{print $1,$2,$colSNP1,...,$colSNPn}' file.dose
#Display doses for chosen SNPs on screen (DON'T DO IT ON WITH LARGE SAMPLES)

#The following function extracts SNP doses from the base file into a new file
_parse()
{
    local fs="$1"
    local ofs="$2"
    shift 2
    local _s=
    local f

    for f; do
        _s="${_s}\$${f},"
    done
    awk -F"$fs" -v OFS="$ofs" "{ print ${_s%,} }"
}

#Remove chromosomes in the for as needed
#SNPs in tyhe chr$i.txt files should be displayed by name on a single column (rsID or chr:bp depending on the used format)
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22
do
cp /san2/home/charmet/COMMON/SHARE/Beata/BeaPourRomain/Imputation_JDRF/CHR$i/results_chr$i.info .
cp /san2/home/charmet/COMMON/SHARE/Beata/BeaPourRomain/Imputation_JDRF/CHR$i/results_chr$i.dose.gz .
gunzip -d results_chr$i.dose.gz
l="$(awk -v ORS='\\|' '{print $1}' chr$i.txt)"
l=${l%\\|}
var="$(grep -n $l results_chr$i.info | awk -F ':' '{print $1}')"
k='SNP\|'$l
grep $k results_chr$i.info > chr$i.info
a=( $var )
n=${#a[@]}
n=$(expr $n - 1)
p="$(for (( c=0; c<=$n; c++ ))
do
echo $(expr ${a[$c]} + 1)
done | awk -v ORS=' ' '{print $1}')"
p=${p%' '}
p='1 2 '$p
_parse ' ' '\t' $p < results_chr$i.dose > chr$i.dose
rm -f results_chr$i.dose results_chr$i.info
done
