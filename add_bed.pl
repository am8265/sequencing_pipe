# done in way too much of a hurry...
use strict; use File::Basename; use Data::Dumper;
use lib "/nfs/central/home/dh2880/.bin";
use lib "/home/dh2880/.bin";
use ARGH;
my$n=$ARGV[0];
my($s_a,    $s_X,   $s_Y,   $c,$outdir)=(0,0,0,0,'/nfs/goldsteindata/refDB/captured_regions/Build37/');
my @f_a;  my@f_X; my@f_Y;   # yuck, give 'em names... 
$outdir.=qq{/$n}; # print qq{argh= $outdir\n};
$ARGV[1]||die qq{usage: $0 <name> <bed_file> [force]\n}; my$fn=$ARGV[1]; my$ref_should_be_fai='/scratch/HS_Build37/BWA_INDEX_hs37d5/hs37d5.dict';
die qq{dir '$outdir' already exists} if(-d $outdir && (!defined($ARGV[3]) || $ARGV[3] ne 'force'));
mkdir($outdir)||die;
my$fb; { my@q=split(/\//,$fn); $fb=pop(@q); $fb=~s{\.bed$}{}; }
open(my$f,'<',$fn)||die; # my($bn,$dn)=(basename($fn),dirname($fn));
my@argh=qw/all autosome X Y/; my@o=map{$outdir.'/'.$fb.'_'.$_.'.bed'}@argh; my%fns;
{ for my$f (0..$#o){
    # print qq{using $f $argh[$f] $o[$f]\n};
    open($fns{$argh[$f]},'>',$o[$f])||die;
}
my%fai; @fai{map{chomp;my@b=split(qq{\t});substr($b[1],3)}`cat $ref_should_be_fai`}=(); # we have no fai file?!?
while(my$bored=<$f>){ chomp($bored); my@x=split("\t",$bored); next if(index($bored,"\t")<0);
    substr($x[0],0,3)=q{} if(substr($x[0],0,3) eq 'chr');
    if(!exists$fai{$x[0]}) {
        print qq{[$c] Fix the scaffold name, '$x[0]', is not a known region from '$ref_should_be_fai'};
        last if(++$c>10);;
    }
    # mainly just doing grep -P '\t' bedfile | awk  '{print $3-$2+1}' | awk '{s+=$1} END {print s}'
    # if($x[0] eq q{X}){      push(@f_X,join(@x[0..2])); $s_X+=$x[2]-$x[1]+1; 
    if($x[0] eq q{X}){      my$fh=$fns{X}; print $fh join("\t",@x[0..2]).qq{\n}; $s_X+=$x[2]-$x[1]+1; 
    }elsif($x[0] eq q{Y}){  my$fh=$fns{Y}; print $fh join("\t",@x[0..2]).qq{\n}; $s_Y+=$x[2]-$x[1]+1; 
    }else{ my$fh=$fns{autosome}; print $fh join("\t",@x[0..2]),qq{\n}; $s_a+=$x[2]-$x[1]+1; }
    my$fh=$fns{all}; print $fh join("\t",@x[0..2]).qq{\n}; 
} } # close files...
# replace into captureKit values
# ('TwistV1.3','TwistV1.3','all','/nfs/goldsteindata/refDB/captured_regions/Build37/TwistV1.3/Twist_Exome_Target_hg19_all.bed',3,'TwistV1.3
# all'
die qq{fix the names please\n} if($c);
print qq{len= $s_a, $s_X, $s_Y\n};

{
# all'
    for my$f (0..$#o){
    ##### yuck. this offends me...?!?
    print qq{using $f $argh[$f] $o[$f] : };
    my$S    =   $argh[$f] eq 'autosome'     ?       $s_a
            :   $argh[$f] eq 'all'          ?       $s_a+$s_X+$s_Y
            :   $argh[$f] eq 'X'            ?       $s_X
            :   $argh[$f] eq 'Y'            ?       $s_Y
            :                                       0;
    print qq{\n};
    print qq{WARNING: Region '$argh[$f]' for file '$o[$f]' is empty\n} if($S==0);
    my$sql=qq{replace into captureKit values ('$n','$n','$argh[$f]','$o[$f]',$S,'}
      .($argh[$f] eq 'autosome'?'Autosomal':$argh[$f] eq 'all'?'Complete':$o[$f].' Chromosome')
      .qq{ $n bed file')};
    # print qq{$sql\n};
    &ARGH::mq($sql);
} }

my$cmd=qq{/nfs/seqscratch09/dsth/BedPatch/Bins2Keep.sh $n $o[0] /nfs/seqscratch09/dsth/BedPatch/ConsolidatedBed_OrionPurge.bed /nfs/seqscratch09/dsth/BedPatch 1kb};
# print qq{using $cmd\n};
system($cmd) && die;
my$base=qq{${n}_${fb}_all_};
die if(!-e $base.'1kbBlocks.bed' || !-e $base.'1kbBlocksIds.txt');
system("mv ${base}1kbBlocks.bed /nfs/seqscratch09/dsth/BedPatch/1KB_ConsolidatedBed_OrionPurge/") && die;
system("mv ${base}1kbBlocksIds.txt /nfs/seqscratch09/dsth/BedPatch/1KB_ConsolidatedBed_OrionPurge/") && die;
# seems there's a new table!?!
# print qq{e.g. insert into blocks_by_kit values ('IDTERPv1mtDNA','/nfs/seqscratch09/dsth/BedPatch/1KB_ConsolidatedBed_OrionPurge/IDTERPv1mtDNA_IDTERPv1mtDNA_all_1kbBlocksIds.txt'\n};
&ARGH::mq("replace into blocks_by_kit values ('$n','/nfs/seqscratch09/dsth/BedPatch/1KB_ConsolidatedBed_OrionPurge/${base}1kbBlocksIds.txt')");
print qq{bye\n};

