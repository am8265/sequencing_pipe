
#Yield_Sum.sh    
#Takes list of sample IDs and seqtype, and outputs run summary information
IFS=$'\n'
INPUT_LIST=$1
SEQTYPE=$2
EXOMEKIT=$3
SAMPLES=`less $INPUT_LIST`

sdb() {
        mysql --defaults-group-suffix=sequencedb "$@"
}


#Deletes previous run summary information
cd /nfs/genotyping/release_scripts/
rm Yield.txt

if [ "$EXOMEKIT" = 'roche' ]; then
    EXOMEKIT="Roche"
fi

#Ensures SEQTYPE is viable input
if [ "$SEQTYPE" = 'Custom' ]; then
    SEQTYPE='Custom Capture'
    abbrevSeqType='cc'
    fastqCheckType='genome'
fi
if [ "$SEQTYPE" = 'custom' ]; then
    SEQTYPE='Custom Capture'
    abbrevSeqType='cc'
    fastqCheckType='genome'

fi
if [ "$SEQTYPE" = 'custom_capture' ]; then
    SEQTYPE='Custom Capture'
    abbrevSeqType='cc'
    fastqCheckType='genome'

fi

if [ "$SEQTYPE" = 'genome' ]; then
    SEQTYPE='Genome'
    abbrevSeqType='wg'
    fastqCheckType='genome'

fi
if [ "$SEQTYPE" = 'exome' ]; then
    SEQTYPE='Exome'
    abbrevSeqType='ex'
    fastqCheckType='exome'

fi
if [ "$SEQTYPE" = 'rnaseq' ]; then
    SEQTYPE='RNASeq'
    abbrevSeqType='rs'
    fastqCheckType='genome'

fi


if [ -z "$SEQTYPE" ]; then
    echo Please enter seqtype
    exit 1
fi

#Loops over all samples in Input list and finds most recent prepID for specified SeqType. Then uses that prepID and gathers run summary information for each sample.
for s in $SAMPLES ; do
	prepIDs=$(sdb -e "SELECT \
		DISTINCT p.prepID FROM SeqType s \
		JOIN prepT p on p.prepID=s.prepID \
		WHERE \
		s.DBID=(select DISTINCT DBID FROM prepT WHERE CHGVID='$s') AND \
		s.SeqType='$SEQTYPE' AND \
		p.failedPrep='0' \
		AND p.exomekit='$EXOMEKIT' \
		AND p.CHGVID = '$s' \
		ORDER BY p.prepDate DESC" -NB)

	echo  "select DISTINCT p.prepID FROM SeqType s JOIN prepT p on p.prepID=s.prepID WHERE s.DBID=(select DISTINCT DBID FROM prepT WHERE CHGVID='$s') AND s.SeqType='$SEQTYPE' AND p.failedPrep='0' AND exomekit='$EXOMEKIT' AND p.CHGVID='$s' ORDER BY p.prepDate DESC"
	#echo $s $prepIDs
	for prepID in $prepIDs ; do
		

		#Read Summary Portion
		FCIDs=`sdb -e "select FCIllumID FROM Flowcell f JOIN Lane l ON f.FCID=l.FCID JOIN prepT p ON l.prepID=p.prepID WHERE (FailR1 IS NULL or FailR2 IS NULL) AND l.prepID='$prepID' AND f.fail=0 GROUP BY f.fcid" -NB  | sort -u`
		for FCID in $FCIDs ; do

			seqsata=`sdb -e "select SeqsataLoc FROM Flowcell WHERE FCIllumID='$FCID'" -NB`
			#echo $seqsata,$prepID,$FCID

			seqsatalocPrefix=$(echo $seqsata | cut -c1)
			#echo $seqsatalocPrefix,$s
			if [ "$seqsatalocPrefix" = 's' ] ; then
				seqDataPath="/nfs/$seqsata/seqfinal/whole_genome"

			elif [ "$seqsatalocPrefix" = 'f' ] ; then
				seqDataPath="/nfs/$seqsata/$(echo $SEQTYPE | tr '[:lower:]' '[:upper:]')"
			elif [ "$seqsatalocPrefix" = 'i' ] ; then
				seqDataPath="/nfs/$seqsata/$(echo $SEQTYPE | tr '[:lower:]' '[:upper:]')"
			else
				echo "Check SeqsataLoc for flowcell $FCID"
				#echo $seqsata,$FCID,$s,$seqsatalocPrefix
				exit 1
			fi

			#echo -e "select f.FCillumID FCID,GROUP_CONCAT(l.LaneNum ORDER BY l.LaneNum SEPARATOR '') LaneNum,CONCAT('$seqDataPath/',s.CHGVID) Loc,IF(l.FailR1 IS NULL AND l.FailR2 IS NULL,'PE','SE') 'READ',f.LenR1,f.LenR2,s.GwasID,s.Topstrandfile,p.exomeKit from Flowcell f JOIN Lane l ON f.FCID=l.FCID JOIN SampleT s ON l.DBID=s.DBID JOIN prepT p ON l.prepID=p.prepID WHERE (FailR1 IS NULL or FailR2 IS NULL) AND l.prepID='$prepID' AND f.FCillumID='$FCID'"
			info=$(sdb -e "select \
				f.FCillumID FCID,GROUP_CONCAT(l.LaneNum ORDER BY l.LaneNum SEPARATOR '') LaneNum,\
				CONCAT('$seqDataPath/',p.CHGVID) Loc,\
				IF(l.FailR1 IS NULL AND l.FailR2 IS NULL,'PE','SE') 'READ'\
				,f.LenR1,f.LenR2,s.GwasID,s.Topstrandfile,p.exomeKit \
				FROM Flowcell f \
				JOIN Lane l ON f.FCID=l.FCID \
				JOIN SampleT s ON s.DBID=l.DBID
				JOIN prepT p ON l.prepID=p.prepID \
				WHERE (FailR1 IS NULL or FailR2 IS NULL) AND \
				l.prepID='$prepID' AND \
				f.FCillumID='$FCID'" -NB | \
			       	sort -u | \
			       	sed 's/\t/,/g' | \
			       	sed 's/ /_/g' \
			)

			for f in $info ; do

				#echo $f
				flowcell=`echo $f | cut -d, -f1`
				all_lane=`echo $f | cut -d, -f2`
				sourc=`echo $f | cut -d, -f3`
				reads=`echo $f | cut -d, -f4`
				len1=`echo $f | cut -d, -f5`
				len2=`echo $f | cut -d, -f6`
				GwasID=`echo $f | cut -d, -f7`
				Topstrand=`echo $f | cut -d, -f8`
				ExomePrepKit=`echo $f | cut -d, -f9`

				if [ "$ExomePrepKit" = 'Roche_SeqCap_EZ_V3' ]; then
				ExomePrepKit='nimble'
				fi
				#echo $s,$flowcell

				if [ "$GwasID" = 'N/A' -o -z "$GwasID" ]; then
				GwasID=na
				fi

				if [ "$Topstrand" = 'N/A' -o -z "$Topstrand" ]; then
				Topstrand=na
				fi

				if [ -z "$ExomePrepKit" ]; then
				ExomePrepKit=na
				fi

				echo $flowcell,$all_lane,$sourc,$reads,$len1,$len2,,,,$GwasID,$Topstrand,$ExomePrepKit | sed 's/,/\t/g' | tee -a Yield.txt
			done
		done

	done
done

summaryOutFile=/nfs/sva01/Summaries/$(date +"%y.%m.%d")_${abbrevSeqType}Run_Summary_LSRC_v1.8.txt
touch $summaryOutFile
cat Yield.txt >> $summaryOutFile
chmod 775 $summaryOutFile
#cat Yield.txt

echo "Checking summary file and fastq files..."

for sample in `cat $INPUT_LIST`; do echo $sample ; perl /nfs/goldstein/goldsteinlab/software/Bioinformatics_Pipeline/production/fastq_check_cassava18.pl -s $sample -t $fastqCheckType -r $summaryOutFile ; done

echo $summaryOutFile
echo "Done"
