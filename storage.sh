getopt --test > /dev/null

if [[ $? != 4 ]]; then
	echo "I’m sorry, `getopt --test` failed in this environment."
	exit 1
fi

SHORT=b:i:s:
LONG=bcl:,input:,seqsata:

PARSED=`getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@"`
if [[ $? != 0 ]]; then
	exit 2
fi

eval set -- "$PARSED"

while true; do
    case "$1" in
        -b|--bcl)
            BCLDrive=$2
            shift 2
            ;;
        -i|--input)
            runFolder=$2
            shift 2
            ;;
        -s|--seqsata)
            Seqsata="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done


MACHINE=$(echo $runFolder | cut -d_ -f2)
DATA_DRIVE=${Seqsata}
FCID=$(echo $runFolder | cut -d_ -f4)

source=$BCLDrive/
seqsata=/nfs/${Seqsata}/

email="jb3816@cumc.columbia.edu dn2404@cumc.columbia.edu"
completed=`grep -o "INFO: all completed" $source/*$FCID*/nohup.sge`
LOG_FILE=$seqsata/summary/GAF_PIPELINE_LOGS/${MACHINE}_${FCID}_${DATA_DRIVE}.log

echo "input: $runFolder, bcl: $BCLDrive, seqsata: $Seqsata" >> $LOG_FILE

#get original sum of fastq.gz file sizes
O_filesize=`du --apparent-size -c $source/*$FCID*/Project*/Sample*/*fastq.gz | tail -1 | cut -f1`
O_fileNum=$(ls $source/*$FCID*/Project*/Sample*/*fastq.gz | wc -l)
#Check if BCL finished successfully
if [ -n "$completed" ] 
then
	echo "BCL completed successfully!"
else
	echo "BCL did NOT complete successfully!"
	exit 0
fi

cd $source
echo "Start of $FCID storage logging" >> $LOG_FILE
echo "================================================================================" >> $LOG_FILE

echo "SUM of fastq.gz files: $O_filesize" >> $LOG_FILE
echo "Number of fastq.gz files: $O_fileNum" >> $LOG_FILE

#Iterate over every sample on the flowcell
for s in $source/*$FCID*/Project*/Sample*; do

	sampleID=$(echo $s | awk -F/ '{print $NF}' | cut -d_ -f2-)
	seqtype=$(mysql -usequence_connect -h 10.73.50.38 -pg3n3t1c5213 sequenceDB -e "SELECT distinct(st.seqtype) FROM Lane l JOIN Flowcell f ON l.fcid=f.fcid JOIN SeqType st ON l.prepid=st.prepid JOIN prepT p ON l.prepid=p.prepid WHERE FCILLUMID='$FCID' AND CHGVID='$sampleID'" -NB | tr '[:lower:]' '[:upper:]' | sed 's/ /_/g')
	echo "mysql -usequence_connect -h 10.73.50.38 -pg3n3t1c5213 sequenceDB -e \"SELECT distinct(st.seqtype) FROM Lane l JOIN Flowcell f ON l.fcid=f.fcid JOIN SeqType st ON l.prepid=st.prepid JOIN prepT p ON l.prepid=p.prepid WHERE FCILLUMID='$FCID' AND CHGVID='$sampleID'\" -NB | tr '[:lower:]' '[:upper:]' | sed 's/ /_/g'" >> $LOG_FILE

	echo "starting transfer of $sampleID" >> $LOG_FILE
	echo "================================================================================" >> $LOG_FILE

	echo mkdir -p $seqsata/$seqtype/$sampleID/$FCID >> $LOG_FILE
	mkdir -p $seqsata/$seqtype/$sampleID/$FCID

	echo chmod 775 $seqsata/$seqtype/$sampleID >> $LOG_FILE
	chmod 775 $seqsata/$seqtype/$sampleID

	echo "Original files list" > $seqsata/$seqtype/$sampleID/$FCID/$sampleID.$FCID.files.txt
	echo ls -ltr $s* \>\>  $seqsata/$seqtype/$sampleID/$FCID/$sampleID.$FCID.files.txt >> $LOG_FILE
	ls -ltr $s* >>  $seqsata/$seqtype/$sampleID/$FCID/$sampleID.$FCID.files.txt	
 
	#add if statement 
	if [ $source == "/nfs/fastq18/BCL/" ] ; then
		echo $(date) mv $s/* $seqsata/$seqtype/$sampleID/$FCID >> $LOG_FILE
		mv $s/* $seqsata/$seqtype/$sampleID/$FCID
	else
		echo $(date) rsync -avP $s/* $seqsata/$seqtype/$sampleID/$FCID >> $LOG_FILE
		rsync -avP --remove-source-files $s/* $seqsata/$seqtype/$sampleID/$FCID
	fi

	echo rsync -avP $source/*$FCID*/Basecall_Stats_$FCID/Demultiplex_Stats.htm $seqsata/$seqtype/$sampleID/$FCID >> $LOG_FILE
	rsync -avP $source/*$FCID*/Basecall_Stats_$FCID/Demultiplex_Stats.htm $seqsata/$seqtype/$sampleID/$FCID

	echo "Transferred file list" >> $seqsata/$seqtype/$sampleID/$FCID/$sampleID.$FCID.files.txt	
	echo ls -ltr $seqsata/$seqtype/$sampleID/$FCID \>\> $seqsata/$seqtype/$sampleID/$FCID/$sampleID.$FCID.files.txt >> $LOG_FILE
	ls -al $seqsata/$seqtype/$sampleID/$FCID >> $seqsata/$seqtype/$sampleID/$FCID/$sampleID.$FCID.files.txt	

done

mv_filesize=`du --apparent-size -c $seqsata/*/*/$FCID/*fastq.gz | tail -1 | cut -f1`
echo du --apparent-size -c $seqsata/*/*/$FCID/*fastq.gz | tail -1 | cut -f1 >> $LOG_FILE

mv_filenum=`ls $seqsata/*/*/$FCID/*fastq.gz | wc -l` 
echo ls $seqsata/*/*/$FCID/*fastq.gz | wc -l >> $LOG_FILE

#Create SAV file of run
tar cvzf /nfs/seqscratch1/Runs/$runFolder/${FCID}_$(echo $runFolder | xargs -n1 basename | cut -d_ -f1,2)_SAV.tar.gz /nfs/seqscratch1/Runs/$runFolder/RunInfo.xml /nfs/seqscratch1/Runs/$runFolder/runParameters.xml /nfs/seqscratch1/Runs/$runFolder/InterOp/ >> $LOG_FILE

#Save bcl2fastq output
zip $source/$FCID.bcl.nohup.zip $source/*$FCID*/nohup.sge

#Copy log files
mv $source/$FCID.bcl.nohup.zip $seqsata/summary/bcl_nohup
rsync -avP /nfs/seqscratch1/Runs/$runFolder/$FCID*_SAV.tar.gz $seqsata/summary/SAV/

echo "SUM of fastq.gz files after move: $mv_filesize" >> $LOG_FILE
echo "Number of fastq.gz files after move: $mv_filenum" >> $LOG_FILE

echo "================================================================================" >> $LOG_FILE

#check if filesizes are the same after the move
if [ $O_filesize != $mv_filesize ] 
then
	echo "BCL move FAILURE for $FCID on `date`" | mail -s "GAF:BCL move FAILURE" $email
	echo failure
else
	touch /nfs/seqscratch1/Runs/$runFolder/rsync_complete.txt
	echo "BCL move SUCCESS for $FCID on `date`" | mail -s "GAF:BCL move SUCCESS" $email
	#rm -rf $seqsata/*$FCID*
	#echo "Removing BCL Unaligned folder" >> $LOG_FILE
	#echo "rm -rf $seqsata/*$FCID*" >> $LOG_FILE
	touch /nfs/seqscratch1/Runs/$runFolder/StorageComplete.txt
	echo "Done"
fi


