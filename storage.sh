#!/bin/bash
#store_seqsata.sh
#Joshua Bridgers
#Stores fastq.gz from Unaligned folder of a run.  Also moves summary files and creates log files in seqsata.

FCID=$1
source=$2
seqsata=/nfs/fastq15
runfolder=$3
email=jb3816@cumc.columbia.edu

completed=`grep -o "INFO: all completed" $source/*$FCID*/nohup.sge`
MACHINE=$(echo $runfolder | cut -d_ -f2)
DATA_DRIVE=$(echo $seqsata | cut -d/ -f3)
LOG_FILE=$seqsata/summary/GAF_PIPELINE_LOGS/${MACHINE}_${FCID}_${DATA_DRIVE}.log

#get original sum of fastq.gz file sizes
O_filesize=`du --apparent-size -c $source/*$FCID*/Project*/Sample*/*fastq.gz | tail -1 | cut -f1`

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

#Iterate over every sample on the flowcell
for s in $source/*$FCID*/Project*/Sample*; do

	sampleID=$(echo $s | awk -F/ '{print $NF}' | cut -d_ -f2-)
	seqtype=$(~/sequenceDB.sh "SELECT distinct(st.seqtype) FROM Lane l JOIN Flowcell f ON l.fcid=f.fcid JOIN SeqType st ON l.prepid=st.prepid JOIN prepT p ON l.prepid=p.prepid WHERE FCILLUMID='$FCID' AND CHGVID='$sampleID'" -NB | tr '[:lower:]' '[:upper:]' | sed 's/ /_/g')
	echo "~/sequenceDB.sh \"SELECT distinct(st.seqtype) FROM Lane l JOIN Flowcell f ON l.fcid=f.fcid JOIN SeqType st ON l.prepid=st.prepid JOIN prepT p ON l.prepid=p.prepid WHERE FCILLUMID='$FCID' AND CHGVID='$sampleID'\" -NB | tr '[:lower:]' '[:upper:]' | sed 's/ /_/g'" >> $LOG_FILE

	mkdir -p $seqsata/$seqtype/$sampleID/$FCID
	chmod 775 $seqsata/$seqtype/$sampleID
	mv $s/* $seqsata/$seqtype/$sampleID/$FCID
	rsync -avP $source/*$FCID*/Basecall_Stats_$FCID/Demultiplex_Stats.htm $seqsata/$seqtype/$sampleID/$FCID
	ls -al $seqsata/$seqtype/$sampleID/$FCID > $seqsata/$seqtype/$sampleID/$FCID/$sampleID.$FCID.files.txt
	
	echo mkdir -p $seqsata/$seqtype/$sampleID/$FCID >> $LOG_FILE
	echo chmod 775 $seqsata/$seqtype/$sampleID >> $LOG_FILE
	echo mv $s/* $seqsata/$seqtype/$sampleID/$FCID >> $LOG_FILE
	echo rsync -avP $source/*$FCID*/Basecall_Stats_$FCID/Demultiplex_Stats.htm $seqsata/$seqtype/$sampleID/$FCID >> $LOG_FILE
	echo ls -al $seqsata/$seqtype/$sampleID/$FCID \> $seqsata/$seqtype/$sampleID/$FCID/$sampleID.$FCID.files.txt >> $LOG_FILE

done

mv_filesize=`du --apparent-size -c $seqsata/*/*/$FCID/*fastq.gz | tail -1 | cut -f1` 

#Create SAV file of run
zip $runfolder/${FCID}_$(echo $runfolder | awk -F/ '{print $NF}' | cut -d_ -f1,2)_SAV.zip $runfolder/RunInfo.xml $runfolder/runParameters.xml $runfolder/InterOp/

#Save bcl2fastq output
zip $source/$FCID.bcl.nohup.zip $source/*$FCID*/nohup.sge

#Copy log files
mv $source/$FCID.bcl.nohup.zip $seqsata/summary/bcl_nohup
rsync -avP $runfolder/$FCID*_SAV.zip $seqsata/summary/SAV/

echo "SUM of fastq.gz files after move: $mv_filesize" >> $LOG_FILE
echo "================================================================================" >> $LOG_FILE

#check if filesizes are the same after the move
if [ $O_filesize != $mv_filesize ] 
then
	echo "BCL move FAILURE for $FCID on `date`" | mail -s "GAF:BCL move FAILURE" $email
	echo failure
else
	touch $runfolder/rsync_complete.txt
	#rm -rf $seqsata/*$FCID*
	#echo "Removing BCL Unaligned folder" >> $LOG_FILE
	#echo "rm -rf $seqsata/*$FCID*" >> $LOG_FILE
	echo "Done"
fi

touch $runfolder/StorageComplete.txt

