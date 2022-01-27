##### this stuff is terrible. was written in a major hurry when certain people went AWOL. should use mysql api...
package GoButton;
use strict;
use Sys::Hostname;
use Data::Dumper;
use File::Basename;
# use GoButton;                                 
use Cwd;
use Exporter qw(import);
our @EXPORT_OK = qw(sample_age);

# deal with syncing?!?! deal with passwords etc..?!?
my($user,$host,$pass,$db);

my$mysql=qq{mysql -u${user} -p${pass} -h${host} ${db} -B -e };

sub mysql {
    my$mysql=qq{mysql -u${user} -p${pass} -h${host} ${db} -B -e };
    return $mysql;
}

sub arsv2 {
    my$q=shift;
    # print qq{USING $q\n};
    my$cmd=qq{$mysql "$q"};
    print qq{USING $cmd\n};
    my@list=`$mysql "$q"`;
    my@out;
    my@h;
    for my$p (@list){ 
        chomp($p);
        my@p=split("\t",$p);
        if(scalar(@h)==0) {
            @h=@p;
        }else{
            my%h;
            for my$i (0..$#p){
                # print Dumper \@p;
                $h{$h[$i]}=$p[$i];
                # printf(qq{%-15s %-50s\n},$h[$i],$p[$i]);
                # &go_through_plate($p[$i]);
            }
            push(@out,+{%h});
        }
    }
# my@list=`$mysql "select * from plates where plateDate > 1510442653 order by plateid desc"`;
    return @out;
}

sub mq {
    my$q=shift;
    my$mysql=&mysql;
    # print qq{USING '$mysql $q'\n};
    my@list=`$mysql "$q"`;
    my@out;
    my@h;
    for my$p (@list){ 
        chomp($p);
        my@p=split("\t",$p);
        if(scalar(@h)==0) {
            @h=@p;
        }else{
            my%h;
            for my$i (0..$#p){
                # print Dumper \@p;
                $h{$h[$i]}=$p[$i];
                # printf(qq{%-15s %-50s\n},$h[$i],$p[$i]);
                # &go_through_plate($p[$i]);
            }
            push(@out,+{%h});
        }
    }
# my@list=`$mysql "select * from plates where plateDate > 1510442653 order by plateid desc"`;
    return @out;
}

sub sample_age {
    my$x=shift; $x||=0;
    my@x=caller(0);
    my$p=qq{$x[3]\t:\t};
    my$wd='/nfs/seqscratch09/informatics/tmp/WalDB.sample.sql'; # my$wd=$ENV{HOME}.'/WalDB.sample.sql';
    # my$wd='/nfs/seqscratch_ssd/informatics/tmp/WalDB.sample.sql'; # my$wd=$ENV{HOME}.'/WalDB.sample.sql';
    # why?!?
    chomp(my$wdf=`stat -c %Y ${wd} 2>/dev/null`);
    $wdf||=0;
    chomp(my$now=`date +%s`);
    my$age=sprintf("%.2f",($now-$wdf)/(60*60));
    print STDERR qq{${p}hours=$age\n} if($x >=0);
    if($x==1||$x==-2||$age>3){ # ?
        die qq{set ATAV_HOST\n} if(!$ENV{ATAV_HOST});
        my$ahost=$ENV{ATAV_HOST};
        die qq{set ATAV_USER\n} if(!$ENV{ATAV_USER});
        my$auser=$ENV{ATAV_USER};
        die qq{set ATAV_PASS\n} if(!$ENV{ATAV_PASS});
        my$apass=$ENV{ATAV_PASS};
        print qq{${p}dumping samples from WalDB\n} if($x>=-1);
        my$cmd=qq{mysqldump --skip-lock-tables --skip-triggers -u$auser -p$apass -h$ahost WalDB sample > }.$wd;
        # print qq{$cmd\n};
        system($cmd) && die;
        print qq{${p}loading samples from WalDB\n} if($x>=-1);
        $cmd=qq{mysql -u${user} -p${pass} -h${host} dsth < }.$wd;
        system($cmd) && die;
    }else{ 
        print STDERR qq{${p}NOT SYNCHING DB AS DUMP IS ONLY $age DAYS OLD\n} if($x>=0); 
    }
}

sub get_by_merged_status {
    my$x=shift;
    my$f=shift;
    $x||die;
    $f||='desc';
    my@wtf21=&mq("select d.pseudo_prepid from dragen_sample_metadata d where $x order by pseudo_prepid $f");
    # my@wtf21=&mq("select d.pseudo_prepid from dragen_sample_metadata d where $x order by priority desc");
    my@x=caller(0);
    my$p=qq{$x[3]\t:\t};
    print qq{${p}got }.scalar(@wtf21).qq{\n};
    return \@wtf21;
}

sub get_pipe_info {
    my$x=shift;
    $x||die;
    # print qq{getting pipe info for $x\n};
    my@wtf21=&mq("select p.*,d.step_name from dragen_pipeline_step p, dragen_pipeline_step_desc d where p.pipeline_step_id=d.id and pseudo_prepid = $x");#  order by");
    my%y;
    my%z;
    for my$i(@wtf21){
        $y{$i->{pipeline_step_id}}=$i;
        $z{$i->{step_status}}++;
    }
    my@x=caller(0);
    my$p=qq{$x[3]\t:\t};
    # print qq{${p}got }.scalar(@wtf21).qq{\n};
    return (\%y,\%z);
}

sub mu {
    my$q=shift;
    $q.=q{; select row_count() as affected;};
    my@arsv2=&mq($q);
    # print Dumper \@arsv2;
    die if(scalar(@arsv2)!=1); 
    die Dumper caller() if($arsv2[0]{affected}!=1);
    # die Dumper caller(0) if($arsv2[0]{affected}!=1);
}

sub get_preps_by_expt{
    my$eid=shift;
    my@x=&mq("select "
      .' s.sample_aka s_sample_aka, s.species s_species, e.subproject s_project, s.BroadPhenotype s_phenotype, p.prepid p_PREPID_NOT_PPID, p.sample_internal_name p_sample_name, '
      .' p.p_prepid p_pseudo_prepid, p.failedprep p_failedprep, dsm.is_external dsm_is_external, dsm.end_point dsm_end_point, ' #### ONLY USE dsm_end_point!?!
      .' p.sample_type p_sample_type, p.exomeKit p_exomeKit, p.experiment_id p_experiment_id, p.status p_status, p.status_time p_status_time, p.poolid p_poolid, p.externalData p_externalData, p.end_point p_end_point, '
      .' dsm.sample_name dsm_sample_name, dsm.sample_type dsm_sample_type, dsm.capture_kit dsm_capture_kit, dsm.pseudo_prepid dsm_pseudo_prepid, '
      .' dsm.priority dsm_priority, dsm.seqscratch_drive dsm_seqscratch_drive, dsm.is_merged dsm_is_merged, dsm.component_bams dsm_component_bams, '
      .' w.sample_name w_sample_name, w.sample_finished w_sample_finished, w.sample_failure w_sample_failure, w.prep_id w_pseudo_prepid, '
      .' w.sample_type w_sample_type, w.sample_id w_sample_id, '
      .' qm.pseudo_prepid qc_pseudo_prepid, qm.AlignSeqFileLoc qc_AlignSeqFileLoc '
      .' from prepT p '
      .' left join Experiment ex on p.experiment_id=ex.id '
      .' left join SampleT s on ex.sample_id=s.sample_id '
      .' left join dragen_sample_metadata dsm on p.p_prepid=dsm.pseudo_prepid '
      .' left join dragen_qc_metrics qm on p.p_prepid=qm.pseudo_prepid '
      .' left join dsth.sample w on p.p_prepid=w.prep_id '
      ." where p.experiment_id = ".$eid);
   
    # die qq{problem with $PP = $pp\n}.Dumper \@x if(scalar(@x)!=1);
    return \@x;
}

sub get_experiment_info{
    my$eid=shift;
    $eid||die;
    # print qq{pulling experiment $eid\n};
    my@x=&mq("select "
      .' s.sample_aka s_sample_aka, s.species s_species, e.subproject s_project, s.BroadPhenotype s_phenotype, '
      .' dsm.err_step, dsm.is_external dsm_is_external, dsm.end_point dsm_end_point, ' #### ONLY USE dsm_end_point!?!
      .' e.sample_type e_sample_type, e.exomeKit e_exomeKit, e.id e_experiment_id, '
      .' dsm.sample_name dsm_sample_name, dsm.sample_type dsm_sample_type, dsm.capture_kit dsm_capture_kit, dsm.pseudo_prepid dsm_pseudo_prepid, '
      .' dsm.priority dsm_priority, dsm.seqscratch_drive dsm_seqscratch_drive, dsm.is_merged dsm_is_merged, dsm.component_bams dsm_component_bams, '
      .' w.sample_name w_sample_name, w.sample_finished w_sample_finished, w.sample_failure w_sample_failure, w.prep_id w_pseudo_prepid, '
      .' w.sample_type w_sample_type, w.sample_id w_sample_id, '
      .' qm.pseudo_prepid qc_pseudo_prepid, qm.AlignSeqFileLoc qc_AlignSeqFileLoc '
      .' from Experiment e '
      # join prepT p on e.id=p.experiment_id '
      .' left join SampleT s on e.sample_id=s.sample_id '
      .' left join dragen_sample_metadata dsm on e.id=dsm.experiment_id'
      .' left join dragen_qc_metrics qm on e.id=qm.pseudo_prepid '
      .' left join dsth.sample w on e.id=w.prep_id '
      ." where e.id = ".$eid);
  # print Dumper \@x;
  # print qq{blarp\n};
    # next if($x[0]{p_sample_type} eq 'RNAseq');
    # die qq{this is a nonsense location\n}.Dumper \@x  if( $x[0]{dsm_seqscratch_drive} eq 'NULL' && $x[0]{p_status} ne 'External Data Submitted');
   
    if(scalar(@x)!=1) {
        die qq{problem with experiment_id:$eid =\n}.Dumper \@x, [caller()];
    }

    # print qq{okay\n};
    # print Dumper $x[0];
    # print qq{done\n};
    return $x[0];
}

sub get_info_by_pp {
    my$pp=shift;
    my$p=shift;
    $p||=0;
    $pp||die;
    my$PP=$p?'p.prepid':'p.p_prepid';
    my@x=&mq("select "
      .' s.sample_aka s_sample_aka, s.species s_species, e.subproject s_project, s.BroadPhenotype s_phenotype, p.prepid p_PREPID_NOT_PPID, p.sample_internal_name p_sample_name, '
      .' p.p_prepid p_pseudo_prepid, p.failedprep p_failedprep, dsm.is_external dsm_is_external, dsm.end_point dsm_end_point, ' #### ONLY USE dsm_end_point!?!
      .' p.sample_type p_sample_type, p.exomeKit p_exomeKit, p.experiment_id p_experiment_id, p.status p_status, p.status_time p_status_time, p.poolid p_poolid, p.externalData p_externalData, p.end_point p_end_point, '
      .' dsm.sample_name dsm_sample_name, dsm.sample_type dsm_sample_type, dsm.capture_kit dsm_capture_kit, dsm.pseudo_prepid dsm_pseudo_prepid, '
      .' dsm.priority dsm_priority, dsm.seqscratch_drive dsm_seqscratch_drive, dsm.is_merged dsm_is_merged, dsm.component_bams dsm_component_bams, '
      .' w.sample_name w_sample_name, w.sample_finished w_sample_finished, w.sample_failure w_sample_failure, w.prep_id w_pseudo_prepid, '
      .' w.sample_type w_sample_type, w.sample_id w_sample_id, '
      .' qm.pseudo_prepid qc_pseudo_prepid, qm.AlignSeqFileLoc qc_AlignSeqFileLoc '
      .' from prepT p '
      # .' left join ${db}.dragen_sample_metadata dsm on p.p_prepid=dsm.pseudo_prepid '
      # .' from ${db}.dragen_sample_metadata dsm join prepT p on dsm.pseudo_prepid=p.p_prepid left join dsth.sample w on dsm.pseudo_prepid=w.prep_id'
      .' left join Experiment ex on p.experiment_id=ex.id '
      .' left join SampleT s on ex.sample_id=s.sample_id '
      .' left join dragen_sample_metadata dsm on p.p_prepid=dsm.pseudo_prepid '
      .' left join dragen_qc_metrics qm on p.p_prepid=qm.pseudo_prepid '
      .' left join dsth.sample w on p.p_prepid=w.prep_id '
      ." where $PP = ".$pp.' group by p.experiment_id');
  # print Dumper \@x;
  # print qq{blarp\n};
    # next if($x[0]{p_sample_type} eq 'RNAseq');
    # die qq{this is a nonsense location\n}.Dumper \@x  if( $x[0]{dsm_seqscratch_drive} eq 'NULL' && $x[0]{p_status} ne 'External Data Submitted');
   
    if(scalar(@x)!=1) {
        print Dumper \@x;
        die qq{problem with $PP = $pp\n}.Dumper \@x, Dumper caller();
    }

    # print qq{okay\n};
    # print Dumper $x[0];
    # print qq{done\n};
    return $x[0];
}

sub get_scratch {
    my$f=shift;
    $f||die;
    return '/nfs/'.$f->{dsm_seqscratch_drive}.'/ALIGNMENT/BUILD37/DRAGEN/'.uc($f->{dsm_sample_type}).'/'.$f->{dsm_sample_name}.'.'.$f->{dsm_pseudo_prepid}.'/';
}

sub get_archive {
    my$f=shift;
    $f||die;
    return $f->{qc_AlignSeqFileLoc}.'/'.$f->{dsm_sample_name}.'.'.$f->{dsm_pseudo_prepid}.'/';
}
    # my@dps=&mq("select * from dragen_pipeline_step,dragen_pipeline_step_desc where pipeline_step_id=id "
    # ." and pipeline_step_id!= 1 and pseudo_prepid = ".$f->{pseudo_prepid});
    # my$s_dir='/nfs/'.$f->{dsm_seqscratch_drive}.'/ALIGNMENT/BUILD37/DRAGEN/'.uc($f->{dsm_sample_type}).'/'.$f->{dsm_sample_name}.'.'.$f->{dsm_pseudo_prepid}.'/';
    # my$dir=$s_dir;
    # my$a_dir=$f->{qc_AlignSeqFileLoc}.'/'.$f->{dsm_sample_name}.'.'.$f->{dsm_pseudo_prepid}.'/';
    # print qq{s_dir= $s_dir\na_dir= $a_dir\n\n};

    # { my$sge_wrapper=qq{$s_dir/sge_wrapper.log};
    # if(-e $sge_wrapper) {
        # print qq{this is done by newer procedure : $sge_wrapper\n};
        # next;
    # } }

###### not state aware - why pull the wrapper jobs if using this?!?
sub is_running {
    my$wjid=shift;
    my@qs=`qstat -j $wjid 2>&1`;                                                                                                                                         
    print("using $wjid\n");
    chomp($qs[0]);
    if(!defined $qs[0]) { die print qq{run me on a submission node!?!?\n};
    }elsif($qs[0] ne 'Following jobs do not exist: ') {
        print qq{THIS SEEMS TO BE RUNNING\n};
        print Dumper \@qs;
        return 1;
    }else{  
        print qq{didn't find job ($wjid) - should prolly check qacct for when it finished\n\n};
        return 0;
    }          
}

sub get_running {
    my%js;
    # ignore deleted and screwed jobs that just persist forever!?!?
    my@l=map{chomp;$_}`qstat -r | perl -pe 's{(\\s+\\d+\\s+)\\n}{\$1}' | grep -v -P '\\sdr\\s'|  grep Full | perl -ne 'm{(gatk_wrapper\\.\\S+)} && print qq{\$1\\n}'`;
    # my@l=map{chomp;$_}`qstat -r | perl -pe 's{(\\s+\\d+\\s+)\\n}{\$1}' | grep -v -P '\\s(Eqw|dr)\\s'|  grep Full | awk '{print \$12}'`;
    print Dumper \@l if(scalar(@l));
    # my@l=map{chomp;$_}`qstat -r | grep Full | grep -v 'merge' | awk '{print \$3}'`;
    my$wp=q{gatk_wrapper_};
    my$wp_new=q{gatk_wrapper.};
    for my$j (@l) {
        if(substr($j,0,length($wp)) eq $wp) { $js{ substr($j,length($wp),(length($j)-length($wp))) }{Wrapper}=$j;
        }elsif(substr($j,0,length($wp_new)) eq $wp_new) { $js{ substr($j,length($wp_new),(length($j)-length($wp_new))) }{Wrapper}=$j;
        }elsif($j=~m{^([A-Z]\w+)\.(\w+\.\d+)*$}) { $js{$2}{$1}=$j;
        }elsif($j eq 'ValidateBAM') { # argh...?!?
        }elsif($j eq 'STDIN') { # argh...?!?
        }else{ print qq{what is '$j'\n}; next; }
    }
    return \%js;
}

sub get_sdir { return q{/nfs/seqscratch09/informatics/logs/gatk/}; }
# sub get_sdir { return q{/nfs/seqscratch_ssd/informatics/logs/gatk/}; }

sub get_uname {
    my$x=shift;
    $x||die;
    return $x->{dsm_sample_name}.'.'.$x->{dsm_pseudo_prepid};
}

sub submit {

    my($dir,$sc,$pp,$state,$olds)=@_;

    die qq{there's no script file $sc} if(!-e $sc);
    die if(!$olds);

    print qq{updating to $state from $olds\n};

    &mu("update dragen_sample_metadata set is_merged = 2 where pseudo_prepid = ".$pp.qq{ and is_merged = $olds})
      if($olds==1);

    print "hmm: qsub -q wrapper.q -pe wrap 1 -l mem_free=20G $sc\n";
    chomp(my$jid=`qsub -q wrapper.q -pe wrap 1 -l mem_free=20G $sc`);

    my$X=qq{$jid\t$sc};
    my$jidn;
    if($X=~m{Your\sjob\s(\d+)\s}){ $jidn=$1; }
    else { 
        die qq{there's not job id available\n}; 
    }

    print qq{using $jidn / $X\n};
    system("echo 'ism:${olds}:$state\t$X' >> ${dir}/sge_wrapper.log") && die;

    return ($X,$jidn);
}

sub set_merged {
    my($pp,$s,$care)=@_;
    $care||='check';
    print qq{UPDATING TO $s\n};
    my@arsv2=&mq("update dragen_sample_metadata set is_merged = $s where pseudo_prepid = $pp ; select row_count() as affected;");
    print \@arsv2;
    die if(scalar(@arsv2)!=1);
    if($care ne 'ignore') { die Dumper \@arsv2 if($arsv2[0]{affected}!=1); }
}

sub relaunch {
    my($w,$dir,$sh,$e,$olds)=@_;
    unlink($w);
    print qq{submit $sh\n};
    my($X,$jid)=&GoButton::submit($dir,$sh,$e->{dsm_pseudo_prepid},3,$olds);
    print qq{got : $X/$jid\n};
    system("echo '$X' >> ${dir}/sge_wrapper.log") && die;
    print qq{jid = $jid\n};
}

sub reset {

    my($F,$who,$sge)=@_;
    unlink($who) if(-e$who);
    # unlink($sge) if(-e$sge);
    my@dpsscrewed=&mq("delete from dragen_pipeline_step where step_status ='failed' and pseudo_prepid = ".$F->{dsm_pseudo_prepid});
    # my@dpsscrewed=&mq("delete from dragen_pipeline_step where pipeline_step_id > 1 and pipeline_step_id < 32 and pseudo_prepid = ".$F->{dsm_pseudo_prepid});
    # my@dpsscrewed=&mq("delete from dragen_qc_metrics where pseudo_prepid = ".$F->{dsm_pseudo_prepid});
    my@dpsscrewed=&mq("update prepT set status = 'Released_to_Pipeline' where p_prepid = ".$F->{dsm_pseudo_prepid});
    my@arsv2=&mq("update dragen_sample_metadata set is_merged = 2 where pseudo_prepid = ".$F->{dsm_pseudo_prepid}.'; select row_count() as affected;');
    die if(scalar(@arsv2)!=1);
    die if($arsv2[0]{affected}!=1);

}

sub entry {                                                                                                                                                                    
    no strict 'refs';
    my@list; 
    ($host,$db,$user,$pass,@list)=@_;
    for my$x (@list) {
        &${x}; 
    }
}    

sub generate_sh {

    my$info=shift;
    my$local=shift;

    my$u_n=&GoButton::get_uname($info);
    my$dir=get_scratch($info);
    die qq{problem with unique name wrt., dir : $dir / $u_n\n} if(index($dir,$u_n)<0);
    # #\$ -pe wrap 1
    my$sdir=q{};
    my$pp=$info->{dsm_pseudo_prepid};
    
    if($local==0) {
        $sdir=&GoButton::get_sdir();
        $local=q{};
    }elsif($local==1){
        $sdir=$dir;
        $local=' --run-locally ';
    }else{ die qq{wtf?!?\n}; }

#\$ -p 199                                  
#\$ -M ${user}\@cumc.columbia.edu            

my$sh=qq{#!/bin/bash                               
#\$ -S /bin/bash                           
#\$ -cwd                                   
#\$ -V                                     
#\$ -l mem_free=10G                        
#\$ -o ${sdir}${u_n}.out
#\$ -e ${sdir}${u_n}.err
#\$ -N gatk_wrapper.${u_n}
export LUIGI_CONFIG_PATH=/nfs/goldstein/software/dragen_pipe/master/dragen/python/luigi.cfg 
export PATH=/nfs/goldstein/software/git-2.5.0/bin:/nfs/goldstein/software/python2.7.7/bin/:\$PATH; 
export PYTHONPATH=/nfs/goldstein/software/dragen_pipe/master/dragen/python 
/nfs/goldstein/software/python2.7.7/bin/luigi --module gatk_pipe ArchiveSample --logging-conf-file /nfs/goldstein/software/dragen_pipe/master/dragen/python/logging.conf --pseudo-prepid $pp --workers 3 --local-scheduler $local};

    my$sc=qq{$sdir/${u_n}.sh};
    print qq{in $sc\n};
    die qq{the scractch location on this is f'd : $sc\n} if(index($sc,'NULL')>=0);
    open(my$fh,'>',$sc)|| die qq{unable to write to script file $sc\n}.Dumper caller();
    print $fh $sh;
    close($fh)||die;
    die qq{this is missing $sc\n} if(!-e $sc);
    system("chmod +x $sc") && die qq{couldn't change permisions on $sc\n};
    print qq{come on $sc\n};
    return $sc;

}

sub most_recent_job {
    my$sge=shift;
    chomp(my$sgem=`cat $sge | grep Your | tail -n1`);
    my$jid=0;
    if($sgem=~/Your\s+job\s+(\d+)/sm) {
        $jid=$1;
        print qq{found jid = $jid\n};
    }else{
        die qq{WTF, unable to find jid from '$sge'\n\nwith='$sgem'};
    }
    return $jid;
}

sub push {

    my$PURGE=1;
    my$which=$ARGV[0]; $which||=q{};

    chomp(my$cur=`qstat | grep gatk -c`);
    print qq{there are currenty $cur\n};

    my($sample,$total,$max_p) = (10000,45000,1000);
    if($cur>$total){ print qq{won't launch more for now ($total)\n}; exit(0); }

    &GoButton::sample_age();

    my$ISM=1;
    my@wtf21=&GoButton::mq('select experiment_id FROM dragen_sample_metadata m ' # .q{ or end_point = 'picu' } # .q{ and (end_point = 'research' or end_point = 'clinical' }
    .q{ where is_merged = }.$ISM.' order by priority asc');

    my$sc=scalar(@wtf21);
    print qq{Checking $sc 'eligible' bams\n};

    for my$S(@wtf21){ # while(1) {

        my$pp=$S->{experiment_id};
        my@pipe=&GoButton::get_pipe_info($pp);
        my$info=&GoButton::get_experiment_info($pp);
        if($info->{err_step}>0) {
            print qq{WTF : we shouldn't ever have this!\n};
            next;
        }

        # for my$k(sort keys%{$info}) { printf("%-20s\t%s\n",$k,$info->{$k}); }

if(0) {
    print qq{running checks to avoid re-release\n};
        if($info->{qc_AlignSeqFileLoc} ne 'NULL') { print qq{This has an archive location - clean it up (exiting)!?!\n}; exit(1); }
        my$pc=scalar(keys%{$pipe[0]});
        # print qq{Step has $pc dps entries\n};
        if($pc==0){ print qq{This has no dps entries - it is 'pre-release'\n}; &GoButton::set_merged($pp,0); next; }
}else{
    print qq{not running checks\n};
}
        my$archiveonwards=0;
        for my$k(keys%{$pipe[0]}) { ++$archiveonwards if($k>30); }

        if(keys%{$pipe[0]}>1) { 
            die if(!$PURGE); 
        }

        my$screwed=exists$pipe[1]{failed}?1:0;
        my$released=exists$pipe[0]{1} && $pipe[0]{1}{step_status} eq 'completed'?1:0;

        if($released==1 && $screwed>=1) { 
            die if(!$PURGE); 
        }elsif($info->{w_sample_finished} ne 'NULL') { 
            print Dumper $info; 
            if( ( $info->{w_sample_finished}==0 && $info->{w_sample_failure}==0 )
                # || ($info->{w_sample_finished}==0 && $info->{w_sample_failure}==1)
            ) {
                print qq{we'll let this slide for now!?!?\n};
            }else{
                die qq{we need to be much more stringent about this\n};
            }

        }elsif($info->{w_sample_failure} ne 'NULL') { die; }

        my$dir=&GoButton::get_scratch($info);
        if(!-d $dir) { 
            # &GoButton::mu("update dragen_sample_metadata set is_merged = 80000 where experiment_id = ".$pp); next;
            die;
        }elsif(-e $dir.'/who.txt') { 
            print q{have a lock file }.$dir.'/who.txt - this should not happen - it was already pushed in?!?'; 
            &GoButton::mu("update dragen_sample_metadata set is_merged = 2 where experiment_id = ".$pp); 
            next;
            die $dir.'/who.txt'; 
        }

        if($released!=1) {
            die; # print Dumper $info;
        }

        die if($screwed>0 && !$PURGE); 
        die if($archiveonwards>0); 
        die if($pp!=$info->{dsm_pseudo_prepid}); 

        # print qq{dir= $dir\n};
        my$sc=&GoButton::generate_sh($info,0);
        my$sc=&GoButton::submit($dir,$sc,$pp,2,1);

        my@s=&GoButton::mq("update prepT set status = 'Pipeline Dispatched' where (failedprep = 0 or failedprep >= 100) and experiment_id = ".$pp."; select row_count() affected");
        print Dumper \@s;
        if($s[0]{affected}==0) {
            die "update prepT set status = 'Pipeline Dispatched' where (failedprep = 0 or failedprep >= 100) and experiment_id = ".$pp."; select row_count() affected";
        }

        if(--$max_p<=0) { print qq{reach max...\n}; exit(0); }

    }

    print qq{\nbye\n\n};
}

sub run_stuff {

    my$lim=1500;

    my$pp=shift;
    my$js=shift;
    my$F=shift;
    $F||=q{};
    my($d,$f,$messd,$mess)=(0,0,0,q{});
    my$e;

    eval { $e=&GoButton::get_experiment_info($pp); };
    if($@) {
        print qq{This is REALLY odd $F (experiment_id=$pp)\n};
        exit(1);
    }
    my$ism=$e->{dsm_is_merged};
   
    return if($e->{dsm_end_point} eq 'locked');
    # return if($e->{dsm_end_point} ne 'research' && $e->{dsm_end_point} ne 'clinical');
    return if( ( $ism < 2 || $ism > 9 ) && ( $ism < 90000 || $ism > 91000 )); # this in theory this is 'in pipe'
    return if($e->{err_step}>0); # don't flog the horse!?!
    # return if($e->{dsm_pseudo_prepid}>=111000&&$e->{dsm_pseudo_prepid}<=111100);

    # print qq{$_\t}.$e->{$_}.qq{\n} for (sort keys %{$e});
    if($ism>=10 && $ism<=40) { print qq{not processing ism=$ism\n}; die; } # paranoia

    my@a=&GoButton::get_pipe_info($pp);

    # print qq{\n-----------------------------------------------\n\n};

    # this was relaunching jobs that have completely finished and setting them to 2/3!?!
    if( $a[1]{completed}==27 && $a[0]{27}{step_status} eq 'completed' 
    # && $e->{p_status}=~m{Pipeline Completed ArchiveSample } # simplest to ignore for now?!?
    ) {
        print qq{hmm\n};
        &GoButton::set_merged($pp,10);
        next;
    }

    my$uname=&GoButton::get_uname($e);
    print qq{checking for $uname\n};

    if(exists$js->{$uname}) {
        print qq{\n\t> HAVE JOBS RUNNING SO IGNORE ($uname)\n\n\t};
        print Dumper $js->{$uname};
        print qq{\ndsm_is_merged= }.$e->{dsm_is_merged}.qq{\n};
        return;
    }else{ print qq{no jobs found for $uname\n}; }

    my$dir=&GoButton::get_scratch($e);

    my$w=qq{$dir/who.txt};
    my$sge=qq{$dir/sge_wrapper.log};

    if(!-e$dir) { 
        print Dumper $e;
        print qq{THIS IS SERIOUSLY WRONG?!?\n};
        my$adir=&GoButton::get_archive($e);
        print qq{checking $adir\n};
        if(-e qq{$adir/md5sums.txt}) {
            print qq{how is this happening!?!?\n};
            # sleep(1);
            &GoButton::mu(q{update dragen_sample_metadata set is_merged = 30 where experiment_id = }.$e->{dsm_pseudo_prepid});
        }else{
            die qq{Hmm: - missing scratch dir!?!\n}; 
        }
    }else{
        if(!-e$w) { 
            print Dumper $e;
            print qq{Hmm: - missing lock file ($w) - are we ignoring queued jobs?!?\n}; 
            # &GoButton::mu(q{update dragen_sample_metadata set is_merged = 1 where experiment_id = }.$e->{dsm_pseudo_prepid});
            sleep(3);
            # die;
        }
        if(!-e$sge) { 
            # &GoButton::mu(q{update dragen_sample_metadata set end_point = 'research',is_merged = 1 where experiment_id = }.$e->{dsm_pseudo_prepid});
            print qq{Hmm: - missing submission file ($sge)?\n}; 
            sleep(10);
            next;
            # die;
        }
    }

    chomp(my$wm=`cat $w`);
    chomp(my$sgem=`cat $sge`);

    print qq{checking $dir\n};
    print q{is_merged= }.$e->{dsm_is_merged}.qq{\n};

    chomp(my$wm=`cat $w`);

    my$jid=&GoButton::most_recent_job($sge);

    # &GoButton::reset($e,$w,$sge);

    my$sh=&GoButton::get_sdir.$uname.'.sh';

    print qq{\n----\nwho= '$wm'\n\n----\nsge= '$sgem'\n\n----\nsh= $sh\n\n----\n};

    # die qq{there's no script file $sh for $dir} if(!-e$sh);
    print qq{there's no script file $sh for $dir} if(!-e$sh);
    
    # paranoid double check of most recent jobs - i.e. we have a list of running jobs?!?
    if(!&GoButton::is_running($jid)) {
        print qq{wrapper jid ($jid) is missing\n};
        &GoButton::relaunch($w,$dir,$sh,$e,$ism);
    }else{ 
        print qq{wrapper jid ($jid) is running\n};
        return;
    }

    if($lim--<=0){ 
        print qq{done for now\n};
        exit(0);
    }

}

sub push_back {

    &GoButton::sample_age();

    # my@crash=map{chomp;$_}`grep java_command /nfs/goldstein/software/dragen_pipe/master/dragen/python/hs_err_pid*log 2>/dev/null | \
    # perl -pe 's{.*(picard\\.jar|\\-T) (\\S+).*}{\$2}' | \
    # head -300 | sort | uniq -c `;
    # what?!?
    if(0) { # scalar(@crash)){
        for my$f (map{chomp;$_}`ls /nfs/goldstein/software/dragen_pipe/master/dragen/python/hs_err_pid*log`) {
            #### could/should chnage thread number?!? must check job isn't running!?!? nah, let the scavanger on dragen/atav do that when it runs it?!?
            chomp(my$pp=`perl -ne 'm{(CUSTOM_CAPTURE|EXOME|GENOME_AS_FAKE_EXOME)\\/+([^\\.]+)\\.(\\d+)} && print qq{\$3\\t\$2\\t\$1\\t\$_\\n}' $f`);
            my@pp=split("\t",$pp);
            next if($pp[0]==111002 || $pp eq q{});
            next;
            my$e=&GoButton::get_info_by_pp($pp[0]); 
            if($pp[1] eq $e->{s_sample_name}) {
                die qq{fix this!?!\n};
            }
        }
        system('rm /nfs/goldstein/software/dragen_pipe/master/dragen/python/hs_err_pid*log') && die;
        sleep(10);
    }

    my$js=&GoButton::get_running();

    my$dbg=1;

    { print qq{using db\n}; sleep(1);
    my$c=&GoButton::get_by_merged_status('(is_merged >= 2 and is_merged <= 9) or (is_merged >= 90000 and is_merged <= 91000)');
    for my$d (@{$c}) { 
        &run_stuff($d->{pseudo_prepid},$js); 
    } }

    { print qq{using wrapper files\n}; sleep(1); 
    # for my$l (qw{/nfs/seqscratch_ssd/ALIGNMENT/BUILD37/DRAGEN /nfs/fastq_temp2/ALIGNMENT/BUILD37/DRAGEN /nfs/fastq_temp/ALIGNMENT/BUILD37/DRAGEN/}) {
    for my$l (qw{/nfs/seqscratch_ssd/ALIGNMENT/BUILD37/DRAGEN}){
        open(my$x,'-|',qq{find $l -name sge_wrapper.log -maxdepth 3})||die;
        while(my$d=<$x>){
            chomp($d);
            my$pp;
            # if($d=~/\.(\d+)\/sge_wrapper\.log/){ $pp=$1; }else{ die; }
            if($d=~/\.(\d+)\/sge_wrapper\.log/){ $pp=$1; }else{ last; }
            &run_stuff($pp,$js,qq{($d)});
        }
    } }
    print qq{\n\nbye\n};
}

sub archive_md5 {

    ++$|;
    my$s=$ARGV[0];
    my$timeout_d=10000;
    my@list=qw/ realn.recal.bam realn.recal.bai analysisReady.annotated.vcf.gz analysisReady.annotated.vcf.gz.tbi g.vcf.gz g.vcf.gz.tbi 
    coverage_bins pipeline_data.tar.gz /;
    my@list2=qw/coverage.tar.gz gq.tar.gz/;
    $s||=q{};
    my$ism=10;
    my@wtf21;

    # print '>>>>>>>>>>>>>>> '.$0.' : '.hostname().' : '.time().' : ';

    if($s) {
        @wtf21=&GoButton::mq("select m.*,q.AlignSeqFileLoc FROM dragen_sample_metadata m "
        ." left join dragen_qc_metrics q on m.pseudo_prepid=q.pseudo_prepid where m.pseudo_prepid = $s");
        $ism=$wtf21[0]{is_merged};
    }else{
        @wtf21=&GoButton::mq("select m.*,q.AlignSeqFileLoc FROM dragen_sample_metadata m, dragen_qc_metrics q "
            ." where m.pseudo_prepid=q.pseudo_prepid and is_merged = $ism and (m.pseudo_prepid < 111000 or m.pseudo_prepid >112000) limit 1");
    }

    if(scalar(@wtf21)==0) { 
        print qq{There's nothing to do\n};
        return;
        # exit(0);
    }elsif($wtf21[0]{AlignSeqFileLoc} eq 'NULL'){
        die qq{There's no archive location!\n};
    }

    print Dumper \@wtf21;

    my$f=$wtf21[0]; # print Dumper $f;

    my$fer="update dragen_sample_metadata set is_merged = 23 where pseudo_prepid = ".$f->{pseudo_prepid}.qq{ and is_merged = $ism};
    print qq{using $fer\n};
    &GoButton::mu($fer);

    print q{WE HAVE };
    my@md5s;
    print Dumper $f;

    my$uname=$f->{sample_name}.'.'.$f->{pseudo_prepid};
    my$sdir='/nfs/'.$f->{seqscratch_drive}.'/ALIGNMENT/BUILD37/DRAGEN/'.uc($f->{sample_type}).'/'.$f->{sample_name}.'.'.$f->{pseudo_prepid}.'/';
    my$adir=$f->{AlignSeqFileLoc}.'/'.$uname.'/';

    print qq{Attempting to check $sdir\n};
    die $sdir if(!-d $sdir);

    print qq{Attempting to check $adir\n};
    die $adir if(!-d $adir);

    print qq{Temporarily locking it sample (21)\n};
    &GoButton::mu("update dragen_sample_metadata set is_merged = 21 where is_merged = 23 and pseudo_prepid = ".$f->{pseudo_prepid});

    my$timeout=$timeout_d;
    if($f->{sample_type} eq 'Genome_As_Fake_Exome') {
        unshift(@list,'bam');
        $timeout=3*$timeout_d; # 15000;
    }

    for my$x (@list) {
        my$y=$adir.$uname.'.'.$x;
        push(@md5s,$uname.'.'.$x);
        print qq{checking $y\n};
        if(!-e $y) {
            my@arsv2=&GoButton::mu("update dragen_sample_metadata set is_merged = 2 where pseudo_prepid = ".$f->{pseudo_prepid});
            die qq{prob with $y};
        }
        die qq{prob with $y} if(!-s $y);
        system("chmod 440 $y") && die;
    }

    for my$x (@list2) {
        my$y=$adir.$x;
        # push(@md5s,$x);
        print qq{checking $y\n};
        die qq{prob with $y} if(!-e $y);
        die qq{prob with $y} if(!-s $y);
        system("chmod 440 $y") && die;
    }

    my$md5out='md5sums.txt';

    # just be paranoid and use abs paths!?!
    my$ARCHIVE_MD5SUMS=qq{$adir$md5out};

    if(-e $ARCHIVE_MD5SUMS) {
        print qq{seems that we're patching so remove '$ARCHIVE_MD5SUMS'\n}; 
        if(-e $ARCHIVE_MD5SUMS) {
            system("chmod 700 $ARCHIVE_MD5SUMS") && die;
            die if($ARCHIVE_MD5SUMS=~m{[\*\s]});
            my$ul=unlink($ARCHIVE_MD5SUMS);
            print qq{we unlink=$ul\n};
            die if($ul!=1);
        }
    }else{ print qq{standard run\n}; }

    # { chdir($adir)||die;
    # my$cmd=q{md5sum }.join(' ',map{$adir.$_}@md5s).qq{ > $ARCHIVE_MD5SUMS};
    # print qq{\nusing $cmd\n}; 
    # system($cmd) && die; 

    for my$i (0..$#md5s){
        my$cmd=qq{md5sum $adir/}.$md5s[$i].q{ }.($f==0?'>':'>>').q{ }.$ARCHIVE_MD5SUMS;
        print qq{\nusing $cmd\n}; 
        &system_timed($cmd,$timeout,$f);
    }

    my$s_md5file=qq{$sdir${md5out}.tmp};

    system("chmod 440 $ARCHIVE_MD5SUMS") && die;

    # { chdir($sdir)||die;
    my$cmd=q{md5sum }.join(' ',map{$sdir.$_}@md5s).qq{ > $s_md5file};
    print qq{\nusing $cmd\n}; 
    &system_timed($cmd,$timeout,$f);

    # my@md5so=`$cmd`;
    # { my$X=`md5sum $adir$md5out`;
    # my$Y=`md5sum $sdir$md5out`;

    { my$X=`awk '{print \$0}' $ARCHIVE_MD5SUMS | md5sum`;
    my$Y=`awk '{print \$0}' $s_md5file | md5sum`;
    if($X eq $Y) { die qq{this is strange!?!\n};
    }else{ print qq{this would be strange!?!\n}; } }

    { 
    my$X=`awk '{print \$1}' $ARCHIVE_MD5SUMS | md5sum`;
    my$Y=`awk '{print \$1}' $s_md5file | md5sum`;

    if($X ne $Y) { die qq{there's something nasty going on:\n$ARCHIVE_MD5SUMS -> $Y\n$s_md5file -> $X\n!?!\n};
    }else{ print qq{yay they match!?!\n}; } }

    print qq{push it through and cleanup $s_md5file\n};
    die if(unlink(qq{$s_md5file})!=1);

    &GoButton::mu("update dragen_sample_metadata set is_merged = 20 where pseudo_prepid = ".$f->{pseudo_prepid});

    print qq{\ndone\n};

    sub system_timed { # http://www.perlmonks.org/?node_id=90508
        my$call=shift;
        my$timer=shift;
        my$f=shift;
        my$sig=qq{timed_out\n};
        eval {  local $SIG{'ALRM'} = sub { die $sig; };
                alarm($timer); system("$call ; echo 'I made it'") && die "system call failed\n"; alarm(0); };
        if ($@) {
            print Dumper $f;
            my@arsv2=&GoButton::mu("update dragen_sample_metadata set is_merged = 22 where pseudo_prepid = ".$f->{pseudo_prepid});
            die if(scalar(@arsv2)!=1);
            die if($arsv2[0]{affected}!=1);
            print Dumper $f;
            if ($@ eq $sig) { print "timed out ($call)\n"; exit(1); }
            else { print "yikes, had a problem ($@)\n"; exit(1) }
        }
        else { print "made it ($call)\n"; }
    }

}


###### this was part of the back-patching of all old WGS to systematically generate full dragen gVCF - now it simply pulls back from the gVCF queue
sub run_gvcf {

    if(!$db) {
        ###### dind't do it this way in general to allow a better procedure upstream with a single set point... this is all pretty terrible
        # ($host, $user, $pass, $db)=@_ 
        die qq{Need LIMS_USER env\n} if(!$ENV{LIMS_USER});
        $user=$ENV{LIMS_USER};
        die qq{Need LIMS_PASS env\n} if(!$ENV{LIMS_PASS});
        $pass=$ENV{LIMS_PASS};
        die qq{Need LIMS_HOST env\n} if(!$ENV{LIMS_HOST});
        $host=$ENV{LIMS_HOST};
        $db='sequenceDB';
    }

    my$outdir=q{/nfs/informatics/production/gvcf};

    print qq{will run a single sample\n};

    # my@x=map{chomp;$_}`ls /nfs/seqscratch_ssd/dsth/staging/Batch_IV/*/*.bam`;
    my@x=&GoButton::mq("select d.*,WGS_Dragen_gVCF from dragen_sample_metadata d left join dragen_qc_metrics q on d.pseudo_prepid=q.pseudo_prepid "
    # ."where is_merged in (200000,200040,200100) order by is_merged asc");
    # get newest going first - i.e. get ipf to bypass quickly?!?
    # ."where is_merged = 200000 order by is_merged asc");
    ."where is_merged = 200000 order by experiment_id desc");
    
    print q{got }.scalar(@x).qq{ samples\n};

    my$delay=0;

    if($delay) {
        print qq{NOT RUNNING FOR NOW AS WILL ALLOW FOR DIRECT REMAPPING OF SAMPLES WITHOUT ARCHIVING!?!\n};
        exit(1);
    }

    # chomp(my$pwd=`pwd -P`);

    for my$x (@x) {

        print Dumper $x;

        if($x->{WGS_Dragen_gVCF} ne 'NULL') {
            die qq{this no longer makes no sense\n};
        }

        my$bam=q{};
        my$out=$outdir.q{/}.$x->{sample_name}.'.'.$x->{experiment_id}.'/';

        if($x->{is_merged}==200000) {
            print qq{this is new wgs sample : };
            $bam=q{/nfs/}.$x->{seqscratch_drive}.q{/ALIGNMENT/BUILD37/DRAGEN/}.uc($x->{sample_type}).'/'
            .$x->{sample_name}.'.'.$x->{experiment_id}.'/'
            .$x->{sample_name}.'.'.$x->{experiment_id}.'.bam';
        }else{
            die qq{this is no longer required\n};
            # print qq{this is catching up with old wgs samples : };
            # $bam=q{/nfs/fastq_temp2/homes/dsth/restaging_wgs/}
            # .$x->{sample_name}.'.'.$x->{experiment_id}.'/'
            # .$x->{sample_name}.'.'.$x->{experiment_id}.'.bam';
        }

        print qq{$bam\n};

        # this should really be in release procedure!?!

        if(0) { # if($x->{sample_name}=~/^(AZ)?IPF\d+$/ || $x->{sample_name}=~/^CGND/ || substr($x->{sample_name},0,3) eq 'Pul') {
                &GoButton::mu("update dragen_sample_metadata set is_merged = 1 where pseudo_prepid = ".$x->{experiment_id});
                &GoButton::mq("update prepT set status = 'Released_to_Pipeline', status_time = unix_timestamp() where experiment_id = ".$x->{experiment_id});
                &GoButton::mq("update dragen_pipeline_step set finish_time = CURRENT_TIMESTAMP(), step_status = 'completed' where pseudo_prepid = ".$x->{experiment_id}." and pipeline_step_id = 1");
                # &GoButton::mu("update prepT set status = 'Released_to_Pipeline', status_time = unix_timestamp() where experiment_id = ".$x->{experiment_id});
                # &GoButton::mu("update dragen_pipeline_step set finish_time = CURRENT_TIMESTAMP(), step_status = 'completed' where pseudo_prepid = ".$x->{experiment_id}." and pipeline_step_id = 1");
                next;
        }

        if(!-e$bam) {

            print qq{the bam ($bam) is missing!?!?\n};
            ###############################################################
            my$bism=$x->{is_merged};
            my$bism2=$bism+3;
            &GoButton::mu("update dragen_sample_metadata set is_merged = $bism2 where pseudo_prepid = ".$x->{experiment_id}." and is_merged = $bism ");
            ###############################################################
            # die qq{WTF?!?\n};
            next;

        }

        my$dir=dirname($bam);
        my$base=basename($bam,'.bam');
        print qq{> have $out, $dir & $base\n\n};

        my$log = qq{${out}/${base}.log};
        my$md5=$out.q{/}.$x->{sample_name}.'.'.$x->{experiment_id}.q{.gvcf.gz.md5sum};

        if(-e $md5) {
            my$gvcf=$md5;
            $gvcf=~s{\.md5sum}{};
            if (!-e $gvcf) {
                die qq{this makes no sense!\n};
            }else{
                print qq{let's do some cleanup ($gvcf, $md5)\n};
            }
            ##### don't do this automatically as we don't protect finished data
            # unlink($gvcf);
            # unlink($md5);
            die qq{this almost certainly means we deprecated a sample and re-released. this means we really must wipe all the files and re-run to avoid picking up old versions\n};
        }

        print qq{bam= $bam\nlog= $log\nout= $out\nmd5= $md5\n\n};
        if (!-d$out) {
            print qq{generate output dir $out\n};
            mkdir($out);
        }
        die qq{couldn't generate output dir\n} if (!-d$out);

        my$yay=qq{${out}/${base}.yay};
        my$fucv2=qq{${out}/${base}.boo};
        
        my$host=$ENV{HOSTNAME};
        print qq{host=$host\n};

        my$started=qq{${out}/${base}.wgs.lock};
        # my$started=qq{${out}/${base}.${host}.lock};

        print qq{started=$started\n};

        # exit(1);

        my$drgn_cmd=q{dragen --ref-dir /staging/REF/b37_decoy/ --dbsnp /nfs/seqscratch_ssd/PIPELINE_DATA/dbsnp_138.b37.vcf.gz }
        .q{--enable-variant-caller true --enable-vcf-indexing true --enable-vcf-compression true --enable-map-align false }
        .q{--vc-emit-ref-confidence GVCF --vc-max-alternate-alleles 6 }
        .qq{--output-directory $out }
        .qq{--bam-input $bam --vc-sample-name $base --output-file-prefix ${base} }
        .qq{--vc-gvcf-gq-bands 5 20 60 2>&1 | tee $log};

        print qq{using\n$drgn_cmd\n};

        if(-e$started) {
            print qq{this is going ($started)\n};
# no longer running with risk of races...
# next;
        }

        # print Dumper $x; 
        # next; 
        # exit(1);

        ###############################################################
        my$bism=$x->{is_merged};
        my$bism2=$bism+2;
        &GoButton::mu("update dragen_sample_metadata set is_merged = $bism2 where pseudo_prepid = ".$x->{experiment_id}." and is_merged = $bism ");
        system(qq{touch $started});
        ###############################################################

        if(system($drgn_cmd)==0){ system(qq{touch $yay});
        }else{ system(qq{touch $fucv2}); }

        if(-e $md5) {

            # dragen_align_se align, check_merge_add-some-tracking-info_and_push-to-gatk.cpp, dragen_pipe thing and this all need merging into single consolidated procedure!?!
            # then get someone else to run it?!?
                
            #### 'should' protect output!?!
            
            if($x->{is_merged}==200000) {

                print qq{release to pipeline\n};
                &GoButton::mu("update dragen_sample_metadata set is_merged = 1 where pseudo_prepid = ".$x->{experiment_id});
                # &GoButton::mq("update dragen_sample_metadata set is_merged = 1 where pseudo_prepid = ".$x->{experiment_id});
                &GoButton::mq("update prepT set status = 'Released_to_Pipeline', status_time = unix_timestamp() where experiment_id = ".$x->{experiment_id});
                # &GoButton::mq("update prepT set status = 'Released_to_Pipeline', status_time = CURRENT_TIMESTAMP() where experiment_id = ".$x->{experiment_id});
                &GoButton::mq("update dragen_pipeline_step set finish_time = CURRENT_TIMESTAMP(), step_status = 'completed' where pseudo_prepid = ".$x->{experiment_id}." and pipeline_step_id = 1");

                ############ need to update dqm directly - just yanking bits from add_in_gvcf_location.pl
                ############ BUT DARNED QC METRICS ENTRY ISN'T CREATE TILL MUCH LATER IN PIPELINE - hence WAS USING ADD IN THING - REVERT TO THAT OR USE INSERT OR JUST PUT INTO PIPELINE!?!
                # $md5=~s{\.md5sum}{};
                # print qq{need to update this one with $md5\n};
                # &GoButton::mu("update dragen_qc_metrics set WGS_Dragen_gVCF = '$md5' where pseudo_prepid = ".$x->{experiment_id});
                # &GoButton::mu("update dragen_qc_metrics set WGS_Dragen_gVCF = '$l' where pseudo_prepid = $exp");

            }else{
                die qq{this is no longer used!?!\n};
                # print qq{wipe bam\n};
                # unlink($bam);
                # &GoButton::mq("update dragen_sample_metadata set is_merged = ".($x->{is_merged}-200002)." where pseudo_prepid = ".$x->{experiment_id});
            }

        }else{ 
            die qq{there's a serious issue\n};
            # exit(1); 
        }

        ##### release dragen wrapper - this is now completely silly but not gonna change it this moment
        exit(0);

    }
}
# &archive_md5();

1;

