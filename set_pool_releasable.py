### Written by Nitin Bhardwaj, IGM
### This program lookf for an unrelead exome pool and if a certain number of samples in that pool meet coverage, the is_releasable is set to 1 for that pool
from __future__ import division

import MySQLdb as mdb
import os
import subprocess
import sys
import json
import smtplib,email,email.encoders,email.mime.text,email.mime.base
from email.mime.multipart import MIMEMultipart as MM

try: 
    if (sys.argv[1] == 'email'):
        print("running with email")
        Email = 1
    rejected_pools = []
except:
    print("running without email")
    Email = 0

def emailit(msg,body):

    smtpserver = 'localhost'
    fromAddr = 'igm-bioinfo@columbia.edu'
    emailMsg = MM('alternative')
    emailMsg['Subject'] = msg
    emailMsg['From'] = fromAddr
    to=['IGM-bioinfo@columbia.edu','igm-hts@columbia.edu'] #,'test@cumc.columbia.edu']
    emailMsg['To'] = ', '.join(to)
    body = body
    emailMsg.attach(email.mime.text.MIMEText(body,'html'))
    server = smtplib.SMTP(smtpserver)
    server.sendmail(fromAddr,to,emailMsg.as_string())
    server.quit()


STAGESEQ_USER = os.getenv('STAGESEQ_USER')
STAGESEQ_PASS = os.getenv('STAGESEQ_PASS')

PRODSEQ_USER = os.getenv('STAGESEQ_USER')
PRODSEQ_PASS = os.getenv('STAGESEQ_PASS')
try:
    #con = mdb.connect('devseq', STAGESEQ_USER, STAGESEQ_PASS, 'sequenceDB')
    con = mdb.connect('prodseq', PRODSEQ_USER, PRODSEQ_PASS, 'sequenceDB')
except:
    print("Could not connect to the database: check login/password") 

with con:

    qualifying = 0

    cur = con.cursor(mdb.cursors.DictCursor)
    cur.execute("select distinct p.id, p.name, p.count from pool p, prepT pt, Experiment e, Lane l, Flowcell f where p.id = pt.poolID and pt.experiment_id = e.id and pt.prepID = l.prepID and l.FCID = f.FCID and (pt.failedprep = 0 or pt.failedprep >= 100) and f.pipelinecomplete = 1 and l.failr1 is null and l.failr2 is null and pt.externaldata is null and e.sample_type = 'Exome' and e.is_released in ('not_released', 'release_rejected') and p.is_releasable=0")
    pools = cur.fetchall()
    for pool in pools:
        poolID = pool["id"]
        pool_name=pool["name"]
        total = pool["count"]
        #print(poolID, total)


        #print("checking to see if pool {0} can be released".format(poolID))
        # query = "select count(*) null_count from prepT pt, Lane l where pt.prepID = l.prepID and rg_metrics_capturemean is null and l.rg_metrics_status != 'low_yield' and pt.poolID = {0};".format(poolID)
        query = "select count(*) null_count from prepT pt, Lane l where pt.prepID = l.prepID and rg_metrics_capturemean is null and l.rg_metrics_status != 'low_yield' and l.rg_metrics_status != 'failed_flowcell' and  pt.poolID = {0};".format(poolID)
	cur.execute(query)
        result = cur.fetchone()
        null_count_and_not_low_yield = result["null_count"]



        if null_count_and_not_low_yield == 0:
           processed = 1
        else:
           processed = 0
           print("pool {0} still processing".format(poolID))

        #print ("processed status = {0}".format(processed))
        qualifying = 0
        if processed == 1:
            query = "select experiment_id, subproject_id from prepT pt, Experiment e where pt.experiment_id = e.id and pt.pool_status = 'pool_members' and pt.poolID = {0} and subproject_id is not NULL;".format(poolID)
            cur.execute(query)
            expt_ids = cur.fetchall()

            for exp in expt_ids:
               exp_id = exp["experiment_id"]
               subproject_id = exp["subproject_id"]
               curl_query = "curl -s -H 'Accept: application/json; indent=4' 'https://core.igm.cumc.columbia.edu/core/api/subproject/?id={0}'".format(subproject_id)
               result = os.popen(curl_query).read()
               y = json.loads(result)
               if y and y[0]["project"]["exome_target_min_coverage"]: # and y[0]["project"]["exome_target_min_coverage"] is not None:
                   wes_min_cov = y[0]["project"]["exome_target_min_coverage"]
               else:
                   wes_min_cov = 50
          
               #print(wes_min_cov) 

               query = "select pt.prepID,sum(rg_metrics_capturemean) sum_rg_metrics_capturemean from prepT pt, Lane l where pt.prepID = l.prepID and pt.experiment_id = {0} group by pt.prepID;".format(exp_id)
               cur.execute(query)
               #print("exp ID {0} wanted {1}, got {2}".format(exp_id, wes_min_cov,cur.fetchone()["sum_rg_metrics_capturemean"])) 

               query = "select count(*) is_sample_qualified from (select pt.experiment_id,sum(rg_metrics_capturemean) sum_rg_metrics_capturemean from prepT pt, Lane l where pt.prepID = l.prepID and pt.experiment_id = {0} group by pt.experiment_id having sum_rg_metrics_capturemean > {1}) tmp;".format(exp_id, wes_min_cov)
               cur.execute(query)
               if cur.fetchone()["is_sample_qualified"] == 1:
                  #print("{} qualifies".format(exp_id))
                  qualifying += 1
        
            if (float(qualifying/total) > 0.5):
                print("Pool {0} qualifies as {1} of {2} samples have reached coverage so updating the pool".format(poolID, qualifying, total))
                query = "select * from pool where id = {0} and is_releasable = 0 and sample_type = 'Exome';".format(poolID)
                cur.execute(query)
                #print("affected = {0}".format(cur.rowcount))
                update_isreleasable = "update pool set is_releasable = 1 where id = {0} and is_releasable = 0 and sample_type = 'Exome';".format(poolID)
                print(update_isreleasable)
                cur.execute(update_isreleasable) 
            elif Email ==1:
                rejected_pools.append(pool_name)
                print("Pool {0},{1} cant be released as only {2} of {3} samples have reached coverage".format(poolID, pool_name, qualifying, total))

if Email ==1 and len(rejected_pools) > 0: 
   body = "The following pools could not be released:\n"
   body += ", ".join(str(pool) for pool in rejected_pools)
   emailit("Pool rejection", body) 

