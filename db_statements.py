GET_QUALIFIED_BAMS = """
    SELECT fcillumid,laneNum FROM prepT p
    JOIN Lane l on p.PREPID=l.PREPID
    JOIN Flowcell f on l.fcid=f.fcid
    WHERE P_PREPID = {pseudo_prepid} AND
    FCILLUMID NOT LIKE 'X%'
    """
    # AND RELEASED = 1
    # Add when release module goes live for HTS

GET_MACHINE_FROM_FCILLUMID_QUERY = """
    SELECT MACHINE
    FROM Flowcell
    WHERE FCILLUMID = "{fcillumid}"
    """

GET_PREPID_FROM_P_PREPID = """
    SELECT PREPID
    FROM prepT
    WHERE P_PREPID = {}
    AND FAILEDPREP = 0
    """
GET_P_PREPID_FROM_PREPID = """
    SELECT DISTINCT P_PREPID
    FROM prepT 
    WHERE PREPID = {prepid}
    """

UPDATE_P_PREPID_IN_PREPT = """
    UPDATE prepT
    SET P_PREPID = "{ppid}"
    WHERE prepID = "{pid}"
    """

GET_YIELD_FROM_PPID = """
    SELECT SUM(LNYIELD) AS LANE_YIELD_SUM
    FROM prepT p
    JOIN Lane l ON p.prepid=l.prepid
    WHERE p.p_prepid={ppid}
    """

GET_MACHINE_FAILED_STATUS_FROM_FCILLUMID_QUERY = """
    SELECT FAIL
    FROM Flowcell
    WHERE FCILLUMID = "{fcillumid}"
    """

GET_MACHINE_COMPLETE_STATUS_FROM_FCILLUMID_QUERY = """
    SELECT COMPLETE
    FROM Flowcell
    WHERE FCILLUMID = "{fcillumid}"
    """

GET_FLOWCELL_RECIPE = """
    SELECT LenR1,LenI1,LenI2,LenR2
    FROM Flowcell
    WHERE FCILLUMID = "{fcillumid}"
    """

GET_POOLID_FROM_DBID = """
    SELECT CHGVID
    FROM SampleT
    WHERE DBID = "{DBID}"
    """

GET_FLOWCELL_PROJECTS = """
    SELECT distinct(s.GAFBIN)
    FROM SampleT s
    JOIN Lane l ON s.DBID=l.DBID
    JOIN Flowcell f ON f.FCID=l.FCID
    WHERE FCILLUMID = "{fcillumid}"
    """

GET_GAFBIN_FROM_SAMPLE_NAME = """
    SELECT GAFBIN
    FROM SampleT
    WHERE CHGVID = "{CHGVID}"
    """
GET_TOTAL_NUM_LANES_FROM_FLOWCELL = """
    SELECT COUNT(DISTINCT lanenum) as TOTAL_NUM_LANES
    FROM Lane l
    JOIN Flowcell f on l.FCID=f.FCID
    WHERE FCILLUMID = "{fcillumid}"
    """
GET_CLUSTER_DENSITY_FOR_LANE = """
    SELECT DISTINCT CLUSTDEN
    FROM Lane l
    JOIN Flowcell f on l.FCID=f.FCID
    WHERE FCILLUMID="{fcillumid}" AND LANENUM={lanenum}
    """

GET_SEQTYPE_FROM_PREPID = """
    SELECT SEQTYPE
    FROM SeqType
    WHERE PREPID={prepid}
    """

GET_FLOWCELL_CREATOR = """
    SELECT REPLACE(name,' ','') as FULLNAME
    FROM Flowcell f
    JOIN users u ON f.userid=u.userid
    WHERE FCILLUMID = "{fcillumid}"
    """

GET_USERID_FROM_UNI = """
    SELECT USERID
    FROM users
    WHERE NETID="{uni}"
    """

GET_FLOWCELL_CHEMVER = """
    SELECT CHEMVER
    FROM Flowcell f
    WHERE FCILLUMID = "{fcillumid}"
    """

GET_PREPID_ON_LANE = """
    SELECT PREPID
    FROM Lane l
    JOIN Flowcell f
    ON l.FCID=f.FCID
    WHERE FCILLUMID="{fcillumid}" AND l.laneNum="{lane_num}"
    ORDER BY DBID
    """

GET_DBID_ON_LANE = """
    SELECT DBID
    FROM Lane l
    JOIN Flowcell f
    ON l.FCID=f.FCID
    WHERE FCILLUMID="{fcillumid}" AND l.laneNum="{lane_num}"
    ORDER BY DBID
    """
