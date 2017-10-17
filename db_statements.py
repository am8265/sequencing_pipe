GET_MACHINE_FROM_FCILLUMID_QUERY = """
    SELECT MACHINE
    FROM Flowcell
    WHERE FCIllumID = "{fcillumid}"
    """

GET_MACHINE_FAILED_STATUS_FROM_FCILLUMID_QUERY = """
    SELECT FAIL
    FROM Flowcell
    WHERE FCIllumID = "{fcillumid}"
    """

GET_MACHINE_COMPLETE_STATUS_FROM_FCILLUMID_QUERY = """
    SELECT COMPLETE
    FROM Flowcell
    WHERE FCIllumID = "{fcillumid}"
    """

GET_FLOWCELL_RECIPE = """
    SELECT LenR1,LenI1,LenI2,LenR2
    FROM Flowcell
    WHERE FCILLUMID = "{fcillumid}"
    """

GET_FLOWCELL_PROJECTS = """
    SELECT distinct(s.GAFBIN)
    FROM SampleT s
    JOIN Lane l ON s.DBID=l.DBID
    JOIN Flowcell f ON f.FCID=l.FCID
    WHERE FCILLUMID = "{fcillumid}"
    """

GET_TOTAL_NUM_LANES_FROM_FLOWCELL = """
    SELECT COUNT(DISTINCT lanenum) as TOTAL_NUM_LANES
    FROM Lane l
    JOIN Flowcell f on l.FCID=f.FCID
    WHERE FCILLUMID = "{fcillumid}"
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
