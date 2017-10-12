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

