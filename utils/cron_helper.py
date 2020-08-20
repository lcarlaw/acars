from datetime import datetime

def timestamp(string):
    """Print out some simple date and time information for log files 
    """
    print("%s  %s : " % (string, datetime.strftime(datetime.now(), "%c")), end="")
