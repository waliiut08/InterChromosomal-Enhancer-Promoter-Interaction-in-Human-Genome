import os
from datetime import datetime, timedelta 

#----------------------------------------------------------------------------------------------
def read_file(filename):
	enhancerFL = open(filename,"r")
	enhancerData = enhancerFL.readlines()
	enhancerFL.close()
	return enhancerData

#----------------------------------------------------------------------------------------------
def write_file(data,filename,mode):
	overlappedData = open(filename,mode)
	overlappedData.write(data)
	overlappedData.close()


#-----------------------------------------------------------------------------------------------
def overlapped(start1,end1,start2,end2):
	return int(end1)>=int(start2) and int(end2) >= int(start1)


#-----------------------------------------------------------------------------------------------


def read_data(fileRawData):
	return [item.replace("\n","").split("\t") for item in fileRawData]
#-----------------------------------------------------------------------------------------------


## find between strings............
def find_between( string, first, last ):
    try:
        start = string.index( first ) + len( first )
        end   = string.rindex( last)
        return string[start:end]
    except ValueError:
        return "no such string found.."

#------------------------------------------------------------------------------------------------


def make_directory(outputFile):
	if not os.path.exists(os.path.dirname(outputFile)):
		try:
			os.makedirs(os.path.dirname(outputFile))
		except OSError as exc:
			if exc.errno != errno.EEXIST:
				raise 

#------------------------------------------------------------------------------------------------

def elapsed_time(start, end):
	elapsedTime = end - start
	m, s = divmod(elapsedTime, 60)
	h, m = divmod(m, 60)
	print "Running time: %02d:%02d:%02d" % (h, m, s)

#------------------------------------------------------------------------------------------------

def expected_finishing_time(totalRequiredTime):
	will_finish_at = datetime.now() + timedelta(minutes = totalRequiredTime)
	print "Expected finishing time:- " + str(format(will_finish_at, '%H:%M:%S')) + "\n"
	
