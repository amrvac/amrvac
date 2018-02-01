import sys,os
#this python script returns the list and the information on the functions
# in the krome_user.f90 file or krome_all.f90
# the details are taken from the block of comments before the function when present
# Grassi+KROME Team 01Jul2015, last edit 27Aug2016

#NO NEED TO MODIFY THIS SCRIPT
listNamesOnly = False
argv = sys.argv
searchAll = searchName = searchDescription = ""
if(len(argv)>1):
	#list of arguments available
	args = ["-h","-s=","-sf=","-sd"]
	#check if arguments are in the list
	for v in argv[1:]:
		aFound = False
		for a in args:
			#if argument found go to next
			if(a in v):
				aFound = True
				break
		#if never found error message add -h to show help
		if(not(aFound)):
			print "argument "+v+" unknown!"
			argv.append("-h")
			break
	#show help
	if("-h" in argv):
		print "usage: python "+argv[0]+" [OPTIONS]"
		print "\t-n list names only"
		print "\t-sf=STRING search string in function name"
		print "\t-sd=STRING search string in description"
		print "\t-s=STRING search string in both"
		print "\t-h show this help and exit"
		sys.exit()
	#store search string to search later on
	for v in argv:
		if("-s=" in v): searchAll = v.replace("-s=","").strip().lower()
		if("-sf=" in v): searchName = v.replace("-sf=","").strip().lower()
		if("-sd=" in v): searchDescription = v.replace("-sd=","").strip().lower()
	listNamesOnly = ("-n" in argv)

fileName = "krome_user.f90"
if(not(os.path.isfile(fileName))): fileName = "krome_all.f90"

fh = open(fileName,"rb")
alltext = ""
for row in fh:
	alltext += row

alltext = alltext.replace("&\n","")

inpartF = ["function","(",")"]
inpartS = ["subroutine","(",")"]
infun = issub = False
readData = False
storecom = ""
flist = []
for row in alltext.split("\n"):
	srow = row.strip()
	if(srow==""): continue
	if("module krome_user" in srow): readData = True
	if("end module krome_user" in srow): break
	if(not(readData)): continue
	if(srow=="contains"): storecom = ""
	isfun = True
	for part in inpartF:
		if(not(part in srow) or srow[0]=="!"): 
			isfun = False
			break
	issub = True
	for part in inpartS:
		if(not(part in srow) or srow[0]=="!"):
			issub = False
			break

	if(isfun or issub):
		if(storecom==""): storecom = "[no comments available]"
		ffname = srow
		#if search options are available apply
		okList = False
		#search in function name
		if(searchName in ffname.lower() and searchName!=""): okList = True
		#search in function description
		if(searchDescription in storecom.lower() and searchDescription!=""): okList = True
		#search in name+description
		if((searchAll in storecom.lower() or searchAll in ffname.lower()) and searchAll!=""): okList = True
		#if search strings are all empty skip
		if(searchAll=="" and searchName=="" and searchDescription==""): okList = True
		#append data to the list
		flist.append([ffname,storecom.strip(),okList])
		continue

	if(("end function" in srow) or ("end subroutine" in srow)):
		storecom = ""

	if(srow[0]=="!" and not("*******" in srow)): storecom += "  "+srow[1:].strip()+"\n"

#sort aplhabetically
flist = sorted(flist,key=lambda x:x[0])
#print resuslts
icount = 0
for x in flist:
	#if not matches searh criteria skip
	if(not(x[2])): continue
	fname = x[0]
	while("  " in fname):
		fname = fname.replace("  "," ")
	print str(icount+1)+") "+fname.replace(", ",",")
	print "  "+x[1]+"\n"
	icount += 1

