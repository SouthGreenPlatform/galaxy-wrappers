from subprocess import check_output
import sys, os, datetime, re

def updateTag() :
    if not os.path.isdir(".git"):
        sys.exit("This does not appear to be a Git repository.")
    out=check_output(["git", "log", "--tags", "--simplify-by-decoration", "--pretty='format=:%aD %d'"]).decode("utf-8")
    vers=""
    outs=out.split("\n")
    with open(".version","w") as verOut :
        print("File '.version' created.")
        for i in outs :
            date=re.search(":\w+,\s(\d+\s\w+\s\d+)\s",i)
            tag=re.search("tag:\s\w(\d+\.\d+\.\d+)",i)
            if date and tag :
                verOut.write(tag.group(1) + " - " +date.group(1)+"\n")

def getVersion() :
    with open(".version","r") as verF :
        f=verF.readlines()
        print("\n"+"Version : "+f[0].split(" - ")[0])
        print("Last update : "+f[0].split(" - ")[1])
    return

