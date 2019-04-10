#sync data back to this system
import sys
import subprocess

outfolder=sys.argv[1]
p=subprocess.call(["rsync","-arhP","subhadra@192.168.1.244:/storage/subhadra/kabir/output/freq/"+outfolder+"/","/data/kabir/output/freq/"+outfolder+"/"])
