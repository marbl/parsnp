import os,sys,string
#parsnp basic INSTALL script
user_home = os.environ["HOME"]
print "<<Welcome to Parsnp utility script install>>"

#check for python version
if (sys.version_info[0] < 2) or (sys.version_info[0] == 2 and sys.version_info[1] < 6):
    
    print "Python version is %s. Parsnp requires at least 2.6"%(sys.version)
    sys.exit(1)

#complete shebang
os.system("python setup.py install_scripts --install-dir=`pwd` build")

scripts = ["parsnp.py"]
#copy to currdir
files = os.listdir(".")
for script in scripts:
    os.system("mv %s %s"%(script,script.replace(".py","")))