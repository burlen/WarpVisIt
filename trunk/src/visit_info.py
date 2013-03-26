#visitdir = "/global/common/hopper2/graphics/visit/2.5.2/linux-x86_64"
#Get the VisIt environment
import subprocess
visitenv = subprocess.check_output("visit -env" , shell=True)
#Get the visit director from the environment
visitdir = [ item for item in visitenv.split("\n") if item.startswith("VISITARCHHOME") ][0]
visitdir = visitdir.strip("VISITARCHHOME=")
print visitdir
#Define the other paths we need
visitdir_bin = visitdir + "/bin"
visitdir_lib = visitdir +"/lib"
visitdir_sitepackages = visitdir_lib + "/site-packages"

